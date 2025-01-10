#include "sssp.h"

#include <fstream>
#include <functional>
#include <numeric>

#include "dijkstra.h"

using namespace std;
using namespace parlay;

template <typename NodeID_, typename rng_t_,
          typename uNodeID_ = typename std::make_unsigned<NodeID_>::type>
class UniDist {
public:
  UniDist(NodeID_ max_value, rng_t_ &rng) : rng_(rng) {
    no_mod_ = rng_.max() == static_cast<uNodeID_>(max_value);
    mod_ = max_value + 1;
    uNodeID_ remainder_sub_1 = rng_.max() % mod_;
    if (remainder_sub_1 == mod_ - 1)
      cutoff_ = 0;
    else
      cutoff_ = rng_.max() - remainder_sub_1;
  }

  NodeID_ operator()() {
    uNodeID_ rand_num = rng_();
    if (no_mod_)
      return rand_num;
    if (cutoff_ != 0) {
      while (rand_num >= cutoff_)
        rand_num = rng_();
    }
    return rand_num % mod_;
  }

private:
  rng_t_ &rng_;
  bool no_mod_;
  uNodeID_ mod_;
  uNodeID_ cutoff_;
};

template <class Algo>
void run(Algo &algo, const Graph &G, int rounds, int sources, bool verify) {
  std::mt19937_64 rng(27491095);
  UniDist<NodeId, std::mt19937_64> udist(G.n - 1, rng);

  double total_time = 0;
  for (int v = 0; v < sources; v++) {
    NodeId s, deg;
    do {
      s = udist();
      deg = G.offset[s + 1] - G.offset[s];
    } while (deg == 0);

    printf("source %d: %-10d\n", v, s);
    if (v == 0) {
      internal::timer t;
      algo.sssp(s);
      t.stop();
      printf("Warmup Round: %f\n", t.total_time());
      fflush(stdout);
    }

    for (int i = 0; i < rounds; i++) {
      internal::timer t;
      algo.sssp(s);
      t.stop();
      printf("Round %d: %f\n", i, t.total_time());
      fflush(stdout);
      total_time += t.total_time();
    }

    // ofstream ofs("sssp.tsv", ios_base::app);
    // ofs << average_time << '\t';
    // ofs.close();

    if (verify) {
      printf("Running verifier...\n");
      internal::timer t;
      auto dist = algo.sssp(s);
      t.stop();
      printf("Our running time: %f\n", t.total_time());
      verifier(s, G, dist);
    }
    printf("\n");
  }
  double average_time = total_time / rounds;
  printf("Average time: %f\n", average_time);
}

int main(int argc, char *argv[]) {

  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-p parameter] [-w] [-s] [-v] [-a "
            "algorithm]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-p,\tparameter(e.g. delta, rho)\n"
            "\t-w,\tweighted input graph\n"
            "\t-s,\tsymmetrized input graph\n"
            "\t-v,\tverify result\n"
            "\t-a,\talgorithm: [rho-stepping] [delta-stepping] [bellman-ford]\n"
            "\t-S,\tnumber of sources\n"
            "\t-n,\tnumber of trials per source\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }

  char c;
  bool weighted = false;
  bool symmetrized = false;
  bool verify = false;
  string param;
  int algo = rho_stepping;
  char const *FILEPATH = nullptr;
  int rounds = 22;
  int sources = 1;

  while ((c = getopt(argc, argv, "i:p:a:wsvn:S:")) != -1) {
    switch (c) {
    case 'i':
      FILEPATH = optarg;
      break;
    case 'p':
      param = string(optarg);
      break;
    case 'a':
      if (!strcmp(optarg, "rho-stepping")) {
        algo = rho_stepping;
      } else if (!strcmp(optarg, "delta-stepping")) {
        algo = delta_stepping;
      } else if (!strcmp(optarg, "bellman-ford")) {
        algo = bellman_ford;
      } else {
        fprintf(stderr, "Error: Unknown algorithm %s\n", optarg);
        exit(EXIT_FAILURE);
      }
      break;
    case 'w':
      weighted = true;
      break;
    case 's':
      symmetrized = true;
      break;
    case 'v':
      verify = true;
      break;
    case 'n':
      rounds = atoi(optarg);
      break;
    case 'S':
      sources = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Error: Unknown option %c\n", optopt);
      exit(EXIT_FAILURE);
    }
  }
  Graph G(weighted, symmetrized);

  printf("Reading graph...\n");
  G.read_graph(FILEPATH);
  if (!weighted) {
    printf("Generating edge weights...\n");
    G.generate_weight();
  }

  fprintf(stdout,
          "Running on %s: |V|=%zu, |E|=%zu, param=%s, num_src=%d, "
          "num_round=%d\n",
          FILEPATH, G.n, G.m, param.c_str(), sources, rounds);

  int sd_scale = G.m / G.n;
  if (algo == rho_stepping) {
    size_t rho = 1 << 20;
    if (param != "") {
      rho = stoull(param);
    }
    Rho_Stepping solver(G, rho);
    solver.set_sd_scale(sd_scale);
    run(solver, G, rounds, sources, verify);
  } else if (algo == delta_stepping) {
    EdgeTy delta = 1 << 15;
    if (param != "") {
      if constexpr (is_integral_v<EdgeTy>) {
        delta = stoull(param);
      } else {
        delta = stod(param);
      }
    }
    Delta_Stepping solver(G, delta);
    solver.set_sd_scale(sd_scale);
    run(solver, G, rounds, sources, verify);
  } else if (algo == bellman_ford) {
    Bellman_Ford solver(G);
    solver.set_sd_scale(sd_scale);
    run(solver, G, rounds, sources, verify);
  }
  return 0;
}
