#include <cerrno>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include "kplex_bits.h"

static void usage(const char *program) {
    std::printf("Usage: %s [-g <dataset-name-or-directory>] -k <k> -q <tau> [-b]\n"
                "  -g       Dataset name under ./dataset (default: jazz), or a directory path\n"
                "  -k       Paper k, including the vertex itself (1 <= k <= 32)\n"
                "  -q       Result size tau; --tau is an alias (2*k-1 <= tau <= 64)\n"
                "  -b       Read b_degree.bin and b_adj.bin instead of edges.txt\n"
                "  -h       Show this help\n", program);
}

static int integer(const std::string &value) {
    char *end = nullptr;
    errno = 0;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if(errno || end == value.c_str() || *end || parsed < 0 || parsed > INT_MAX)
        throw std::invalid_argument("invalid integer: " + value);
    return static_cast<int>(parsed);
}

int main(int argc, char **argv) {
    try {
        std::string dataset = "jazz";
        int paper_k = -1, tau = -1;
        bool binary = false;
        for(int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if(option == "-h" || option == "--help") { usage(argv[0]); return 0; }
            if(option == "-b") { binary = true; continue; }
            if(option != "-g" && option != "-k" && option != "-q" && option != "--tau")
                throw std::invalid_argument("unknown argument: " + option);
            if(i + 1 == argc) throw std::invalid_argument("missing value for " + option);
            const std::string value = argv[++i];
            if(option == "-g") dataset = value;
            else if(option == "-k") paper_k = integer(value);
            else tau = integer(value);
        }
        if(dataset.empty() || paper_k < 1 || paper_k > 32 || tau < 1 || tau > 64 || tau < 2*paper_k-1)
            throw std::invalid_argument("requires 1 <= k <= 32 and 2*k-1 <= tau <= 64");
        if(!__builtin_cpu_supports("popcnt") || !__builtin_cpu_supports("bmi2"))
            throw std::runtime_error("this implementation requires an x86-64 CPU with POPCNT and BMI2");

        const std::string directory = dataset.find('/') == std::string::npos
            ? "dataset/" + dataset : dataset;
        Graph graph;
        if(binary) graph.read_binary(directory); else graph.read_text(directory);
        std::printf("Graph: %s, n: %u, m: %u (undirected)\n",
                    directory.c_str(), graph.n, graph.pstart[graph.n]/2);

        const auto start = std::chrono::steady_clock::now();
        KPlexBits solver(graph, paper_k-1, tau);
        solver.release_input = &graph;
        solver.run();
        const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now()-start).count();
        std::printf("Fixed-size K-Plex Count: %llu, tau: %d, k: %d, Total Time: %lld (microseconds)\n",
                    solver.answer, tau, paper_k, static_cast<long long>(elapsed));
        return 0;
    } catch(const std::exception &error) {
        std::fprintf(stderr, "Error: %s\n", error.what());
        return 1;
    }
}
