## Efficient Enumeration of $k$-Plexes with Arbitrary Size

This repository provides the implementation of **ArbKPEnum** for enumerating $k$-plexes of a specified size $\tau$.

## Code Structure

- `src/main.cpp` – Command-line arguments, graph loading, timing, and result reporting.
- `src/kplex_base.h` – Degeneracy ordering, global peeling, root subgraph construction, iterative pruning, and pool combination enumeration.
- `src/kplex_bits.h` – Bit-matrix search, relaxed vertex pool identification and updates, and constrained pool enumeration.
- `src/graph.h` – Graph input and CSR storage.
- `src/linear_heap.h` – Bucket-based minimum heap for degeneracy ordering and peeling.
- `src/cuckoo_hash.h` – Adjacency lookup during root subgraph construction.
- `dataset/` – Example datasets.

---

## Compile the codes

From the repository root, run:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j
cd ..
```

The executable `kplex` will be generated in the `build/` directory.

---

## Datasets

- **Included datasets:** `jazz` (198 vertices, 2,742 undirected edges) and `caida` (as-caida; 26,475 vertices, 53,381 undirected edges), located in `dataset/`.
- **Additional datasets** are available from [SNAP](https://snap.stanford.edu/data/index.html) and [Network Repository](https://networkrepository.com/). Convert them to one of the supported input formats before running the algorithm.

Each included dataset provides both formats:

- `edges.txt` – The first line contains `n m`, followed by `m` undirected edges `u v`. Vertex IDs range from `0` to `n-1`. Each edge needs to appear only once.
- `b_degree.bin` and `b_adj.bin` – Binary CSR input using native-endian 32-bit unsigned integers. `b_degree.bin` contains the integer width (`4`), vertex count, total adjacency-entry count (`2m`), and vertex degrees. `b_adj.bin` stores each vertex's neighbors consecutively, with both directions of every undirected edge present. The included binary files are little-endian.

---

## Run the procedure

Run the following commands from the **repository root** so that the relative path `dataset/` resolves correctly.

### Example Commands

- **Enumerate $(2,6)$-plexes on jazz:**

```bash
./build/kplex -g jazz -k 2 -q 6
```

- **Enumerate $(2,8)$-plexes on caida:**

```bash
./build/kplex -g caida -k 2 -q 8
```

- **Enumerate $(3,8)$-plexes on caida using binary input:**

```bash
./build/kplex -g caida -k 3 --tau 8 -b
```

### Parameters

| Parameter | Description |
| --- | --- |
| `-g NAME` | Read `dataset/NAME/`. Defaults to `jazz`. |
| `-g PATH` | Read a dataset directory directly when the value contains `/`, e.g., `dataset/jazz`. |
| `-k K` | Required. Paper parameter $k$, including the vertex itself; each vertex may have at most `K-1` other non-neighbors. |
| `-q TAU` / `--tau TAU` | Required. Result size $\tau$. Supported range: `1 <= K <= 32` and `2*K-1 <= TAU <= 64`. |
| `-b` | Use binary input. The default input is `edges.txt`. |
| `-h` / `--help` | Show usage information. |

The program reports the number of results and running time in microseconds. Timing includes preprocessing and enumeration, but excludes graph loading.
