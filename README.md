# Parallel — DMDGP (CPU baseline) with OpenMP search over sign masks

This repository is a modular C implementation of a **DMDGP-like** pipeline where each configuration is encoded by a bitmask `k` controlling the signs of the dihedral sine terms. The goal is to efficiently search over all masks and find **one feasible embedding** (any solution) using **CPU parallelism (OpenMP)**, with a design that is intentionally GPU-friendly for future work.

## Problem summary (high level)

Given a DMDGP instance with vertices in the fixed order `1..n`:

- Distances `d(u,v)` are provided for a set of edges.
- Distances and the vertex order define:
  - bond lengths `d(k-1,k)`
  - bond angles `theta[k]` for `k = 3..n`
  - dihedral cosines `cos(omega[k])` for `k = 4..n`
  - absolute dihedral sines `|sin(omega[k])|` for `k = 4..n`

For each bitmask `k` with `m = n-3` bits:

- Bit `(n - t)` determines the **sign of** `sin(omega[t])` for `t = 4..n`.
- `h(k)` builds an embedding (points in `R^3`) using a 4×4 homogeneous transform chain.
- `g(h(k))` scores constraint violations using **no square roots**:
  \[
    g(x) = \sum_{(u,v,d)\in E} (\|x_u-x_v\|^2 - d^2)^2
  \]
- We search for any mask such that `g(h(k)) <= delta`.

This repository currently returns the **first found** solution (not necessarily the smallest `k`), optimized for speed.

---

## Input format

The input file is:

```txt

n m
u1 v1 d1
u2 v2 d2
...
um vm dm

```

Example:

```txt

7 18
5 6 1.5260000000
2 6 3.9491605627
...

````

### DMDGP validation assumptions

We assume the instance is a valid DMDGP in the order `1..n`. The program will fail early (with a clear message) if required distances are missing, in particular the backbone distances needed to compute:

- `theta[k]` from `{d(k-2,k-1), d(k-1,k), d(k-2,k)}`
- `cos(omega[k])` from `{d(k-3,k-2), d(k-2,k-1), d(k-1,k), d(k-2,k), d(k-3,k)}`

---

## Build

Requirements:

- `cc` (clang or gcc)
- OpenMP runtime (`-fopenmp`)
- `libm`

Build:

```bash
make
````

Binaries are emitted into `./build/`.

---

## Executables

### 1) `precompute`

Validates the instance and prints precomputed arrays (`theta`, `cos(omega)`, `|sin(omega)|`).

```bash
./build/precompute data/7_18.in
```

### 2) `points`

Computes and prints the embedding `h(k)` (points `1..n`) for a given decimal mask `k`.

```bash
./build/points data/7_18.in 0
./build/points data/7_18.in 1
```

Output format:

```txt
i  xi  yi  zi
```

### 3) `search`

Searches masks `k in [0, 2^(n-3))` in parallel and returns the first mask such that `g(h(k)) <= delta`.

```bash
OMP_NUM_THREADS=10 ./build/search data/30_168.in 1e-3
```

Recommended CPU pinning (often improves stability/perf):

```bash
OMP_PROC_BIND=true OMP_PLACES=cores OMP_NUM_THREADS=10 ./build/search data/30_168.in 1e-3
```

---

## Performance notes

- The loop is parallelized over `k` (the `2^(n-3)` configurations).
- The expensive part is building the transform chain and scoring `g`.
- Several values are precomputed to reduce per-`k` cost:

  - `bond[k] = d(k-1,k)`
  - `cos(theta[k])`, `sin(theta[k])`
  - `cos(omega[k])` and `|sin(omega[k])|`
  - (optionally) per-`k` transform matrices `A_plus[t]` / `A_minus[t]` for `t>=4`

Because the search returns the first found solution, runtime can vary depending on scheduling and when the solution’s chunk is evaluated.

---

## Project layout

```txt
include/
  instance.h     # parsing, validation, precompute (theta, omega, etc.)
  mat4.h         # 4x4 homogeneous matrix ops
  geom.h         # h(k): build points via transform chain
  score.h        # g(x): score embedding
  search.h       # OpenMP search API
src/
  instance.c
  mat4.c
  geom_mat4.c
  score.c
  search_omp.c
  precompute_main.c
  points_main.c
  search_main.c
data/
  *.in           # instances
build/
  precompute
  points
  search
```

---

## Notes for future GPU port

The code is structured so that the core pipeline:

1. `geom_build_points_mat4(I, k, x)`
2. `score_g_no_sqrt(I, x)`

is cleanly separable from CPU scheduling and can be mapped to a GPU kernel where each thread handles one `k` (or a block of `k` values), followed by a reduction/compaction step to pick a feasible solution.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
