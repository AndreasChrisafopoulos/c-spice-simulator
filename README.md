# C SPICE Simulator

*A SPICE-like circuit simulator written from scratch in C — Modified Nodal Analysis, direct and iterative solvers (dense and sparse), DC / DC-sweep / transient analysis, exercised against real IBM power-delivery-network benchmarks.*

## Overview

This project implements a SPICE-like circuit simulator in C, built around Modified Nodal Analysis (MNA). It parses SPICE-style netlists and solves DC operating point, DC sweep, and transient analyses. Every analysis is available through two parallel numerical paths: a dense path (GSL-backed, and a hand-written from-scratch LU / Cholesky / CG / BiCG implementation) and a sparse path for large circuits. The simulator is exercised against real IBM power-delivery-network benchmark netlists ranging from tens of thousands to several million circuit elements.

## Features

- **Netlist Parsing** — SPICE-like, case-insensitive syntax parser with hash-table-based node/element lookup.
- **Multiple Analyses** — DC operating point (`.op`), DC sweep (`.dc`, up to 32 sweeps per netlist), and transient analysis (`.tran`).
- **Dual Numerical Paths** — every solve is available through GNU Scientific Library (GSL) routines *and* a from-scratch implementation (LU with partial pivoting, Cholesky, hand-written forward/backward substitution), selectable with `.options custom`.
- **Iterative Solvers** — Jacobi-preconditioned Conjugate Gradient (CG, for SPD systems) and BiConjugate Gradient (BiCG), hand-implemented for both dense and sparse matrices.
- **Sparse Support** — triplet/CSC sparse matrices for large-scale networks (tested up to millions of elements, e.g. the IBM power-grid benchmarks).
- **Direct vs. Iterative Comparison** — every run writes two parallel output sets so the solver families can be compared directly (see Output Files).

## Supported Circuit Elements

| Element | Syntax | Status |
|---|---|---|
| Resistor | `R<name> <n+> <n-> <value>` | Fully simulated |
| Capacitor | `C<name> <n+> <n-> <value>` | Fully simulated (transient) |
| Inductor | `L<name> <n+> <n-> <value>` | Fully simulated |
| Voltage Source | `V<name> <n+> <n-> <dc> [tran_spec]` | Fully simulated (DC + transient) |
| Current Source | `I<name> <n+> <n-> <dc> [tran_spec]` | Fully simulated (DC + transient) |
| Diode | `D<name> <n+> <n-> <model> [area=<v>]` | Parsed only — not yet stamped |
| MOSFET | `M<name> <d> <g> <s> <b> <model> L=<v> W=<v>` | Parsed only — not yet stamped |
| BJT | `Q<name> <c> <b> <e> <model> [area=<v>]` | Parsed only — not yet stamped |

Diode, MOSFET, and BJT lines are fully parsed — nodes, model name, and `area`/`L`/`W` are read and stored — but are not yet stamped into the MNA matrix, so they don't currently affect simulation results. `.MODEL` parameters are likewise parsed but not used. Extending the stamping stage with device models is the natural next step; see Known Limitations.

Transient waveforms available on V/I sources: `EXP`, `SIN`, `PULSE`, `PWL` (up to 20 time/value pairs). Full syntax in [`NETLIST.md`](./NETLIST.md).

## Supported Analyses

- `.op` — DC operating point
- `.dc <source> <start> <end> <step>` — DC sweep (up to 32 per netlist)
- `.tran <time_step> <final_time>` — transient analysis, trapezoidal (default) or backward Euler (`.options method=be`)
- `.plot` / `.print V(n1) V(n2) ...` — up to 32 plotted nodes, context-sensitive to the preceding `.dc`/`.tran` block

## Numerical Methods

**Direct**
- Dense LU (partial pivoting) and Cholesky — via GSL, or a from-scratch implementation with hand-written forward/backward substitution (`.options custom`)
- Sparse LU and Cholesky — via a compressed-column sparse factorization kernel (`.options sparse`)

**Iterative** (`.options iter`)
- Conjugate Gradient (CG) for SPD systems (`.options spd`), BiConjugate Gradient (BiCG) otherwise
- Jacobi (diagonal) preconditioning; converges when the relative residual norm drops below `itol` (default `1e-3`, set with `.options itol=<value>`)
- Implemented for both dense (GSL BLAS-based) and sparse (CSC-based) matrices

## Third-Party Components

Two well-known libraries are bundled to keep the simulator dependency-light; neither is implemented by this project:

- **[uthash](https://troydhanson.github.io/uthash/)** by Troy D. Hanson (`include/uthash.h`) — hash table used for node-name lookup during parsing.
- **CSparse** by Timothy A. Davis (`csparse/`) — sparse LU/Cholesky factorization kernel used by the sparse solver path.

Everything else — the netlist parser, MNA stamping, dense/sparse direct and iterative solvers, transient integration, and I/O — is implemented from scratch.

## Project Architecture

```
       [Netlist File (.cir / .spice)]
                    │
                    ▼
               [ Parser ]  <── Checks syntax, creates data structures
                    │
                    ▼
     [ Modified Nodal Analysis ]
                    │
                    ▼
            [ Matrix Assembly ]
                    │
       ┌────────────┼────────────┐
       ▼            ▼            ▼
    [Direct]   [Iterative]  [Sparse] <── Solvers
       │            │            │
       └────────────┼────────────┘
                    ▼
          [ Transient Analysis ] <── Time Integration (TR/BE)
                    │
                    ▼
             [ Output Files ]
```

## Dependencies

- GCC with C11 support
- [GSL – GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [OpenBLAS](https://www.openblas.net/)

On Ubuntu/Debian:

```bash
sudo apt-get install libgsl-dev libopenblas-dev
```

## Building

```bash
make
```

Produces the `project` executable, linked against GSL and OpenBLAS (see `Makefile`).

## Running

```bash
./project <part_number> <netlist_file>
```

- `part_number` selects the netlist folder: `1` → `Part1_Netlists/`, `3` → `Part3_Netlists/`, `6` → `Part6_Netlists/`
- `netlist_file` is the file name inside that folder

Example:

```bash
./project 6 ibmpg1t.spice
```

## Output Files

Every run writes results under `OUT/<netlist_name>_outputfiles0/` (direct solver) and `OUT/<netlist_name>_outputfiles1/` (iterative solver, populated when `.options iter` is set) — so the two solver families can be compared directly on the same circuit.

## Validating Results

`cmpr/compare_results.py` compares simulator output against a reference solution and reports mean absolute/relative error and pass rate. See [`cmpr/README.md`](./cmpr/README.md) for usage. Reference solutions for the IBM power-grid benchmarks are in `IBM_SOLS/`.

## Repository Layout

| Path | Contents |
|---|---|
| `src/`, `include/` | Simulator source — parser, MNA stamping, solvers, transient analysis |
| `csparse/` | Sparse-matrix factorization library (third-party, see above) |
| `Part1_Netlists/` | Test netlists covering basic element types and DC analysis |
| `Part3_Netlists/` | Test netlists for sparse + iterative solvers, including the IBM power-grid DC benchmarks (`ibmpg1`–`ibmpg6`) |
| `Part6_Netlists/` | Test netlists for transient analysis, including the IBM power-grid transient benchmarks |
| `IBM_SOLS/` | Reference solutions for the IBM benchmark netlists |
| `cmpr/` | Output comparison / validation script |
| `NETLIST.md` | Full netlist syntax specification |

## Known Limitations

- Diode, MOSFET, and BJT elements are parsed but not yet stamped into the MNA system (see Supported Circuit Elements); `.MODEL` parameters are read but not used.
- Controlled sources (VCVS, VCCS, CCVS, CCCS) are not supported.
- Maximum netlist line length is 256 characters; node/element/model names up to 32 characters.
- PWL sources support up to 20 (time, value) pairs.

Full details in [`NETLIST.md`](./NETLIST.md).


*Developed as a university project, October 2025 – January 2026.*