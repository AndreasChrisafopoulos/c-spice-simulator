# C SPICE Simulator

A SPICE-like analog circuit simulator implemented in C using Modified Nodal Analysis, dense and sparse matrix representations, direct and iterative numerical solvers, DC analysis, DC sweep, and transient simulation.

The project combines custom numerical implementations with established numerical libraries. Dense LU and Cholesky solvers can run through either GSL or custom factorization code, iterative CG and BiCG solvers are implemented for dense and sparse matrices, and sparse direct factorization uses CSparse.

## Highlights

- Sparse MNA simulation validated on IBM power-grid benchmarks
- Tested up to 127,236 circuit nodes and 127,565 MNA unknowns
- 100% of validated non-zero DC reference voltages within 1%
- 20,020 transient reference samples validated with 100% within 1%
- Dense and sparse direct and iterative solver implementations

## Overview

The simulator parses SPICE-style netlists and constructs the corresponding Modified Nodal Analysis system.

It supports:

- DC operating-point analysis
- DC source sweeps
- transient simulation
- dense and sparse matrix representations
- direct and iterative linear-system solvers
- configurable solver selection through netlist options
- transient voltage and current source waveforms
- comparison against reference solutions for large benchmark circuits

The repository includes IBM power-delivery-network benchmark netlists and corresponding reference data that can be used to exercise the sparse solver path on circuits ranging from tens of thousands to millions of elements.

## Features

- **SPICE-style netlist parsing** with case-insensitive element and command handling
- **Modified Nodal Analysis** matrix construction
- **DC operating-point analysis**
- **DC sweep analysis** with support for multiple sweeps in one netlist
- **Transient analysis** using trapezoidal integration or backward Euler
- **Dense direct solvers** using GSL or custom LU and Cholesky implementations
- **Dense iterative solvers** using custom Conjugate Gradient and BiConjugate Gradient algorithms
- **Sparse matrix support** using triplet and compressed-column representations
- **Sparse direct solvers** using CSparse LU and Cholesky factorization
- **Sparse iterative solvers** using custom CG and BiCG implementations
- **Jacobi preconditioning** for iterative methods
- **Configurable numerical options** through `.options`
- **Transient source support** for `EXP`, `SIN`, `PULSE`, and `PWL`
- **Reference-result comparison** through the included Python validation utility
- **Large benchmark support** through the included IBM power-grid netlists

## Supported circuit elements

| Element | Syntax | Status |
|---|---|---|
| Resistor | `R<name> <n+> <n-> <value>` | Simulated |
| Capacitor | `C<name> <n+> <n-> <value>` | Simulated in transient analysis |
| Inductor | `L<name> <n+> <n-> <value>` | Simulated |
| Voltage source | `V<name> <n+> <n-> <dc> [tran_spec]` | Simulated in DC and transient analyses |
| Current source | `I<name> <n+> <n-> <dc> [tran_spec]` | Simulated in DC and transient analyses |
| Diode | `D<name> <n+> <n-> <model> [area=<v>]` | Parsed but not stamped |
| MOSFET | `M<name> <d> <g> <s> <b> <model> L=<v> W=<v>` | Parsed but not stamped |
| BJT | `Q<name> <c> <b> <e> <model> [area=<v>]` | Parsed but not stamped |

Diode, MOSFET, and BJT descriptions are parsed and stored, including their nodes and model information, but these devices are not currently stamped into the MNA system and therefore do not affect simulation results.

`.MODEL` statements are also parsed but their model parameters are not currently used during simulation.

Adding nonlinear device models and the corresponding nonlinear solution flow is a natural future extension.

## Transient sources

Independent voltage and current sources can use the following transient waveforms:

- `EXP`
- `SIN`
- `PULSE`
- `PWL`

PWL sources support up to 20 time-value pairs.

The complete syntax is documented in [`NETLIST.md`](./NETLIST.md).

## Supported analyses

### DC operating point

```text
.op
```

Constructs and solves the DC MNA system for the circuit.

### DC sweep

```text
.dc <source> <start> <end> <step>
```

Sweeps an independent source over the requested range and solves the circuit at each point.

Up to 32 DC sweeps can be stored from one netlist.

### Transient analysis

```text
.tran <time_step> <final_time>
```

Runs time-domain simulation using:

- trapezoidal integration by default
- backward Euler with `.options method=be`

### Output selection

```text
.plot V(n1) V(n2)
.print V(n1) V(n2)
```

Up to 32 nodes can be selected for output.

## Numerical methods

The simulator provides several numerical paths so the same MNA formulation can be solved using different algorithms and matrix representations.

### Dense direct solvers

The dense path supports:

- GSL LU factorization
- GSL Cholesky factorization for suitable SPD systems
- custom LU factorization with partial pivoting
- custom Cholesky factorization
- custom forward substitution
- custom backward substitution

The custom dense direct path can be selected with:

```text
.options custom
```

The custom LU and Cholesky factorization logic is implemented in this project while GSL provides the alternative library-backed implementation.

### Dense iterative solvers

The dense iterative path includes custom implementations of:

- Conjugate Gradient for SPD systems
- BiConjugate Gradient for general systems
- Jacobi diagonal preconditioning

Iterative solving is enabled with:

```text
.options iter
```

SPD behavior can be selected with:

```text
.options spd
```

The convergence tolerance can be configured with:

```text
.options itol=<value>
```

The default iterative tolerance is `1e-3`.

GSL and BLAS routines are used for supporting vector and matrix operations in the dense numerical path.

### Sparse direct solvers

Sparse matrices are stored using triplet form during assembly and compressed-column form during solution.

Sparse operation can be enabled with:

```text
.options sparse
```

Sparse LU and Cholesky factorization use the bundled **CSparse** implementation.

The project code performs the MNA assembly, sparse-path integration, solver selection, right-hand-side handling, and simulation flow around the CSparse factorization routines.

### Sparse iterative solvers

The sparse iterative path contains project implementations of:

- sparse matrix-vector multiplication
- sparse transposed matrix-vector multiplication
- Jacobi preconditioning
- Conjugate Gradient
- BiConjugate Gradient

These algorithms operate directly on CSC matrices.

## Solver architecture

```text
        SPICE-style netlist
                |
                v
             Parser
                |
                v
     Circuit data structures
                |
                v
      Modified Nodal Analysis
                |
                v
          Matrix assembly
                |
        +-------+-------+
        |               |
        v               v
      Dense           Sparse
        |               |
   +----+----+     +----+----+
   |         |     |         |
   v         v     v         v
 Direct   Iterative Direct  Iterative
   |         |       |         |
 GSL /     Custom  CSparse    Custom
 Custom    CG/BiCG  LU/Chol   CG/BiCG
   |         |       |         |
   +---------+-------+---------+
                |
                v
       DC / Sweep / Transient
                |
                v
            Output files
```

## Third-party components

The simulator uses several external numerical and utility components. These components are not claimed as original implementations of this project.

### GNU Scientific Library

GSL is used for matrix and vector data structures, library-backed dense factorizations, and supporting numerical operations.

### OpenBLAS

OpenBLAS provides optimized BLAS routines used by the numerical implementation.

### uthash

[`include/uthash.h`](include/uthash.h) provides the hash-table implementation used for efficient name-based lookup while parsing circuit nodes and elements.

uthash was developed by Troy D. Hanson and is third-party code.

### CSparse

The [`csparse/`](csparse/) directory contains sparse matrix routines derived from CSparse by Timothy A. Davis.

CSparse provides the sparse LU and Cholesky factorization kernels used by the sparse direct solver path.

The project-specific sparse simulation code around these routines, including MNA assembly, solver integration, iterative methods, and simulation control, is separate from the bundled CSparse implementation.

## Dependencies

Required build dependencies:

- GCC with C11 support
- GNU Scientific Library
- OpenBLAS
- standard C math library
- Git LFS for the large IBM benchmark and reference files

On Ubuntu or Debian:

```bash
sudo apt-get install libgsl-dev libopenblas-dev git-lfs
git lfs install
git lfs pull
```

## Building

From the repository root:

```bash
make
```

The Makefile compiles the project source files under `src/` together with the bundled CSparse source.

The resulting executable is:

```text
project
```

To remove generated build files:

```bash
make clean
```

## Running

The simulator is executed as:

```bash
./project <part_number> <netlist_file>
```

The first argument selects the corresponding netlist directory.

Examples include:

```text
1 -> Part1_Netlists/
3 -> Part3_Netlists/
6 -> Part6_Netlists/
```

Example:

```bash
./project 6 ibmpg1t.spice
```

The program constructs the path to the selected netlist, parses the simulation options, assembles the required matrices, runs the requested analyses, and writes the resulting node data to the output directory.

## Solver selection

Solver behavior is controlled by the options contained in each netlist.

Important options include:

```text
.options custom
.options sparse
.options iter
.options spd
.options itol=<value>
.options method=be
```

These options allow a circuit to select between dense and sparse representations, direct and iterative methods, custom and GSL-backed dense direct solvers, and the supported transient integration methods.

Different solver configurations can therefore be run on the same circuit for numerical comparison.

## Output files

Simulation results are written under the `OUT/` directory in a subdirectory associated with the input netlist.

The exact output files depend on the analyses and solver options enabled by the netlist.

Results can be inspected directly or compared with reference solutions using the validation utility under `cmpr/`.

## Validation

The repository contains:

```text
cmpr/compare_results.py
```

for comparing simulator results against reference data.

The comparison utility reports numerical error statistics such as:

- mean absolute error
- mean relative error
- pass rate

Usage information is available in:

[`cmpr/README.md`](./cmpr/README.md)

Reference solutions for the included IBM benchmark circuits are stored under:

```text
IBM_SOLS/
```


### Validation results

The simulator was rebuilt from a fresh clone on Ubuntu and evaluated against the included IBM reference data.

| Benchmark | Analysis | Circuit scale | Analysis time | Peak memory | Numerical agreement |
|---|---|---:|---:|---:|---|
| IBM PG1 | Sparse LU DC operating point | 30,636 nodes, 44,943 MNA unknowns | 0.360 s | 53 MB | MAE 2.77e-7 V, 100% of non-zero reference nodes within 1% |
| IBM PG2 | Sparse LU DC operating point | 127,236 nodes, 127,565 MNA unknowns | 14.383 s | 338 MB | MAE 4.79e-8 V, 100% of non-zero reference nodes within 1% |
| IBM PG1T | Sparse transient, trapezoidal | 39,681 nodes, 54,265 MNA unknowns | 8.865 s | 87 MB | MAE 4.27e-6 V across 20,020 reference samples, 100% within 1% |

Analysis time is the simulator-reported interval measured after netlist parsing and includes the enabled numerical analyses and their output generation.

For DC validation, relative-error statistics exclude reference voltages effectively equal to zero, where relative error is not meaningful. Absolute-error statistics include the matched node voltages.

The PG2 benchmark was evaluated in sparse mode by adding `.options sparse` to a temporary local copy of the original benchmark netlist. The original benchmark file in the repository was not modified.

For the PG1T transient comparison, the 20 requested nodes were compared over the 1,001 reference time points from 0 to 10 ns.

These benchmark files and reference solutions are used as validation inputs and are not claimed as original benchmark data created by this project.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Simulator implementation |
| `include/` | Project headers and bundled utility headers |
| `csparse/` | Bundled third-party sparse matrix routines |
| `Part1_Netlists/` | Basic circuit and DC-analysis netlists |
| `Part3_Netlists/` | Sparse and iterative solver test netlists, including IBM power-grid cases |
| `Part6_Netlists/` | Transient-analysis netlists, including IBM benchmark cases |
| `IBM_SOLS/` | Reference solutions used for benchmark validation |
| `cmpr/` | Numerical output comparison tools |
| `NETLIST.md` | Detailed supported netlist syntax |
| `Makefile` | Build configuration |

## Implementation scope

The following functionality is implemented by the project code under `src/` and the project-specific headers under `include/`:

- SPICE-style parsing and internal circuit representation
- Modified Nodal Analysis assembly
- resistor, capacitor, inductor, voltage-source, and current-source handling
- DC operating-point analysis
- DC sweep control
- transient integration flow
- transient source evaluation
- custom dense LU factorization
- custom dense Cholesky factorization
- custom forward and backward substitution
- dense CG and BiCG solver logic
- sparse CG and BiCG solver logic
- Jacobi preconditioning
- output generation and analysis control

Third-party numerical and utility components are identified separately in the Third-party components section.

## Known limitations

- Diode, MOSFET, and BJT elements are parsed but are not stamped into the MNA system
- `.MODEL` parameters are parsed but are not used during simulation
- Controlled sources such as VCVS, VCCS, CCVS, and CCCS are not supported
- Nonlinear Newton-Raphson iteration is not implemented
- Maximum netlist line length is 256 characters
- Node, element, and model names support up to 32 characters
- PWL sources support up to 20 time-value pairs

Full syntax details are available in [`NETLIST.md`](./NETLIST.md).

## Academic context and provenance

This simulator was developed for **ECE513 - Circuit Simulation Algorithms (Αλγόριθμοι Προσομοίωσης Κυκλωμάτων)** at the **University of Thessaly** between **October 2025 and January 2026**.

The repository contains the simulator implementation together with circuit netlists, benchmark data, validation utilities, and third-party numerical components used during development and evaluation.

The project-specific implementation should be distinguished from the external components and benchmark material identified above. The included IBM benchmark netlists and reference solutions are validation data and are not claimed as original benchmark designs produced by this project.

The repository was later organized into a clearer public structure for portfolio and educational review.

## License and third-party material

A root MIT License is included for the original project code.

Bundled third-party components and external benchmark material remain subject to their respective upstream terms. The root project license does not replace or supersede the licensing terms, copyright notices, or redistribution requirements of third-party files.

Users redistributing the repository should preserve the applicable upstream notices and licensing information for bundled dependencies and benchmark material.