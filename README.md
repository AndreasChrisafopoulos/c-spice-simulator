# C SPICE Simulator

## Overview
This project implements a SPICE-like circuit simulator written entirely in C. It is designed to parse standard netlist files and perform core circuit simulations, emphasizing algorithmic implementation from scratch, including both direct and iterative matrix solvers.

## Features
*   **Netlist Parsing:** Robust, SPICE-like syntax parser with case-insensitive processing.
*   **Multiple Analyses:** Supports DC Operating Point, DC Sweeps, and Transient Analysis.
*   **Custom Solvers:** Built-in direct solvers and iterative solvers to evaluate performance trade-offs.
*   **Sparse Matrix Support:** Optimized memory and computation for large-scale node networks.
*   **Output Generation:** Generates segregated output data for easy plotting and direct vs. iterative solver comparisons.

## Supported Circuit Elements
*   **Passive:** Resistor (`R`), Capacitor (`C`), Inductor (`L`)
*   **Active:** Diode (`D`), BJT Transistor (`Q`), MOS Transistor (`M`)
*   **Sources:** Independent Voltage (`V`) and Current (`I`) Sources (DC and Transient specs: EXP, SIN, PULSE, PWL)

## Supported Analyses
*   **DC Operating Point** (`.op`)
*   **DC Sweep** (`.dc`): Up to 32 parallel sweep analyses per netlist.
*   **Transient Analysis** (`.tran`)

## Numerical Methods
*   **LU Factorization**
*   **Cholesky Decomposition**
*   **BiCG (BiConjugate Gradient)**
*   **CG (Conjugate Gradient)**
*   **Sparse LU** (Enabled via `.options sparse`)

## Project Architecture

```text
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