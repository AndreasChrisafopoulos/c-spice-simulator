# C SPICE Simulator

## Introduction
This project implements a SPICE-like circuit simulator written entirely in C. It is designed to parse standard netlist files and perform core circuit simulations, emphasizing algorithmic implementation from scratch, including both direct and iterative matrix solvers.

## Features
*   **Netlist Parsing:** Robust, SPICE-like syntax parser with case-insensitive processing.
*   **Multiple Analyses:** Supports DC Operating Point, DC Sweeps, and Transient Analysis.
*   **Custom Solvers:** Built-in direct solvers and iterative solvers to evaluate performance trade-offs.
*   **Sparse Matrix Support:** Optimized memory and computation for large-scale node networks.
*   **Output Generation:** Generates segregated output data for easy plotting and direct vs. iterative solver comparisons.

## Supported Analyses & Components
*   **Analyses:** DC Operating Point (`.op`), DC Sweep (`.dc`), Transient Analysis (`.tran`).
*   **Components:** Resistors, Capacitors, Inductors, Diodes, BJTs, MOS Transistors.
*   **Sources:** Independent Voltage and Current Sources (DC and Transient specs: EXP, SIN, PULSE, PWL).

*For a complete and detailed breakdown of syntax and component limitations, please see [NETLIST.md](NETLIST.md).*

## Project Architecture

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

## Build Instructions

Compile the project using:

    make

## Running the Simulator

Run the simulator with two arguments:

    ./project <netlist_folder_id> <netlist_file>

**Example:**

    ./project 3 part3_simple.cir

**Folder ID Mapping:**
*   `1` → Part 1 netlists
*   `3` → Part 3 netlists
*   `6` → Part 6 netlists

## Output Files

All simulation results are written under the `OUT/` directory. For each executed netlist, two output folders are created:

*   `outputfiles0/` : Results produced using **direct solvers**.
*   `outputfiles1/` : Results produced using **iterative solvers** (enabled via the `.options iter` directive).

This separation allows easy comparison between direct and iterative solution methods.

## Academic Context
This project was developed to explore Electronic Design Automation (EDA) algorithms, focusing on the numerical methods (LU, Cholesky, BiCG, CG) required to reliably solve large, non-linear systems.