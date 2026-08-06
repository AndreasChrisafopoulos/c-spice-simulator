# Netlist Format and Specifications

This document describes the supported netlist format, circuit elements, analyses, options, and known limitations of the circuit simulator.

## General Netlist Format
*   **SPICE-like syntax**: One element or directive per line.
*   **Comments**: Lines starting with `*` are treated as comments.
*   **Case**: Parsing is case-insensitive.
*   **End**: Parsing stops at the `.end` directive.
*   **Ground**: Ground node must be named `0`.

---

## Passive Components

### Resistor
    R<name> <node1> <node2> <value>

### Capacitor
    C<name> <node1> <node2> <value>

### Inductor
    L<name> <node1> <node2> <value>

---

## Active Components & Sources

### Independent Voltage Source
    V<name> <node+> <node-> <dc_value> [transient_spec]

### Independent Current Source
    I<name> <node+> <node-> <dc_value> [transient_spec]

### Diode
*(The `area` parameter is optional)*
    D<name> <node+> <node-> <model> [area=<value>]

### MOS Transistor
    M<name> <drain> <gate> <source> <bulk> <model> L=<value> W=<value>

### BJT Transistor
*(The `area` parameter is optional)*
    Q<name> <collector> <base> <emitter> <model> [area=<value>]

---

## Supported Analyses

### DC Operating Point
    .op

### DC Sweep
    .dc <source_name> <start> <end> <step>

### Transient Analysis
    .tran <time_step> <final_time>

---

## Plot / Print Commands
Plot and print commands may include multiple nodes on the same line. Up to **32 plotted nodes** are supported.

    .plot V(n1) V(n2) V(n3)
    .print V(out)

---

## Transient Source Specifications
Sources with transient specifications must also include a DC value (used for DC operating point and sweeps). 

*   **EXP (Exponential)**: `EXP(i1 i2 td1 tc1 td2 tc2)`
*   **SIN (Sinusoidal)**: `SIN(i1 ia freq td df phase)`
*   **PULSE (Pulse)**: `PULSE(i1, i2, td, tr, tf, pw, per)`
*   **PWL (Piecewise Linear)**: `PWL(t1 i1) (t2 i2) ... (tn in)` *(Up to 20 pairs supported)*

---

## Options (`.options`)
Multiple options may be specified on a single line.

    .options sparse spd iter method=be itol=1e-6

*   `sparse`: Enables sparse matrix representation
*   `spd`: Enables Cholesky factorization or CG
*   `iter`: Enables iterative solution methods
*   `custom`: Enables custom solver mode
*   `method=tr` / `method=be`: Trapezoidal (default) or Backward Euler method
*   `itol=<value>`: Sets the iterative solver tolerance