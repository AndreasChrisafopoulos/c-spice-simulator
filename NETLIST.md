# Netlist Specification and Limitations

This document describes the supported netlist format, circuit elements, analyses, options, and known limitations of the circuit simulator.

## General Netlist Format

- SPICE-like syntax
- One element or directive per line
- Lines starting with `*` are treated as comments
- Parsing is case-insensitive
- Parsing stops at the `.end` directive
- Node names are alphanumeric strings
- Ground node must be named `0`

## Supported Circuit Elements

### Resistor
```
R<name> <node1> <node2> <value>
```

### Capacitor
```
C<name> <node1> <node2> <value>
```

### Inductor
```
L<name> <node1> <node2> <value>
```

### Independent Voltage Source
```
V<name> <node+> <node-> <dc_value> [transient_spec]
```

### Independent Current Source
```
I<name> <node+> <node-> <dc_value> [transient_spec]
```

### Diode
```
D<name> <node+> <node-> <model> [area=<value>]
```
The `area` parameter is optional. **Parsed only** — not yet stamped into the MNA matrix (see Notes and Limitations).

### MOS Transistor
```
M<name> <drain> <gate> <source> <bulk> <model> L=<value> W=<value>
```
**Parsed only** — not yet stamped into the MNA matrix (see Notes and Limitations).

### BJT Transistor
```
Q<name> <collector> <base> <emitter> <model> [area=<value>]
```
The `area` parameter is optional. **Parsed only** — not yet stamped into the MNA matrix (see Notes and Limitations).

## Supported Analyses

### DC Operating Point
```
.op
```

### DC Sweep
```
.dc <source_name> <start> <end> <step>
```

### Transient Analysis
```
.tran <time_step> <final_time>
```

## Plot / Print Commands

Plot and print commands may include multiple nodes on the same line.

```
.plot V(n1) V(n2) V(n3)
.print V(out)
```

Context-sensitive behavior:
- After a `.dc` directive, plot commands refer to the sweep variable
- After a `.tran` directive, plot commands refer to time

Each DC sweep or transient analysis supports up to 32 plotted nodes.

## Transient Source Specifications

Independent voltage and current sources may include one transient specification in addition to their DC value.

### EXP (Exponential)
```
EXP(i1 i2 td1 tc1 td2 tc2)
```
- `i1` — initial value
- `i2` — pulsed value
- `td1` — rise delay time
- `tc1` — rise time constant
- `td2` — fall delay time
- `tc2` — fall time constant

### SIN (Sinusoidal)
```
SIN(i1 ia freq td df phase)
```
- `i1` — DC offset
- `ia` — amplitude
- `freq` — frequency (Hz)
- `td` — delay time
- `df` — damping factor
- `phase` — phase (degrees)

### PULSE (Pulse)
```
PULSE(i1, i2, td, tr, tf, pw, per)
```
- `i1` — initial value
- `i2` — pulsed value
- `td` — delay time
- `tr` — rise time
- `tf` — fall time
- `pw` — pulse width
- `per` — period

### PWL (Piecewise Linear)
```
PWL(t1 i1) (t2 i2) ... (tn in)
```
- `ti` — time point
- `ii` — source value at time `ti`

Notes:
- Up to 20 (time, value) pairs are supported.
- If the simulation time is outside the defined PWL range, the value at the nearest endpoint is used.

General notes:
- Sources with transient specifications must also include a DC value, used for DC operating point and DC sweep analyses.
- Each source may include only one transient specification.

## Options (`.options`)

Multiple options may be specified on a single line.

```
.options sparse spd iter method=be itol=1e-6
```

| Option | Effect |
|---|---|
| `sparse` | Enables sparse matrix representation |
| `spd` | Enables Cholesky factorization (direct) / CG (iterative) |
| `iter` | Enables iterative solution methods |
| `custom` | Enables the from-scratch solver implementation instead of GSL |
| `method=tr` | Trapezoidal method for transient analysis (default) |
| `method=be` | Backward Euler method for transient analysis |
| `itol=<value>` | Sets the iterative solver tolerance |

## Notes and Limitations

### General
- Ground node (`0`) always has node ID 0
- Controlled sources (VCVS, VCCS, CCVS, CCCS) are not supported
- Diode, MOSFET, and BJT elements are parsed (nodes, model name, `area`/`L`/`W`) but are **not stamped into the MNA matrix** — they do not currently affect simulation results
- `.MODEL` parameters are parsed but not validated or used

### Netlist Parsing
- Maximum netlist line length is 256 characters. Long `.plot` or `.print` commands must be split across multiple lines
- Maximum name length for nodes, elements, and models is 32 characters

### Analysis Limits
- Up to 32 DC sweep analyses (`.dc`) are supported per netlist
- Each DC sweep or transient analysis supports up to 32 plotted nodes

### Transient Sources
- Supported transient source types: `EXP`, `SIN`, `PULSE`, `PWL`
- PWL sources support up to 20 (time, value) pairs
- Each source may include only one transient specification

### Supported Elements Summary

| Symbol | Element | Simulated |
|---|---|---|
| `R` | Resistor | Yes |
| `C` | Capacitor | Yes |
| `L` | Inductor | Yes |
| `V` | Independent Voltage Source | Yes |
| `I` | Independent Current Source | Yes |
| `D` | Diode | Parsed only |
| `M` | MOS Transistor | Parsed only |
| `Q` | BJT Transistor | Parsed only |