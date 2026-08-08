# IBM Power-Grid Benchmarks

IBM power-grid benchmark netlists and reference solutions were used as external validation data for this project.

The benchmark files are not redistributed with this repository.

Official source:

https://web.ece.ucsb.edu/~lip/PGBenchmarks/ibmpgbench.html

## Benchmarks used for reported validation

- `ibmpg1.spice` for IBM PG1 DC operating-point validation
- `ibmpg2.spice` for IBM PG2 DC operating-point validation
- `ibmpg1t.spice` for IBM PG1T transient validation
- the corresponding official reference solution files

Downloaded DC netlists can be placed under `Part3_Netlists/`.

Downloaded transient netlists can be placed under `Part6_Netlists/`.

Reference solutions can be placed in a local `IBM_SOLS/` directory for use with the comparison utility under `cmpr/`.

## Validation notes

IBM PG2 was evaluated using a temporary local copy with `.options sparse` added so that the sparse solver path was selected. The original downloaded benchmark was not modified.

IBM PG1T was compared over the 1,001 official reference timestamps from 0 to 10 ns.

The benchmark files and reference solutions remain external material and are subject to the terms of their original source.
