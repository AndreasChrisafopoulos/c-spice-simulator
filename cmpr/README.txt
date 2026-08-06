Comparison Tool for Simulation Outputs
======================================

This folder contains a Python script used to compare two simulation
output files and evaluate their numerical agreement.

The comparison can be used for:
- DC Operating Point analysis
- DC Sweep analysis
- Transient analysis


--------------------------------------------------
Folder Contents
--------------------------------------------------

compare_results.py   : Python comparison script
my_output.txt        : Output produced by our simulator
solution.txt         : Reference output
README.txt           : Usage instructions


--------------------------------------------------
Preparing the Input Files
--------------------------------------------------

To perform a comparison, copy the values you want to compare into:

- my_output.txt       (results from our simulator)
- solution.txt        (reference results)

The two files must correspond to the same analysis and contain
values in the same order.


--------------------------------------------------
Usage
--------------------------------------------------

Navigate to the cmpr directory and run:

For DC Operating Point (OP) analysis:

    python3 compare_results.py my_output.txt solution.txt op

For DC Sweep or Transient analysis:

    python3 compare_results.py my_output.txt solution.txt

--------------------------------------------------
Output
--------------------------------------------------

The script reports:

- Mean Absolute Error
- Mean Relative Error
- Pass Rate based on relative error threshold

Example:

    Mean Absolute Error : 1.68e-06
    Mean Relative Error : 0.003%
    Pass Rate (≤1.0%): 100.0%