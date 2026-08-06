# Circuit Simulator

This project implements a SPICE-like circuit simulator supporting
DC, DC sweep, and transient analysis.

--------------------------------------------------
Building
--------------------------------------------------

Compile the project using:

    make

--------------------------------------------------
Running
--------------------------------------------------

Run the simulator with two arguments:

    ./project <netlist_folder_id> <netlist_file>

Example:

    ./project 3 part3_simple.cir

The first argument selects the netlist folder:
    1 → Part 1 netlists
    3 → Part 3 netlists
    6 → Part 6 netlists

The second argument is the name of the netlist file
inside the selected folder.

--------------------------------------------------
Output Files
--------------------------------------------------

All simulation results are written under the OUT/ directory.

For each executed netlist, two output folders are created:

- outputfiles0  
  Results produced using direct solvers.

- outputfiles1  
  Results produced using iterative solvers
  (enabled via the `.options iter` directive).

This separation allows easy comparison between
direct and iterative solution methods.

--------------------------------------------------
Netlist Specification
--------------------------------------------------

The supported netlist format, elements, analyses,
options, and limitations are documented in:

    NETLIST.md
