# Bachelor's thesis: Visualization of the process of knot formation in the trajectory of (bio)polymers folding
PL:

Repozytorium zawiera pliki będące częścią pracy licencjackiej, realizowanej na wydziale Matematyki, Informatyki i Mechaniki Uniwersytetu Warszawskiego na kierunku Bioinformatyka i Biologia Systemów. 

Temat Pracy: Wizualizacja procesu powstawania węzłów w trajektorii zwijania się (bio)polimerów

Autor: Aleksandra Sekuła 

Promotor: Dr Wanada Niemyska 

Rok akademicki: 2022/2023 

W folderze traj_analysis znajdują się wszystkie, niezbędne pliki do użycia skryptu. Dokumentacja programu: docs/html/index.html

W folderze examples znajduje się przykładowy wykres, wygenerwoany za pomocą skryptu dla jednej z analizowanych trajektorii oraz kilka plików z trajektoriami do testów.

ENG:
The repository contains files that are part of a bachelor's thesis conducted at the Faculty of Mathematics, Informatics, and Mechanics of the University of Warsaw, in the field of Bioinformatics and Systems Biology.

The subject of the thesis: Visualization of the process of knot formation in the trajectory of (bio)polymers folding

Author: Aleksandra Sekuła 

Supervisor:  Dr Wanada Niemyska

Academic year: 2022/2023

In the 'traj_analysis' folder, you will find all the necessary files to use the script. Program documentation: docs/html/index.html

The 'examples' folder contains a sample plot generated using the script for one of the analyzed trajectories and several files with trajectories for testing.

# Traj_analysis
Script analyzes the trajectory of the (bio)polymers. It finds frames of the trajectory in which knot forms based on the given conditions, analyzes the process of knot formation and plots the range chart of the knot core throughout the entire trajectory of the analyzed structure.

# Usage 
Currently, usage is only possible after downloading the 'traj_analysis' folder and directly utilizing its code. The script consists of two main functions: analyze_trajectory and calculate_pdb_knotcore (for more information, refer to the documentation: docs/html/index.html). They can be imported into your own program:
```python
import traj_analysis
import calculate_knotcore

# returns dictionary with the results of the analysis:  {402: ['3_1', None, (10, 80), 0, 1]}
analyze_trajectory("examples/traj.pdb", nterminus=True, nat_knotcore=(13, 80))

# calculate knot core value for the given structure: (13,80)
calculate_pdb_knotcore("examples/2efv.pdb")
```
or use from the command line:
```python
Analysis of the trajectory.

positional arguments:
  file                  Path to the structure file in accepted format: .pdb, .xyz or .xtc.
  nterminus             The end of a structure that goes through a loop when a knot is formed. True if closer to N-terminus, False if closer to C-terminus.

optional arguments:
  -h, --help            show this help message and exit
  -o TOP_FILE, --top_file TOP_FILE
                        Path to a PDB file, a trajectory, or a topology to supply information for non-PDB formats of the main file.
  -n NAT_KNOTCORE [NAT_KNOTCORE ...], --nat_knotcore NAT_KNOTCORE [NAT_KNOTCORE ...]
                        The knot core range in the native form of the structure. Usually the knot core range is given as a tuple, but for the program to work correctly, the first
                        value must be given, followed by a space and the second value. Example: for knot core range (9, 87), program must get: -n 9 87
  -g MIN_GAP, --min_gap MIN_GAP
                        The minimum number of frames before a knot is considered to have formed, in which 80 percent of those frames do not contain a knot.
  -s SCOPE, --scope SCOPE
                        The minimum consecutive number of frames in which a knot occurs in order to consider that a knot has actually formed.
  -k MIN_KNOT, --min_knot MIN_KNOT
                        The minimum number of frames over which a node is present to consider it a stable node in the analysis.
  -c CLOSURE, --closure CLOSURE
                        The method to close the chain. Viable options are parameters of the Closure class (in topoly.params).
  -t TRIES, --tries TRIES
                        Number of tries for stochastic closure methods.
  -m MAX_CROSS, --max_cross MAX_CROSS
                        Maximal number of crossings after reduction to start polynomial calculation.
  -d, --draw_plot       Generate a plot.
  -p PLOT_FILENAME, --plot_filename PLOT_FILENAME
                        Name of the plot file.
  -l PLOT_SCOPE, --plot_scope PLOT_SCOPE
                        The step value determines how frequently the knot core value will be calculated and plotted on the graph.
  -e, --debug           Enable debug mode.
  -f, --full_output     Display full analysis results.

```

```python
calculate_knotcore.py -h
Calculate knot core range for a given structure.

positional arguments:
  file                  Path to the structure file in PDB format.

optional arguments:
  -h, --help            show this help message and exit
  -i CHAIN_ID [CHAIN_ID ...], --chain_id CHAIN_ID [CHAIN_ID ...]
                        If main file is in PDB format. List of chain IDs to be used.
  -a ATOM_LIST [ATOM_LIST ...], --atom_list ATOM_LIST [ATOM_LIST ...]
                        If main file is in PDB format. List of atom names to be used.
  -c CLOSURE, --closure CLOSURE
                        The method to close the chain. Viable options are parameters of the Closure class (in topoly.params).
  -t TRIES, --tries TRIES
                        Number of tries for stochastic closure methods.
  -m MAX_CROSS, --max_cross MAX_CROSS
                        Maximal number of crossings for polynomial calculation.
```

