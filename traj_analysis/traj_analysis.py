#!/usr/bin/python3
# -*- coding: utf-8 -*-

from packages.traj import *
from packages.plot import Plot
import argparse


def analyze_trajectory(file, nterminus, top_file=None, nat_knotcore=None, min_gap=10, scope=10, min_knot=100,
                       closure=1, tries=20, max_cross=15, draw_plot=False, plot_filename="knotcore_plot",
                       plot_scope=100, debug=False, full_output=False):
    """
    Function finds frames in which knot forms based on the given conditions. It evaluates how the knot was
    formed (via slipknot/normally) and whether the loop was +/- in its place at the moment, when the knot was formed.
    Additionally, it can plot the range of knot core values for the entire trajectory of the molecule.

    Args:
        file (str):
                The path to the structure in accepted format: .pdb, .xyz or .xtc.
        nterminus (bool):
                The end of a structure that goes through a loop when a knot is formed.
                True: if closer to the N-terminus (N-terminus is the start of an amino acid chain (protein or
                polypeptide)).
                False: if closer to C-terminus (C-terminus is the end of an amino acid chain (protein or polypeptide)).
        top_file (str, optional):
                If the file is not given in .pdb format, it is required to specify an extra file in order to conduct
                the analysis. This is because the .xyz and .xtc format does not contain topology information.
                Pass in either the path to a .pdb file, a trajectory, or a topology to supply this information.
                Default: None.
        nat_knotcore: (tuple of (int, int), optional)
                The knot core range in the native form of the structure. When given, the program will determine the
                position of the loop after knotting relative to the position in the native frame.
                Default: None.
        min_gap (int, optional):
                The minimum number of frames before a knot is considered to have formed, in which 80 % of those frames
                do not contain a knot.
                Default: 10.
        scope (int, optional):
                The minimum number of consecutive frames in which a node occurs to be considered as having actually
                formed a node.
                Default: 10.
        min_knot (int, optional):
                The minimum number of frames over which a node is present to consider it a stable node in the analysis.
                Default: 100.
        closure (str, optional):
                The method to close the chain. Viable options are the parameters of the Closure class (in
                topoly.params).
                Default: Closure.MASS_CENTER.
        tries (int, optional):
                The number of tries for stochastic closure methods.
                Default: 20.
        max_cross (int, optional):
                The maximal number of crossings after reduction to start
                the polynomial calculation. Default: 15.
        draw_plot (bool, optional):
                If to generate the plot. Default: False.
        plot_filename (str, optional):
                The name of the file with the plot results.
                Default: 'knotcore_plot'. The resulting plot is returned to source.
        plot_scope (int, optional):
                The step value determines how frequently the knotcore value will be calculated and plotted on the graph.
                Default: 100.
        debug (bool, optional):
                The debug mode.
                Default: False.
        full_output (bool, optional):
                How the analysis results are displayed.
                True, if full information with text.
                False, if just the result.
                Default: False

    Returns:
    Dictionary of frames, when a knot is tied as keys and as value the result of the analysis. The result
    is a list where the following values are in order: knot type, number of unknotting frame, knot core value in
    frame, where the knot was knotted, number that indicates way of knotting, and number that indicates the position
    od the loop.

    e.g. {402: ['3_1', None, (10, 80), 0, 1]}
    e.g. if full_output = True: {402: {'Knot type': '3_1', 'Unknotting frame': None, 'Knot core range': (10, 80),
                                       'Knotting via slipknot': True, 'Loop behavior': 'loop is in place'}}
    Symbols used in the code for specific cases:
    - Ways of knotting: 0 - knotting via slipknot
                       1 - knotting direct.
    - Behavior of the loop: 0 - loop tightens.
                            1 - loop is in place.
                            2 - loop expands.

    If plot=True, then plot of the knot core range for the entire trajectory of the molecule. If the 'plot_filename'
    parameter has not been changed, the file will be created in the current directory. Plot is saved in html format.
    """
    if debug:
        print('Analyzing the trajectory with parameters:\n' + str(locals()))

    t = load_structure(file, top_file)

    try:
        lx = list(t.xyz[::])
    except AttributeError as e:
        print("Error occurred during loading data: ", e, ".")
        return None

    trajectory = Traj(lx, t.n_atoms - 1, len(lx) - 1, min_gap, scope, min_knot, nterminus, nat_knotcore, closure, tries,
                      max_cross, debug)

    knot_dict = trajectory.calculate(full_output)

    if draw_plot:
        if len(knot_dict) != 0:
            traj_plot = Plot(trajectory, plot_filename, plot_scope, debug)
            traj_plot.draw_plot()
        elif debug:
            print("The program did not detect any knots in the molecule. \n"
                  "Nothing to plot.")

    return knot_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analysis of the trajectory.')
    parser.add_argument('file', type=str, help='Path to the structure file in accepted format: .pdb, .xyz or .xtc.')
    parser.add_argument('nterminus', type=bool, help='The end of a structure that goes through a loop when a knot is '
                                                     'formed. True if closer to N-terminus, False if closer to '
                                                     'C-terminus.')
    parser.add_argument('-o', '--top_file', type=str, default=None,
                        help='Path to a PDB file, a trajectory, or a topology'
                             ' to supply information for non-PDB formats of the main file.')
    parser.add_argument('-n', '--nat_knotcore',  nargs='+', type=int, default=None,
                        help='The knot core range in the native form of the structure. Usually the knot core range is'
                             ' given as a tuple, but for the program to work correctly, the first value must be given,'
                             ' followed by a space and the second value. '
                             'Example: for knot core range (9, 87), program must get: -n 9 87')
    parser.add_argument('-g', '--min_gap', type=int, default=10, help='The minimum number of frames before a knot is'
                                                                      ' considered to have formed, in which 80 percent of'
                                                                      ' those frames do not contain a knot.')
    parser.add_argument('-s', '--scope', type=int, default=10, help='The minimum consecutive number of frames in which'
                                                                    ' a knot occurs in order to consider that a knot'
                                                                    ' has actually formed.')
    parser.add_argument('-k', '--min_knot', type=int, default=100,
                        help='The minimum number of frames over which a node is present to consider it a stable node'
                             ' in the analysis.')
    parser.add_argument('-c', '--closure', type=int, default=1,
                        help='The method to close the chain. Viable options are parameters of the Closure'
                             ' class (in topoly.params).')
    parser.add_argument('-t', '--tries', type=int, default=20, help='Number of tries for stochastic closure methods.')
    parser.add_argument('-m', '--max_cross', type=int, default=15, help='Maximal number of crossings after reduction '
                                                                        'to start polynomial calculation.')
    parser.add_argument('-d', '--draw_plot', action='store_true', help='Generate a plot.')
    parser.add_argument('-p', '--plot_filename', type=str, default='knotcore_plot', help='Name of the plot file.')
    parser.add_argument('-l', '--plot_scope', type=int, default=100, help='The step value determines how frequently'
                                                                          ' the knot core value will be calculated and'
                                                                          ' plotted on the graph.')
    parser.add_argument('-e', '--debug', action='store_true', help='Enable debug mode.')
    parser.add_argument('-f', '--full_output', action='store_true', help='Display full analysis results.')

    args = parser.parse_args()
    nat_tuple = tuple(args.nat_knotcore)

    res = analyze_trajectory(args.file, args.nterminus, args.top_file, nat_tuple, args.min_gap, args.scope,
                             args.min_knot, args.closure, args.tries, args.max_cross, args.draw_plot,
                             args.plot_filename, args.plot_scope, args.debug, args.full_output)
    print(res)
