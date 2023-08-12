#!/usr/bin/python3
# -*- coding: utf-8 -*-
from packages.knotcore import process_structure_and_calculate
import argparse


def calculate_pdb_knotcore(file, chain_id=None, atom_list=None, closure=1, tries=20, max_cross=30):
    """
    Function calculates the knot core value for the given structure in .pdb, .xyz or .nxyz format. If file is in PDB
    format, then is possible to choose chain (or chains) and atom (or atoms), which should be taken into account during
    the calculation.

    Args:
        file (str):
                The path to the structure in one of accepted format.
        chain_id (list of strings, optional):
                Example: ['A', 'B']
                If inputting .pdb file, which chain from the structure should be used.
                If none, the first chain is taken into account.
                Default: None.
        atom_list (list of strings, optional):
                Example: ['C', 'N']
                If inputting .pdb file, which atoms from the structure/ selected chain should be used.
                If none, the first atom is taken into account.
                Default: None.
        closure (str, optional):
                The method to close the chain. Viable options are the parameters of the Closure class (in
                topoly.params).
                Default: Closure.MASS_CENTER.
        tries (int, optional):
                The number of tries for stochastic closure methods.
                Default: 20.
        max_cross (int, optional):
                The maximal number of crossings after reduction to start
                the polynomial calculation. Default: 30.

    Returns: The knot core value
             None, if failed to calculate
    """
    return process_structure_and_calculate(file, chain_id, atom_list, closure, tries, max_cross)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate knot core range for a given structure.')
    parser.add_argument('file', type=str, help='Path to the structure file in PDB format.')
    parser.add_argument('-i', '--chain_id', nargs='+', default=None,
                        help='If main file is in PDB format. List of chain IDs to be used.')
    parser.add_argument('-a', '--atom_list', nargs='+', default=None,
                        help='If main file is in PDB format. List of atom names to be used.')
    parser.add_argument('-c', '--closure', type=int, default=1,
                        help='The method to close the chain. Viable options are parameters of the Closure'
                             ' class (in topoly.params).')
    parser.add_argument('-t', '--tries', type=int, default=20, help='Number of tries for stochastic closure methods.')
    parser.add_argument('-m', '--max_cross', type=int, default=30,
                        help='Maximal number of crossings for polynomial calculation.')

    args = parser.parse_args()
    if args.chain_id is not None:
        chain_id = list(args.chain_id)
    else:
        chain_id = None

    if args.atom_list is not None:
        atom_list = list(args.atom_list)
    else:
        atom_list = None

    knotcore_value = calculate_pdb_knotcore(args.file, chain_id, atom_list, args.closure, args.tries,
                                            args.max_cross)
    print(knotcore_value)
