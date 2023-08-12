from topoly import alexander
import mdtraj as md
from Bio.PDB import PDBParser, PDBExceptions
import os


def check_file_extension(filename):
    """
    Function check the file extension.
    Args:
        filename (str):
            Name of the file to examine.

    Returns:
        The file extension.
    """
    _, file_extension = os.path.splitext(filename)
    return file_extension


def process_structure_and_calculate(file, chain_id, atom_list, closure, tries, max_cross):
    if check_file_extension(file) == ".pdb":
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("pdb", file)
        except PDBExceptions.PDBConstructionWarning:
            raise PDBExceptions.PDBConstructionWarning

        # looking for the name of first chain
        model = structure[0]
        first_chain = next(model.get_chains())
        first_chain_id = first_chain.id

        if chain_id is None:
            chain_id = [first_chain_id]

        # looking for the name of first atom
        if atom_list is None:
            model = structure[0]
            chain = model[first_chain_id]
            residue = next(chain.get_residues())
            atom_id = next(residue.get_atoms())

        chain_atoms = []
        for model in structure:
            for chain in model:
                if chain.id in chain_id:
                    for residue in chain:
                        for atom in residue:
                            if atom_list is None:
                                if atom.name == atom_id.element:
                                    chain_atoms.append(atom)
                            else:
                                if atom.name in atom_list:
                                    chain_atoms.append(atom)

        with open("output.nxyz", "w") as f:
            i = 0
            for atom in chain_atoms:
                coords = atom.get_coord()
                x, y, z = coords[0], coords[1], coords[2]

                f.write("{} {} {} {}\n".format(i, x, y, z))
                i += 1
        f.close()
        output_file = "output.nxyz"
        pdb = True
    else:
        output_file = file
        pdb = False

    knotcore_result = count_knotcore(output_file, closure=closure, tries=tries, max_cross=max_cross)
    if pdb:
        os.remove("output.nxyz")

    return knotcore_result


def load_structure(file, top_file):
    """
    Function reads the structure in one of the accepted formats.
    Args:
        file (str):
                The path to the structure in accepted format: .pdb, .xyz or .xtc.
        top_file (str):
                If the file is not given in PDB format, it is required to specify an extra file in order to conduct
                the analysis. This is because the XYZ and XTC format does not contain topology information.
                Pass in either the path to a PDB file, a trajectory, or a topology to supply this information.

    Returns:
        The resulting trajectory, as a md.Trajectory object (trajectorymd.Trajectory).
    """

    extension = check_file_extension(file)
    t = ''
    if extension == ".pdb":
        t = md.load_pdb(file)
    if extension == ".xyz":
        if top_file is None:
            raise ValueError("This format of file requires an additional file 'top_file'.")
        else:
            t = md.load_xyz(file, top=top_file)
    if extension == ".xtc":
        if top_file is None:
            raise ValueError("This format of file requires an additional file 'top_file'.")
        else:
            t = md.load(file, top=top_file)

    return t


def get_lider_from_dict(knot_dict):
    kn, pr = '0_1', 0
    for (k, p) in knot_dict.items():
        if p > pr:
            kn, pr = k, p
    return [kn, pr]


def find_knotcore_simple(chain, gap=1, closure=1, tries=30, cutoff=0.42, max_cross=20):
    """
    This function is a slightly modified version of the code from original function authored by Dr Wanda Niemyska.
    The function is used with the author's permission. In the future, there are plans to include the knot core value
    calculation function in the Topoly package, and then it will be more appropriate to use it.

    Function finds the knotcore with a simple approach:
        1. checks the knot type of the whole structure (main_knot), returns None if it is unknotted;
        2. cuts atoms from one end of the chain as long as the shortened chain is still knotted
           (more precisely - forms main_knot type);
           allows the situation that a certain chain is unknotted, but after cutting off another atom,
           the knot returns again (gap parameter); we get cut_beg, cut_end (how much we can cut);
        3. checks if after cutting cut_beg and cut_end atoms from the beginning and end of the structure,
           respectively, there is a main_knot formed - if yes, returns this as a knotcore, if not decreases cut_beg
           and cut_end by one, and tries again;

    It is possible to specify indexes of first and last atoms for .xyz and .nxyz files.

    Args:
        chain (str):
                Structure for which the knot core value is to be calculated, in the .xyz or .nxyz format.
        gap (int, optional):
                The maximum number of frames in which there is no knot, or for some reason, the knot core value cannot
                be calculated, for which the computations will not be interrupted (i.e., the interruption will not be
                considered as the end of the structure, but rather as a discontinuity).
                Default: 1.
        cutoff (float, optional):
                The parameter used during the application of random closures. It determines the threshold above which
                we consider that a knot has formed.
                Default: 0.42.

    Returns: None or
            (begin_of_knotcore, end_of_knotcore), where these are ids from the file (not necessarily
             starting from 0/1)
    """
    id_beg = 0
    id_end = 0

    def find_subknot(beg, end):
        kn = alexander(chain, chain_boundary=[[beg, end]], closure=closure, tries=tries, max_cross=max_cross,
                       run_parallel=False)
        kn = kn[(beg, end)]
        if closure > 1:
            kn, prob = get_lider_from_dict(kn)
        else:
            prob = 1
        return kn, prob

    if chain.endswith('.nxyz'):
        res = []
        with open(chain) as file_nxyz:
            for line in file_nxyz:
                line_res = line.split()[0]
                res.append(line_res)
        id_beg, id_end = int(res[0]), int(res[-1])
    if chain.endswith('.xyz'):
        id_beg = 0
        res = []
        with open(chain) as file_xyz:
            for line in file_xyz:
                res.append(line)
        id_end = len(res) - 1

    main_knot, prob = find_subknot(id_beg, id_end)
    if main_knot == '0_1' or prob < cutoff:
        return None

    cut_beg = 0
    act_kn, prob, act_gap = main_knot, 1, 0
    while id_beg + cut_beg < id_end - 5 and ((act_kn == main_knot and prob >= cutoff) or act_gap <= gap):
        cut_beg += 1
        act_kn, prob = find_subknot(id_beg + cut_beg, id_end)
        if act_kn == main_knot and prob >= cutoff:
            act_gap = 0
        else:
            act_gap += 1
    cut_end = 0
    act_kn, prob, act_gap = main_knot, 1, 0
    while id_beg + 5 < id_end - cut_end and ((act_kn == main_knot and prob >= cutoff) or act_gap <= gap):
        cut_end += 1
        act_kn, prob = find_subknot(id_beg, id_end - cut_end)
        if act_kn == main_knot and prob >= cutoff:
            act_gap = 0
        else:
            act_gap += 1

    act_kn = '0_1'
    while (act_kn != main_knot or prob < 0.8 * cutoff) and cut_beg + cut_end > 0:
        if cut_beg > 0:
            cut_beg -= 1
        if cut_end > 0:
            cut_end -= 1
        act_kn, prob = find_subknot(id_beg + cut_beg, id_end - cut_end)

    res_list = (id_beg + cut_beg, id_end - cut_end)
    return res_list


def count_knotcore(file, closure, tries, max_cross):
    res = find_knotcore_simple(file, closure=closure, tries=tries, max_cross=max_cross)
    return res
