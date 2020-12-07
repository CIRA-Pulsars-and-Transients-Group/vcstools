

def read_rmsynth_out(filename):
    """
    Reads the ouput file from rm_synth_pipe()

    Parameters:
    -----------
    filename:
        The name of the file ouput from rm_synth_pipe()

    Returns:
    --------
    rm_dict: dictionary
        Contains the following keys:
        i: dictionary
            There are i entries where i is the number of different phase ranges
            Contains the following keys:
            rm: float
                The rotation measure of this run
            rm_e: float
                The uncertainty in rm
            phase_ranges: tuple
                The range of phases used for this run
    """
    rm_dict={}
    f = open(filename, "r")
    lines = f.readlines()
    j=0
    for line in lines:
        if line.split()[0] == "Phase:":
            rm_dict[str(j)] = {}
            rm_dict[str(j)]["phase_range"] = (float(line.split()[1]), float(line.split()[-1]))
        if line.split()[0] == "Rotation":
            rm_dict[str(j)]["rm"] = float(line.split()[2])
            rm_dict[str(j)]["rm_e"] = float(line.split()[-1])
            j += 1
    return rm_dict