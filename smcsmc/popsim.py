'''
This submodule interacts with the PopSim consortium's analysis respository

found here: https://github.com/popgensims/analysis

These functions are not explicitly exported for that reason.
'''


import tskit
import os
from .model import *
import pdb
import numpy as np
import pandas as pd
import glob


def convert_smcsmc_output(results_file, output, generation_time, iter):
    out_fp = open(output, "w")
    in_fp = open(results_file, "r")
    header = in_fp.readline()
    out_fp.write("label,x,y,plot_type,plot_num")
    for line in in_fp:
        result = line.split()
        if int(result[0]) == iter:
            if result[4] == "Coal":
                if result[12] == "-1":
                    time_years = float(result[2]) * generation_time
                    size = float(result[10])
                    out_fp.write(f"pop0,{time_years},{size},path,0\n")
    out_fp.close
    return(None)

