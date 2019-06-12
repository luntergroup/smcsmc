# this trick is from python-guide-pt-br.readthecods.io/en/latest/writing_structure

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../python')))

import smcsmc.populationmodels as populationmodels
import smcsmc.processrecombination as processrecombination
import smcsmc.execute as execute


