import sys
import os
from os import path
from zipfile import ZipFile

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from baseutil import *
import logging

import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobrakbase.core.kbasefba import FBAModel
from cobra.io import write_sbml_model, read_sbml_model
from modelseedpy import AnnotationOntology, MSPackageManager, ModelSEEDBiochem,MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
from modelseedpy.core.msprobability import MSProbability
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msensemble import MSEnsemble
from mscommunity import MSCommunity
from modelseedpy.helpers import get_template

class SludgeCommUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self,"Sludge")
        # self.msseedrecon()
    
    def create_phenotypeset_from_compounds(
        self,
        compounds,
        base_media=None,
        base_uptake=0,
        base_excretion=1000,
        global_atom_limits={},
        type="growth"
    ):
        cpd_hash = {}
        for cpd in compounds:
            cpd_hash[cpd] = 10
        return MSGrowthPhenotypes.from_compound_hash(
            cpd_hash,
            base_media=base_media,
            base_uptake=base_uptake,
            base_excretion=base_excretion,
            global_atom_limits=global_atom_limits,
            type=type
        )

util = SludgeCommUtil() 