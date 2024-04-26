import numpy as np
import pandas as pd
import networkx as nx
from scorers.SAScore import sascorer
from scorers.base_scorer import Scorer
from scorers.Sea_like_Tanimoto import sea_like_tanimoto_cal_rdkit
from rdkit.Chem import MolFromSmiles, MolToSmiles, Descriptors, rdmolops

class TanimotoSeaLikeCoef(Scorer):
    
    def __init__(self, params):
        self.params = params

    def score(self, population):

        smiles = population['smiles']
        smiles_rdkit = []
        for s in smiles:
            if pd.isnull(s):
                s = ''
            mol = MolFromSmiles(s)
            smi = MolToSmiles(mol,isomericSmiles=False)
            smiles_rdkit.append(smi)
        
        avg_sea_like_TC_values = []
        smi2_file = "/mnt/projects/RAS-CompChem/static/Mayukh/FNL_JTVAE/H5-B2_pocket/Raw-Data/H5-B2_masterlist_021324_salts_removed_smi.smi"
        threshold = self.params.tc_threshold
        for i in range(len(smiles_rdkit)):
            TC_sum = sea_like_tanimoto_cal_rdkit.sea_like_calc(smiles_rdkit[i],smi2_file,threshold)[0]
            if TC_sum > 0.0:
                avg_sea_like_TC_values.append(TC_sum/float(sea_like_tanimoto_cal_rdkit.sea_like_calc(smiles_rdkit[i],smi2_file,threshold)[1]))
            else:
                avg_sea_like_TC_values.append(TC_sum)
        
        SA_scores = []
        for i in range(len(smiles_rdkit)):
            SA_scores.append(-sascorer.calculateScore(MolFromSmiles(smiles_rdkit[ i ])))
        
        cycle_scores = []
        
        for i in range(len(smiles_rdkit)):
            cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(MolFromSmiles(smiles_rdkit[ i ]))))
            
            if len(cycle_list) == 0:
                cycle_length = 0
            else:
                cycle_length = max([ len(j) for j in cycle_list ])
            if cycle_length <= 6:
                cycle_length = 0
            else:
                cycle_length = cycle_length - 6
            cycle_scores.append(-cycle_length)
        
        SA_scores_normalized = (np.array(SA_scores) - np.mean(SA_scores)) / np.std(SA_scores)
        avg_sea_like_TC_values_normalized = (np.array(avg_sea_like_TC_values) - np.mean(avg_sea_like_TC_values)) / np.std(avg_sea_like_TC_values)
        cycle_scores_normalized = (np.array(cycle_scores) - np.mean(cycle_scores)) / np.std(cycle_scores)

        #STB: Added below three lines to handle nan's in cycle score. Any time you attempt np.nan + 5 (or 
        # any real number like that), the result will be nan, so summing SA score, logP values, and cycle scores 
        # was giving a fitness score of nan. Seems to only happen with certain molecules and not all?
        SA_scores_normalized[np.isnan(SA_scores_normalized)] = 0.0
        avg_sea_like_TC_values_normalized[np.isnan(avg_sea_like_TC_values_normalized)] = 0.0
        cycle_scores_normalized[np.isnan(cycle_scores_normalized)] = 0.0

        targets = SA_scores_normalized - avg_sea_like_TC_values_normalized + cycle_scores_normalized
        population['fitness'] = targets
        population['avg_sea_like_TC'] = avg_sea_like_TC_values
        
        return population
