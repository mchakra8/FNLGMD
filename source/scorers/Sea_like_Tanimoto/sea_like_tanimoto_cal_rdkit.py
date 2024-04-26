
import sys,os
from scorers.Sea_like_Tanimoto import fp_rdkit_lib as rfpgen
from scorers.Sea_like_Tanimoto import tanimoto_tversky_cal_axon_lib as tccalc 

## Writen by Trent Balius in the Shoichet Group
## calculates the tanamoto matrics. 
## fingerprint are calculated with chemaxon
## this uses a simular chemaxon comand as sea.  
## bit comparisons are calculated in python

def sea_like_calc(smiles_string,smilesfile2,threshold):

    fpvec1 = rfpgen.get_fp_single(smiles_string)
    fpvec2 = rfpgen.get_fp_from_file(smilesfile2)
    
    sumTC = 0.0
    count = 0
    Mcount = 0

    for fp1 in fpvec1:
        flag_frist = True 
        for fp2 in fpvec2:
            TC = tccalc.tanimoto(fp1,fp2)
            if (flag_frist):
                flag_frist = False
            else:
                pass
            if TC >= threshold:
                sumTC = sumTC + TC
                count = count + 1
            Mcount = Mcount + 1
  #print("sumTC = %f ; count = %d ; avgTC = %f"%(sumTC,count,sumTC/float(count)))
    return sumTC,count,Mcount

