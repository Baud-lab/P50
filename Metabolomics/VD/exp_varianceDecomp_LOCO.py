import os
import sys
sys.path.insert(0,'/homes/abaud/P50_HSrats/code/LIMIX')
from social_data_P50_HSrats import SocialData
sys.path.insert(0,'/homes/abaud/CFW/code/reproduce/LIMIX')
from dirIndirVD_wSS import DirIndirVD
import pdb
import h5py
import re
import gc

if __name__=='__main__':
    
    task = 'metabo_counts_notqned_65792'

    col=int(sys.argv[1])-1

    #will be pruned_dosages
    kinship_type = sys.argv[2]
    
    DGE = "DGE"
    IGE = None
    IEE = None
    cageEffect = None #change None to "cageEffect" to fit a LMM with random DGE and random cage effects
    
    if len(sys.argv) == 4:
        subset = sys.argv[3]
    else:
        subset = None #None for P50 rats

    
    print("Subset:")
    print(subset)

    for chr in range(1,21) :

        print('chr is ' + str(chr))
        data = SocialData(task,kinship_type,subset, chr)
        doto = data.get_data(col)
    
        trait=doto['trait']
        print('trait is ' + trait)

        if chr == 1:
            covar_outfile_dir = "".join(["/homes/abaud/P50_HSrats/output/null_covars_LOCO/",task,"/",kinship_type,"_".join(filter(None, ['',subset, DGE, IGE, IEE, cageEffect])),'/',trait])
            if not os.path.exists(covar_outfile_dir):
                os.makedirs(covar_outfile_dir, exist_ok = True)

        covar_outfile_name="".join([covar_outfile_dir,'/',trait,'_chr',str(chr),'.h5'])
   
        vc = DirIndirVD(pheno = doto['pheno'], pheno_ID = doto['pheno_ID'], covs = doto['covs'], covs_ID = doto['covs_ID'], covariates_names = doto['covariates_names'], kinship_all = doto['kinship_full'], kinship_all_ID = doto['kinship_full_ID'], cage_all = doto['cage_full'], cage_all_ID = doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=False, subset_IDs = doto['subset_IDs'], SimplifNonIdableEnvs = False, vc_init_type = None, vc_init = None)

        toSave = vc.getToSave()
        toSave_file = h5py.File(covar_outfile_name,'w')
        toSave_file.create_dataset(name = 'sampleID',data = toSave['sampleID'])
        toSave_file.create_dataset(name = 'pheno',data = toSave['pheno'])
        toSave_file.create_dataset(name = 'covs',data = toSave['covs'])
        toSave_file.create_dataset(name = 'covar_mat',data = toSave['covar_mat'])
        toSave_file.close()

        gc.collect()


