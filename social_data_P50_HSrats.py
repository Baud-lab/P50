import numpy as np
import h5py
#re for regular expressions
import re
import pdb

class SocialData():
    
    def __init__(self, task = None, kinship_type = "", subset = None, chr = None):
        assert task is not None, 'Specify task!'
        self.task=task
        self.kinship_type=kinship_type
        self.subset = subset
        self.chr = chr
        self.load()
    
    
    def load(self):
        
        in_file = '/homes/abaud/P50_HSrats/data/P50_rats_round8.h5'

        f = h5py.File(in_file,'r')

        self.measures = f[self.task]['col_header']['phenotype_ID'].asstr()[:]
        self.all_pheno = f[self.task]['matrix'][:].T
        self.pheno_ID = f[self.task]['row_header']['sample_ID'][:]
        self.all_covs = None
        self.covs_ID = None
        self.covariates = None
        self.all_covs2use = None

        self.cage_full = f[self.task]['row_header']['cage'].asstr()[:]
        self.cage_full_ID = f[self.task]['row_header']['sample_ID'][:]
        print('theres cage info')
        
        if self.chr is not None:
            print('social data in LOCO')
            self.kinship_full = f['GRMs_LOCO'][self.kinship_type][''.join(['chr',str(self.chr)])]['matrix'][:]
            self.kinship_full_ID = f['GRMs_LOCO'][self.kinship_type][''.join(['chr',str(self.chr)])]['row_header']['sample_ID'][:]
        else:    
            self.kinship_full = f['GRM'][self.kinship_type]['matrix'][:]
            self.kinship_full_ID = f['GRM'][self.kinship_type]['row_header']['sample_ID'][:]
 
        if self.subset is None:
            self.subset_IDs = self.kinship_full_ID
        else:
            self.subset_IDs = f['subsets'][self.subset][:]

        
    def get_data(self,col,col_MT = None):
        
        self.trait = self.measures[col]

        self.pheno = self.all_pheno[:,col]
        self.trait_MT = None
        self.track_trait = None

        if self.all_covs2use is None:
            self.covs = None
            self.covs_ID = None
            covariates_names = None
            print("Covariates in social_data: ")
            print(covariates_names)
        else:
            covs2use = self.all_covs2use[col].split(',')
            Ic = np.zeros(self.covariates.shape[0],dtype=bool)
            for cov in covs2use:
                Ic = np.logical_or(Ic,self.covariates==cov)
            covariates_names = self.covariates[Ic]
            print("Covariates in social_data: " + covariates_names)
            self.covs = self.all_covs[:,Ic]
    
        return {'trait' : self.trait,
                'trait_MT' : self.trait_MT,
                'pheno' : self.pheno,
                'pheno_ID' : self.pheno_ID,
                'track_trait' : self.track_trait,
                'covs' : self.covs,
                'covs_ID' : self.covs_ID,
                'covariates_names' : covariates_names,
                'kinship_type' : self.kinship_type,
                'kinship_full' : self.kinship_full,
                'kinship_full_ID' : self.kinship_full_ID,
                'cage_full' : self.cage_full,
                'cage_full_ID' : self.cage_full_ID,
                'subset_IDs' : self.subset_IDs}






