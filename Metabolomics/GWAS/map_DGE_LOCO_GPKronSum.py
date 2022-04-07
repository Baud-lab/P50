#mem usage 20GB

import time
start_time = time.time()

import sys
import os
import h5py
import numpy as np
sys.path.insert(0,'/users/abaud/abaud/micturition/code')
from social_data import SocialData
sys.path.insert(0,'/users/abaud/abaud/LIMIX')
# runs GWAS without estimating the covariance structure: uses the total covariance matrix provided (built from DGE SGE DEE SEE and cage effects)
import pdb
import copy
import gc

from limix_core.gp import GP2KronSum
from limix_core.covar import FreeFormCov
from limix_lmm import LMM

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

task = 'phenotypes'

kinship_type = sys.argv[2]

if len(sys.argv) == 4:
	subset = sys.argv[3]
else:
	subset = None

#chr below (1) could be any - only used to get trait
data = SocialData(task,kinship_type,subset,1)
col=int(sys.argv[1])-1
doto = data.get_data(col)
print("Different covariates actually used")
trait = doto['trait']
print(trait)

DGE = "DGE"
IGE = None
IEE = None
cageEffect = None

#try opening pvalues file early so that lmm doesnt run if file opening is going to fail...
pvalues_file_dir = "".join(['/users/abaud/abaud/micturition/output/pvalues_LOCO/',task,"/",kinship_type,"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])

if not os.path.exists(pvalues_file_dir):
	os.makedirs(pvalues_file_dir, exist_ok=True)
pvalues_file_name = "".join([pvalues_file_dir,'/',trait,'.h5'])
#if os.path.exists(pvalues_file_name):
#	sys.exit(0)

covar_outfile_dir = "".join(["/users/abaud/abaud/micturition/output/null_covars_LOCO/",task,"/",kinship_type,"_".join(filter(None, ['',subset, DGE, IGE, IEE, cageEffect])),'/',trait])
chr = 1
covar_outfile_name = "".join([covar_outfile_dir,"/",trait,"_chr",str(chr),".h5"])
covar_outfile = h5py.File(covar_outfile_name,'r')
saved = {}
saved['sampleID'] = covar_outfile['sampleID'][:]
saved['pheno'] = covar_outfile['pheno'][:]
saved['covs'] = covar_outfile['covs'][:]
saved['covar_mat'] = covar_outfile['covar_mat'][:]
covar_outfile.close()

#import genotypes
input_file_name = '/users/abaud/abaud/micturition/data/datas.h5'
input_file = h5py.File(input_file_name,'r')
geno  = input_file['direct']
geno_matrix = geno['matrix'][:].T
geno_sample_ID = geno['row_header']['sample_ID'][:]
position = {
	"chr" : geno['col_header']['chr'][:],
	"pos"   : geno['col_header']['pos'][:],
	"cumpos" : geno['col_header']['cumpos'][:]
}
input_file.close()

#match genotypes with the rest
Imatch = np.nonzero(saved['sampleID'][:,np.newaxis]==geno_sample_ID)
print("Number of individuals in sampleID and with genotypes: " + str(len(Imatch[0])))
saved['sampleID'] = saved['sampleID'][Imatch[0]]
saved['pheno'] = saved['pheno'][Imatch[0]]
saved['covs'] = saved['covs'][Imatch[0],:]
saved['covar_mat'] = saved['covar_mat'][Imatch[0],:][:,Imatch[0]]
geno_matrix = geno_matrix[Imatch[1],:]
saved_geno_matrix = copy.copy(geno_matrix)

pvalues_file = h5py.File(pvalues_file_name,'w')
pvalues_file.create_dataset(name = 'chr',data = position['chr'])
pvalues_file.create_dataset(name = 'pos',data = position['pos'])
pvalues_file.create_dataset(name = 'cumpos',data = position['cumpos'])

#will define one LMM per chr as covar matrix changes for each chr (LOCO)
for chr in range(1, 20):
    print('chr is ' + str(chr))
    
    if chr != 1:
        covar_outfile_name = "".join([covar_outfile_dir,"/",trait,"_chr",str(chr),".h5"])
        covar_outfile = h5py.File(covar_outfile_name,'r')
        saved['covar_mat'] = covar_outfile['covar_mat'][:][Imatch[0],:][:,Imatch[0]]
        covar_outfile.close()
    
    K = saved['covar_mat']
    Cg = FreeFormCov(1)
    Cn = FreeFormCov(1)
    A = np.eye(1)
 
    gp = GP2KronSum(Y=saved['pheno'], F=saved['covs'], A=A, Cg=Cg, Cn=Cn, R=K)    
    gp.covar.Cr.setCovariance(0.5 * np.ones((1, 1)))
    gp.covar.Cn.setCovariance(0.000001 * np.ones((1, 1)))
    info_opt = gp.optimize(verbose=False)

    lmm = LMM(saved['pheno'], saved['covs'], gp.covar.solve)

    print('additional noise VC is ')
    print(Cn.K())

    geno_matrix = saved_geno_matrix[:,position['chr']==chr]
    lmm.process(geno_matrix)
    pvalues = lmm.getPv()
    #beta = lmm.getBetaSNP()
    pvalues_file.create_dataset(name = "".join(['pvalues_chr',str(chr)]),data = pvalues)
    gc.collect()

pvalues_file.close()

print("My program took", time.time() - start_time, "to run")

