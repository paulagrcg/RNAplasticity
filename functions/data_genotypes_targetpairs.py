import numpy as np 
import pickle
import sys 
import os
import numpy as np 
import pickle
import sys 
import os
from functions.targetpairs import *
from functions.genotype_populations import *


if __name__ == "__main__": 
	seed = int(sys.argv[1])
	minprob = int(sys.argv[2])
	maxprob = int(sys.argv[3])
	minprob = float(minprob/100.)
	maxprob = float(maxprob/100.)

	with open("/rds/user/pg520/hpc-work/folddictt.pkl","rb") as f: folddict = pickle.load(f)

	with open("/rds/user/pg520/hpc-work/NDsetsize.pkl","rb") as f:
		neutralsets= pickle.load(f)
	
	complexcomplex, complexsimple, simplesimple = phenosizescale(neutralsets, seed)
	cc, cc_dist = hammingdistance_sets(complexcomplex)
	cs, cs_dist  = hammingdistance_sets(complexsimple)
	ss, ss_dist = hammingdistance_sets(simplesimple)
	path = "/rds/user/pg520/hpc-work/targetflipping/seed_"+str(seed)+"/"

	np.save(path + "cc_cs_ss_pairs_seed_"+str(seed)+".npy", np.array([cc,cs,ss]))
	np.save(path + "cc_cs_ss_distances_seed_"+str(seed)+".npy", np.array([cc_dist,cs_dist,ss_dist]))
	
	for paircc in range(0,len(cc)):
    
		pheno0 = cc[paircc][0]
		pheno1 = cc[paircc][1]
		pathcc = path + "cc_"+str(paircc)+"/"
		seqs0,seqs1,probs0,probs1 = seqs_probs_total(folddict,pheno0,pheno1)
		
		seqscommon,probs_common_0,probs_common_1,indices_common_0,indices_common_1 = seqs_probs_common(seqs0,seqs1,probs0,probs1)
		np.save(pathcc + "seqscommon.npy",seqscommon)

		
		seqstarget,probs_target_0,probs_target_1 = seqs_probs_targetgenos(minprob,maxprob,probs_common_0,probs_common_1,seqscommon)
		np.save(pathcc + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy",seqstarget)
		
		seq_soft_extreme_0,seq_soft_extreme_1,prob_0_softextreme_0,prob_1_softextreme_0,prob_0_softextreme_1,prob_1_softextreme_1 = seqs_probs_soft_extremes(probs_common_0,probs_common_1,seqscommon)
		np.save(pathcc + "seqssoftextreme0.npy",seq_soft_extreme_0)
		np.save(pathcc + "seqssoftextreme1.npy",seq_soft_extreme_1)

		
		seqsextreme0, seqsextreme1 = seqs_probs_true_extremes(seqs0,seqs1,probs0,probs1,seqscommon,maxextreme=True)
		np.save(pathcc + "/seqstrueextreme1.npy",seqsextreme1)
		np.save(pathcc + "/seqstrueextreme0.npy",seqsextreme0)



	for paircs in range(0,len(cs)):
        
		pheno0 = cs[paircs][0]
		pheno1 = cs[paircs][1]
		pathcs = path + "cs_"+str(paircs)+"/"
		seqs0,seqs1,probs0,probs1 = seqs_probs_total(folddict,pheno0,pheno1)
		
		seqscommon,probs_common_0,probs_common_1,indices_common_0,indices_common_1 = seqs_probs_common(seqs0,seqs1,probs0,probs1)
		np.save(pathcs + "seqscommon.npy",seqscommon)

		seqstarget,probs_target_0,probs_target_1 = seqs_probs_targetgenos(minprob,maxprob,probs_common_0,probs_common_1,seqscommon)
		np.save(pathcs + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy",seqstarget)

		seq_soft_extreme_0,seq_soft_extreme_1,prob_0_softextreme_0,prob_1_softextreme_0,prob_0_softextreme_1,prob_1_softextreme_1 = seqs_probs_soft_extremes(probs_common_0,probs_common_1,seqscommon)
		np.save(pathcs + "seqssoftextreme0.npy",seq_soft_extreme_0)
		np.save(pathcs + "seqssoftextreme1.npy",seq_soft_extreme_1)

		seqsextreme0, seqsextreme1 = seqs_probs_true_extremes(seqs0,seqs1,probs0,probs1,seqscommon,maxextreme=True)
		np.save(pathcs + "/seqstrueextreme1.npy",seqsextreme1)
		np.save(pathcs + "/seqstrueextreme0.npy",seqsextreme0)


	for pairss in range(0,len(ss)):
        
		pheno0 = ss[pairss][0]
		pheno1 = ss[pairss][1]
		pathss = path + "ss_"+str(pairss)+"/"
		seqs0,seqs1,probs0,probs1 = seqs_probs_total(folddict,pheno0,pheno1)
		
		seqscommon,probs_common_0,probs_common_1,indices_common_0,indices_common_1 = seqs_probs_common(seqs0,seqs1,probs0,probs1)
		np.save(pathss + "seqscommon.npy",seqscommon)

		seqstarget,probs_target_0,probs_target_1 = seqs_probs_targetgenos(minprob,maxprob,probs_common_0,probs_common_1,seqscommon)
		np.save(pathss + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy",seqstarget)

		seq_soft_extreme_0,seq_soft_extreme_1,prob_0_softextreme_0,prob_1_softextreme_0,prob_0_softextreme_1,prob_1_softextreme_1 = seqs_probs_soft_extremes(probs_common_0,probs_common_1,seqscommon)
		np.save(pathss + "seqssoftextreme0.npy",seq_soft_extreme_0)
		np.save(pathss + "seqssoftextreme1.npy",seq_soft_extreme_1)

		
		seqsextreme0, seqsextreme1 = seqs_probs_true_extremes(seqs0,seqs1,probs0,probs1,seqscommon,maxextreme=True)
		np.save(pathss + "/seqstrueextreme1.npy",seqsextreme1)
		np.save(pathss + "/seqstrueextreme0.npy",seqsextreme0)
