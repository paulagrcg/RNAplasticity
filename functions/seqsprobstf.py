import numpy as np 
import pickle
import sys 
import os
#Collect all seqs and probs data for fitness landscape and target flipping evodyn
def seqs_probs_total(folddict,pheno0,pheno1):    
	seqsprobs0 =  np.array(folddict[pheno0])
	seqs0 = list(seqsprobs0[:,0])
	probs0 = list(np.float_(list(seqsprobs0[:,1])))

	seqsprobs1 = np.array(folddict[pheno1])
	seqs1 = list(seqsprobs1[:,0])
	probs1 = list(np.float_(list(seqsprobs1[:,1])))
                 
	return seqs0,seqs1,probs0,probs1

def seqs_probs_common(seqs0,seqs1,probs0,probs1):
	seqscommon = list(set(seqs0) & set(seqs1))
	indices_common_0 = [seqs0.index(x) for x in seqscommon]
	probs_common_0 = [probs0[ind0] for ind0 in indices_common_0]
	indices_common_1 = [seqs1.index(x) for x in seqscommon]
	probs_common_1 = [probs1[ind1] for ind1 in indices_common_1]
	return seqscommon,probs_common_0,probs_common_1,indices_common_0,indices_common_1

def seqs_probs_true_extremes(seqs0,seqs1,probs0,probs1,seqscommon,maxextreme=True):
	seqsextreme0 = list(set(seqs0) - set(seqscommon))
	seqsextreme1 = list(set(seqs1) - set(seqscommon))
	#indices_extreme_0 = [seqs0.index(x) for x in seqsextreme0]
	#probs_extreme_0 = [probs0[ind0] for ind0 in indices_extreme_0]
	#indices_extreme_1 = [seqs1.index(x) for x in seqsextreme1]
	#probs_extreme_1 = [probs1[ind1] for ind1 in indices_common_1]
	#indexmaxextreme0 = np.argmax(np.array(probs_extreme_0))
	#indexmaxextreme1 = np.argmax(np.array(probs_extreme_1))
	#seq_extreme_0_max = seqsextreme0[indexmaxextreme0]
	#seq_extreme_1_max = seqsextreme1[indexmaxextreme1]
	#probs_extreme_0_max = probs_extreme_0[indexmaxextreme0]
	#probs_extreme_1_max = probs_extreme_1[indexmaxextreme1]
	#return seqsextreme0,seqsextreme1,probs_extreme_0,probs_extreme_1,seq_extreme_0_max,seq_extreme_1_max,probs_extreme_0_max,probs_extreme_1_max
	return seqsextreme0, seqsextreme1

def seqs_probs_soft_extremes(probs_common_0,probs_common_1,seqscommon):
	p0 = np.array(probs_common_0)
	p1 = np.array(probs_common_1)
	diff = p0-p1
	indmaxp0 = diff.argmax()#extreme soft for pheno0
	indmaxp1 = diff.argmin()#extreme soft for pheno1
	seq_soft_extreme_0 = seqscommon[indmaxp0]
	seq_soft_extreme_1 = seqscommon[indmaxp1]
	prob_0_softextreme_0 = probs_common_0[indmaxp0]
	prob_1_softextreme_0 = probs_common_1[indmaxp0]
	prob_0_softextreme_1 = probs_common_0[indmaxp1]
	prob_1_softextreme_1 = probs_common_1[indmaxp1]
	return seq_soft_extreme_0,seq_soft_extreme_1,prob_0_softextreme_0,prob_1_softextreme_0,prob_0_softextreme_1,prob_1_softextreme_1

def seqs_probs_targetgenos(minprob,maxprob,probs_common_0,probs_common_1,seqscommon):
	listtupple_seq_p0_p1 = [(seq,p0,p1) for seq,p0,p1 in zip(seqscommon,probs_common_0,probs_common_1) if  p0 < maxprob and p0 > minprob and p1 < maxprob  and p1 > minprob]
	seqstarget =list(list(zip(*listtupple_seq_p0_p1))[0])
	probs_target_0 = list(list(zip(*listtupple_seq_p0_p1))[1])
	probs_target_1 = list(list(zip(*listtupple_seq_p0_p1))[2])
	return seqstarget,probs_target_0,probs_target_1

if __name__ == "__main__": 

	pair = int(sys.argv[1])
	minprob = int(sys.argv[2])
	maxprob = int(sys.argv[3])
	minprob = float(minprob/100.)
	maxprob = float(maxprob/100.)
	with open("/rds/user/pg520/hpc-work/folddictt.pkl","rb") as f: folddict = pickle.load(f)
	with open("/home/pg520/target_flipping/phenopairs.pkl","rb") as f: phenopairs = pickle.load(f)
	path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"
	if not os.path.isdir(path):os.mkdir(path)	
	phenopair = phenopairs[pair]
	
	pheno0 = phenopair[0]
	pheno1 = phenopair[1]
	del phenopair
	seqs0,seqs1,probs0,probs1 = seqs_probs_total(folddict,pheno0,pheno1)
	del folddict
	seqscommon,probs_common_0,probs_common_1,indices_common_0,indices_common_1 = seqs_probs_common(seqs0,seqs1,probs0,probs1)
	"""
	np.save(path + "seqscommon.npy",seqscommon)
	np.save(path + "probs_common_0.npy",probs_common_0)
	np.save(path + "probs_common_1.npy",probs_common_1)
	"""
	seqstarget,probs_target_0,probs_target_1 = seqs_probs_targetgenos(minprob,maxprob,probs_common_0,probs_common_1,seqscommon)
	np.save(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy",seqstarget)
	np.save(path + "probs_target_0_"+str(minprob)+"_"+str(maxprob)+".npy",probs_target_0)
	np.save(path + "probs_target_1_"+str(minprob)+"_"+str(maxprob)+".npy",probs_target_1)
	"""
	seq_soft_extreme_0,seq_soft_extreme_1,prob_0_softextreme_0,prob_1_softextreme_0,prob_0_softextreme_1,prob_1_softextreme_1 = seqs_probs_soft_extremes(probs_common_0,probs_common_1,seqscommon)
	np.save(path + "seqssoftextreme0.npy",seq_soft_extreme_0)
	np.save(path + "seqssoftextreme1.npy",seq_soft_extreme_1)
	np.save(path + "prob_0_softextreme_0.npy",prob_0_softextreme_0)
	np.save(path + "prob_1_softextreme_0.npy",prob_1_softextreme_0)
	np.save(path + "prob_0_softextreme_1.npy",prob_0_softextreme_1)
	np.save(path + "prob_1_softextreme_1.npy",prob_1_softextreme_1)
	"""
	seqsextreme0, seqsextreme1 = seqs_probs_true_extremes(seqs0,seqs1,probs0,probs1,seqscommon,maxextreme=True)
	np.save(path + "/seqstrueextreme1.npy",seqsextreme1)
	np.save(path + "/seqstrueextreme0.npy",seqsextreme0)
	#np.save(path + "/probs_trueextreme_0.npy",probs_extreme_0)
	#np.save(path + "/probs_trueextreme_1.npy",probs_extreme_1)
	#np.save(path + "/seq_trueextreme_0_max.npy",seq_extreme_0_max)
	#np.save(path + "/probs_trueextreme_0_max.npy",probs_extreme_0_max)
	#np.save(path + "/seq_trueextreme_1_max.npy",seq_extreme_1_max)
	#np.save(path + "/probs_trueextreme_1_max.npy",probs_extreme_1_max)


