import numpy as np 
import pickle
import sys 
import os
from functions.targetpairs import *

def seqs_probs_total(folddict,pheno0,pheno1):    
	seqsprobs0 =  np.array(folddict[pheno0])
	seqs0 = list(seqsprobs0[:,0])
	probs0 = list(np.float_(list(seqsprobs0[:,1])))

	seqsprobs1 = np.array(folddict[pheno1])
	seqs1 = list(seqsprobs1[:,0])
	probs1 = list(np.float_(list(seqsprobs1[:,1])))
                 
	return seqs0,seqs1,probs0,probs1
 
def seqs_probs_common(seqs0,seqs1,probs0,probs1): #plastic genotypes
	seqscommon = list(set(seqs0) & set(seqs1))
	indices_common_0 = [seqs0.index(x) for x in seqscommon]
	probs_common_0 = [probs0[ind0] for ind0 in indices_common_0]
	indices_common_1 = [seqs1.index(x) for x in seqscommon]
	probs_common_1 = [probs1[ind1] for ind1 in indices_common_1]
	return seqscommon,probs_common_0,probs_common_1,indices_common_0,indices_common_1

def seqs_probs_true_extremes(seqs0,seqs1,probs0,probs1,seqscommon,maxextreme=True):
	seqsextreme0 = list(set(seqs0) - set(seqscommon))
	seqsextreme1 = list(set(seqs1) - set(seqscommon))

	return seqsextreme0, seqsextreme1

def seqs_probs_soft_extremes(probs_common_0,probs_common_1,seqscommon): #the low plasticity phenotypes
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



