import numpy as np
import pickle 
from gpmapfunctions import *
import sys
from collections import defaultdict
import random
from seqsprobstf import *

def fitness_hamming(fold,target):
    hamdist = 0
    for j in range(0,len(fold)):
        if fold[j]!=target[j]: hamdist += 1/(len(fold))
        fitness = 1. - hamdist
    return fitness

def plastic_fitness_hamming(target,seq,gpmap,L):
   phsprobs =extractnormalisedprobs(gpmap[seq],L)
   folds = list(phsprobs.keys())
   probs = list(phsprobs.values())
   plasticfitness = 0
   for f,p in zip(folds,probs):
       plasticfitness += fitness_hamming(f,target)*p
   return plasticfitness

def fitnessgeno(geno,gpmap,foldfitdict,fitlands,pair,samplenum,targetpheno,targetpheno_old,L=12,plastic=False):
	if plastic ==True: fitnessplastic = []
	else: fitnessmeansamples = []
	phvsprobseq = extractnormalisedprobs(gpmap[geno],L)
	folds = list(phvsprobseq.keys())
	probs = list(phvsprobseq.values())
	for f in folds: #fitness for each fold in the geno ensemble
		if f not in foldfitdict.keys(): #we assume that we are calling this function with consistent fitlands option to fill up the dictionary
			if fitlands ==1: #random fitness
				if f == targetpheno: foldfitdict[targetpheno] = 1.0
				elif f == '.'*L: foldfitdict['.'*L] = 0.0
				else: foldfitdict[f] = np.random.random_sample() #uniform between 0 and 1
			if fitlands == 2: #hamming fitness
				hamdist = 0
				if f == '.'*L: foldfitdict['.'*L] = 0.0
				for j in range(0,len(f)):
					if f[j]!=targetpheno[j]: hamdist += 1/(len(f))
				foldfitdict[f] = 1. - hamdist #fitness is 1 - hamdist since as structure is more different than TP then lower fitness
			if fitlands == 3:#hamming fitness with power 2 
				hamdist = 0
				if f == '.'*L: foldfitdict['.'*L] = 0.0
				for j in range(0,len(f)):
					if f[j]!=targetpheno[j]: hamdist += 1/(len(f))
				foldfitdict[f] = (1. - hamdist)**2
			if fitlands == 4:#random fitness with targets inverse max/min
				if f == targetpheno: foldfitdict[targetpheno] = 1.0
				elif f == '.'*L: foldfitdict['.'*L] = 0.0
				elif f == targetpheno_old: foldfitdict[targetpheno_old] = 0.0
				else: foldfitdict[f] = np.random.random_sample()
	if plastic == True: 
		plasticfitness = 0
		for f,p in zip(folds,probs): 
			plasticfitness += foldfitdict[f]*p	
		return foldfitdict, plasticfitness #return the updated foldfitdict and the fitness of that geno
	else: 
		maxfitlist = []
		for i in range(0,100):
			foldschoices = random.choices(folds,weights=probs,k=samplenum)
			maxfitness = max([foldfitdict[k] for k in foldschoices])
			maxfitlist.append(maxfitness)
		meanmax = np.mean(maxfitlist)
		return foldfitdict, meanmax #return the updated foldfitdict and the fitness of that geno


if __name__ == "__main__":
	
	pair = int(sys.argv[1])
	fitlands = int(sys.argv[2])
	samplenum = int(sys.argv[3])
	plasticoption = int(sys.argv[4])
	minprob = float(int(sys.argv[5])/10.)
	maxprob = float(int(sys.argv[6])/10.)
	#with open("/rds/user/pg520/hpc-work/folddictt.pkl","rb") as f:
	#     folddict = pickle.load(f)
	with open("/home/pg520/target_flipping/phenopairs.pkl","rb") as f:phenopairs = pickle.load(f)

	path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"
	if not os.path.isdir(path): os.mkdir(path)
	phenopair = phenopairs[pair]
	pheno0 = phenopair[0]
	pheno1 = phenopair[1]
	with open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl","rb") as f:
		gpmap = pickle.load(f)
	
	seqscommon = np.load(path + "seqscommon.npy")
	probs_common_0 = np.load(path + "probs_common_0.npy")
	probs_common_1 = np.load(path + "probs_common_1.npy")
	
	seqstarget = np.load(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy")

	seq_soft_extreme_0 = np.load(path + "seqssoftextreme0.npy")
	seq_soft_extreme_1 = np.load(path + "seqssoftextreme1.npy")

	#seq_extreme_0_max = np.load(path + "/seq_trueextreme_0_max.npy")
	#seq_extreme_1_max = np.load(path + "/seq_trueextreme_1_max.npy")
	#probs_extreme_0_max = np.load(path + "/probs_trueextreme_0_max.npy")
	#probs_extreme_1_max = np.load(path + "/probs_trueextreme_1_max.npy")

	L=12
	fitnessgenos = defaultdict(list)
	fitnesstargets = defaultdict(list)
	fitsoftextreme = []
	fitnessmaxextreme = defaultdict(float)
	for targetphenonum in range(0,2):
		foldfitdict = defaultdict(float)
		if targetphenonum == 0: targetphenonum_old = 1
		if targetphenonum == 1: targetphenonum_old = 0
		targetpheno = phenopairs[pair][targetphenonum]
		targetpheno_old = phenopairs[pair][targetphenonum_old]
		
		print("targetseqs: ", seqstarget)
		print("targetseqs fraction:",float(len(seqstarget)/len(probs_common_0)))
		
		if targetphenonum == 1: 
			seqsoftextreme = seq_soft_extreme_0 #extremeoption=0: #p0>p1
			#if maxextreme =True: seqprobsextreme = (seq_extreme_0_max,probs_extreme_0_max, 0.0)
		elif targetphenonum == 0: 
			seqsoftextreme = seq_soft_extreme_1  #extremeoption=1:#p1>p0 extreme
			#if maxextreme = True: seqprobsextreme = (seq_extreme_1_max,0.0,probs_extreme_1_max)
			
		for g,p0,p1 in zip(seqscommon,probs_common_0,probs_common_1):
			foldfitdict, fitness = fitnessgeno(g,gpmap,foldfitdict,fitlands,pair,samplenum,targetpheno,targetpheno_old,L=12,plastic=plasticoption)	
			if g in seqstarget: 
				fitnesstargets[targetphenonum].append((p0,p1,fitness))
			if str(g) == str(seqsoftextreme): 
				fitsoftextreme.append((p0,p1,fitness))
				print(p0,p1,fitness)
			fitnessgenos[targetphenonum].append(fitness)
		#foldfitdict, fitness = fitnessgeno(seqprobsextreme[0],gpmap,foldfitdict,fitlands,pair,targetpheno,targetpheno_old,L=12,plastic=plasticoption)
		#fitnessextreme[targetphenonum].append((seqprobsextreme[0],seqprobsextreme[1],fitness))
	np.savetxt(path + "fitness_plastic"+str(plasticoption)+"_fitlands"+str(fitlands)+"_samplenum"+str(samplenum)+".txt",np.c_[np.array(probs_common_0),np.array(probs_common_1),np.array(fitnessgenos[0]),np.array(fitnessgenos[1])])
	np.savetxt(path + "fitness_targets_plastic"+str(plasticoption)+"_fitlands"+str(fitlands)+"_samplenum"+str(samplenum)+".txt",np.c_[np.array(list(list(zip(*fitnesstargets[0]))[0])),np.array(list(list(zip(*fitnesstargets[0]))[1])),np.array(list(list(zip(*fitnesstargets[0]))[2])),np.array(list(list(zip(*fitnesstargets[1]))[2]))])
	np.savetxt(path + "fitness_softextreme_plastic"+str(plasticoption)+"_fitlands"+str(fitlands)+"_samplenum"+str(samplenum)+".txt",np.c_[np.array(list(list(zip(*fitsoftextreme))[0])),np.array(list(list(zip(*fitsoftextreme))[1])),np.array(list(list(zip(*fitsoftextreme))[2]))])

