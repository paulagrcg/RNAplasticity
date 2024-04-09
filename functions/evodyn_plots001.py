import numpy as np
from collections import defaultdict, Counter
import random
import pickle
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import scipy
from scipy import stats

stb = {'A': '0', 'U': '1', 'C': '2', 'G': '3'}

def stringtoint(seqs):
	seqsint = []
	for seq in seqs:
		seqsint.append(int(''.join([stb[n] for n in seq])))
	return seqsint

def countseqtypes(generations, popgenos, seqscommon, extremes0, extremes1, seqstarget):
	countseqscommon = [0]*generations
	countseqstarget = [0]*generations
	countseqsex0 = [0]*generations
	countseqsex1 = [0]*generations
	countother = [0]*generations
	seqsex0Counter = Counter(extremes0)
	seqsex1Counter = Counter(extremes1)
	seqscommonCounter = Counter(seqscommon)
	seqstargetCounter = Counter(seqstarget)

	for gen in range(1,generations+1):
		freq = Counter(popgenos[gen])
		for seq in freq.keys():
			n = seqstargetCounter[seq]
			m = seqscommonCounter[seq]
			k = seqsex0Counter[seq]
			l = seqsex1Counter[seq]
			if n == 1: countseqstarget[gen-1] += freq[seq]
			elif m == 1: countseqscommon[gen-1] += freq[seq]
			elif k == 1: countseqsex0[gen-1] += freq[seq]
			elif l == 1: countseqsex1[gen-1] += freq[seq]
			else: countother[gen-1] += freq[seq]
	return countseqstarget,countseqscommon,countseqsex0,countseqsex1,countother

def countseqtypes1(popgenos, seqscommon, extremes0, extremes1, seqstarget):
	countseqscommon = 0
	countseqstarget = 0
	countseqsex0 = 0
	countseqsex1 = 0
	countother = 0
	seqsex0Counter = Counter(extremes0)
	seqsex1Counter = Counter(extremes1)
	seqscommonCounter = Counter(seqscommon)
	seqstargetCounter = Counter(seqstarget)
	freq = Counter(popgenos)
	for seq in freq.keys():
		n = seqstargetCounter[seq]
		m = seqscommonCounter[seq]
		k = seqsex0Counter[seq]
		l = seqsex1Counter[seq]
		if n == 1: countseqstarget += freq[seq]
		elif m == 1: countseqscommon += freq[seq]
		elif k == 1: countseqsex0 += freq[seq]
		elif l == 1: countseqsex1 += freq[seq]
		else: countother += freq[seq]
	return countseqstarget,countseqscommon,countseqsex0,countseqsex1,countother

def extremes_with_flip(extremeoption,generations, gengap, countseqsex0,countseqsex1):
	extremes = []
	now = 1-extremeoption 
	for gen in range(1,generations+1): 
		if gen % gengap == 0:
			now = 1 - now
		if now == 0: extremes.append(countseqsex0[gen-1])
		if now == 1: extremes.append(countseqsex1[gen-1])
	return extremes
 	
if __name__ == "__main__":
	
	parameters = {}
	pair = int(sys.argv[1])
	fitlands = int(sys.argv[2])#1.random 2.hamming 3.random with inverse
	samplenum = int(sys.argv[3])
	plasticoption = int(sys.argv[4])
	gengap = int(sys.argv[5])
	generations = int(sys.argv[6])
	mu = float(sys.argv[7])
	minprob = float(sys.argv[8])/100.
	maxprob = float(sys.argv[9])/100.
	initoption = int(sys.argv[10]) #extremeoption -> 0: soft extremes, 1: true extremes, 2: random from seqscommon
	extremeoption = int(sys.argv[11]) #e.g. if we want to start simulation with target 0, the initial population should be at extreme of target 1
	Npop = int(sys.argv[12])
	start = int(sys.argv[13])
	end = int(sys.argv[14])
	parametersweep = int(sys.argv[15])
	run = int(sys.argv[16])
	if parametersweep == 1:
		path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"
		pathh = path + "fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_plasticoption"+str(plasticoption)+"_generations"+str(generations)+"_Npop"+str(Npop)+ "/"
		path1 = pathh + "run_001"+str(run) + "/"
		df = pd.read_csv('/rds/user/pg520/hpc-work/targetflipping/parametersweep_0.001.csv')
		mu = float(df.iloc[run]['mu'])
		samplenum = int(df.iloc[run]['samplenum'])
		gengap = int(df.iloc[run]['gengap'])
	else: 
		path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"
		path1 = path + "fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+ "/"

	seqstrueextreme0 = stringtoint(np.load(path + "seqstrueextreme0.npy"))
	seqstrueextreme1 = stringtoint(np.load(path + "seqstrueextreme1.npy"))
	seqscommon = stringtoint(np.load(path + "seqscommon.npy"))
	seqstarget = stringtoint(np.load(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy"))
	
	with open(path1 + "popgenos.pkl", "rb") as f: 
		popgenos = pickle.load(f)
	df = pd.read_csv(path1 + 'data_fit_tgs_ev_entr_robust.txt', sep=" ", header=None)
	df.columns = ["gen", "fitness","TGs","evolv","entroy","robust"]
	
	countseqstarget,countseqscommon,countseqsex0,countseqsex1,countother = countseqtypes(generations, popgenos, seqscommon, seqstrueextreme0, seqstrueextreme1, seqstarget)

	extremes = extremes_with_flip(extremeoption,generations, gengap, countseqsex0,countseqsex1)

	############################ PLOTS ####################################

	
	f = plt.figure()
	ax = plt.gca()
	f.set_figwidth(40)
	f.set_figheight(10)

	plt.plot(np.array(df["gen"])[start:end], countseqsex0[start:end], label = "seqs true extremes 0")
	plt.plot(np.array(df["gen"])[start:end], countseqsex1[start:end], label = "seqs true extremes 1")
	plt.plot(np.array(df["gen"])[start:end], countseqscommon[start:end], label = "seqs common (Pearson r: "+str(stats.pearsonr(countseqscommon,df["fitness"])[0])+")")
	plt.plot(np.array(df["gen"])[start:end], countseqstarget[start:end], label = "seqs target [" + str(minprob) + "," + str(maxprob)+ "]. (Pearson r: "+str(stats.pearsonr(countseqstarget,df["fitness"])[0])+")")
	plt.plot(np.array(df["gen"])[start:end], countother[start:end], label = "other seqs (Pearson r: "+str(stats.pearsonr(countother,df["fitness"])[0])+")")
	plt.scatter(np.array(df["gen"])[start:end],extremes[start:end], marker = "+", label = "extremes with flip (Pearson r: "+str(stats.pearsonr(extremes,df["fitness"])[0])+")", color = "k")
	ax.set_ylabel("N")
	ax2 = ax.twinx()
	color = 'm'
	ax2.plot(np.array(df["gen"])[start:end], np.array(df["fitness"])[start:end], label= "fitness", color = color)
	ax.set_title(path1) 
	ax2.set_ylabel("Fitness", color = color)
	ax2.tick_params(axis='y',labelcolor = color)
	ax.legend(loc='upper center',  ncol = 3)	
	f.tight_layout()
	plt.grid()
	plt.savefig(path1+"plot_genos_vs_fitness_"+str(start)+"_"+str(end)+".eps")	
