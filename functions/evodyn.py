import numpy as np
from collections import defaultdict
import random
from quantitiesND4 import *
from gpmapfunctions import *
from seqsprobstf import *
from fitnesslandscape import fitness_hamming
import pandas as pd
from evodyn_plots001 import stringtoint
from collections import Counter

stn = {'A': '0', 'U': '1', 'C': '2', 'G': '3'}
def evodyn(rhogND,evgND,entropy,targets,initpop, genostargetset, seqscommon,extreme0, extreme1, samplenum = 1,plasticoption = False, gengap = 10, L=12, gpmap = {}, mu = 0.05, generations = 10000, Npop = 100, fitlands = 2):
	popgenos = {}
	#initpop is a tupple with initial target pheno and initial population list
	targetpheno = targets[initpop[0]]

	initpopulation = []
	for i in range(0,Npop): initpopulation.append(str(initpop[1]))
	popgenos[0] = initpopulation

	probRWS = defaultdict(float)
	probstarg0dict = defaultdict(dict)
	probstarg1dict = defaultdict(dict)
	meanfitness = defaultdict(float)
	switchnum = 0
	listseqprobs = {}
	targetgenosgen = {}
	adaptedphenonum = {}
	meanentropygenos = {}
	meanevolvgnd = {}
	meanrobustg = {}
	seqstargetlist = []
	seqscommonlist = []
	extreme0list = []
	extreme1list = []
	otherlist = []
	foldfitdict = defaultdict(float)
	seqscommonCounter = Counter(seqscommon)
	seqstargetCounter = Counter(seqstarget)
	seqsex0Counter = Counter(extreme0)
	seqsex1Counter = Counter(extreme1)

	for i in range(1,generations+1):
		listseqprobs[i] = defaultdict(dict)
		targetgenosgen[i] = 0
		adaptedphenonum[i] = 0
		meanevolvgnd[i] = 0
		meanentropygenos[i] = 0
		meanrobustg[i] = 0
	#adapted regime
	countadapt = 0 

	for i in range(1,generations+1): # i is generation number
		targetgenos = 0
		#analyse how each member of the population changes its prob0 and prob1 respectively
		prob0 = [0]*(Npop) 
		prob1 = [0]*(Npop)

		probs = defaultdict(float)
		if i==1: genos = initpopulation #initalise population
		       
		else: genos = newpop 

		#for geno in genos:
		#   listgenosint = [] 
		#	listgenosint.append(int(''.join([stn[n] for n in geno])))
		popgenos[i] = [int(''.join([stn[n] for n in geno])) for geno in genos] #taken from i-1 mutated population
		#print(popgenos[i])
		#print(seqscommon)
		#types of genos data collection (in frequencies)
		freq = Counter(popgenos[i])
		countseqstarget = sum(freq[seq] for seq in freq if seqstargetCounter[seq] == 1)
		countseqscommon = sum(freq[seq] for seq in freq if seqscommonCounter[seq] == 1)
		countseqsex0 = sum(freq[seq] for seq in freq if seqsex0Counter[seq] == 1)
		countseqsex1 = sum(freq[seq] for seq in freq if seqsex1Counter[seq] == 1)
		countother = sum(freq[seq] for seq in freq if seqstargetCounter[seq] != 1 and seqscommonCounter[seq] != 1 and seqsex0Counter[seq] != 1 and seqsex1Counter[seq] != 1)

		seqstargetlist.append(countseqstarget/Npop)
		seqscommonlist.append(countseqscommon/Npop)
		extreme0list.append(countseqsex0/Npop)
		extreme1list.append(countseqsex1/Npop)
		otherlist.append(countother/Npop)

		folds = {}
		probs = {}
		#data of population phenotype ensembles and TPs probs
		for j in range(0,Npop):
			phvsprobseq = extractnormalisedprobs(gpmap[genos[j]],L)
			listseqprobs[i][j] = phvsprobseq
			for f,p in phvsprobseq.items():
				if f==targets[0][0:L]: prob0[j] = p
				if f== targets[1][0:L]: prob1[j]= p
			folds[j] = list(phvsprobseq.keys())
			probs[j] = list(phvsprobseq.values())

		#data collection of target genotypes and evolvability
		for g in genos:
			if g in genostargetset:targetgenos +=1
			meanentropygenos[i] += entropy[g]/len(genos)
			meanevolvgnd[i] += evgND[g]/len(genos)
			meanrobustg[i] += rhogND[g]/len(genos)
		   #meanevolvgnd[i]+=evolvabilitygND(gpmap,g,4,12)[g]/len(genos)

		targetgenosgen[i] = targetgenos

		#1) random fitness
		#2) hamming fitness 
		#3) random fitness with targets inverse gax/min

		if i%gengap==0: #change target every gengap generations
			switchnum +=1
			if targetpheno == targets[0]:
				targetpheno = targets[1]
				targetpheno_old = targets[0]
			else:
				targetpheno = targets[0]
				targetpheno_old = targets[1]

 
			if fitlands ==1 or fitlands == 3: #random fitness
				for phs in folds.values():
					for f in phs: #don't change all the random values from before just change the targetphenos 
						if f == targetpheno: foldfitdict[targetpheno] = 1.0
						elif f == '.'*L: foldfitdict['.'*L] = 0.0
						elif fitlands == 3 and f == targetpheno_old: foldfitdict[targetpheno_old] = 0.0 #random fitness with targets inverse max/min
						elif f not in foldfitdict.keys(): foldfitdict[f] = np.random.random_sample() #uniform between 0 and 1
			else:#hamming fitness new fitness dict
				foldfitdict = defaultdict(float)
				for phs in folds.values():
					for f in phs:foldfitdict[f] = 0 
				for f in foldfitdict.keys():
					if f == '.'*L: foldfitdict['.'*L] = 0.0
					if fitlands==2: foldfitdict[f] = fitness_hamming(f,targetpheno)
					if fitlands==3: foldfitdict[f] = fitness_hamming_power(f,targetpheno)

		else: #no target flip, just add new phenotypes into fitness landscape 
			for phs in folds.values():
				for f in phs:
					if f not in foldfitdict.keys():
						if fitlands==1: #random fitness 
							if f == targetpheno: foldfitdict[targetpheno] = 1.0
							elif f == '.'*L: foldfitdict['.'*L] = 0.0
							else: foldfitdict[f] = np.random.random_sample()
						if fitlands==2: #hamming fitness
							if f == '.'*L: foldfitdict['.'*L] = 0.0
							else: foldfitdict[f] = fitness_hamming(f,targetpheno)
						if fitlands==3:
							if f== '.'*L: foldfitdict['.'*L] = 0.0
							else: foldfitdict[f] = fitness_hamming_power(f,targetpheno)
						if fitlands == 4:#random fitness with targets inverse max/min
							if f == targetpheno: foldfitdict[targetpheno] = 1.0
							elif f == '.'*L: foldfitdict['.'*L] = 0.0
							elif f == targetpheno_old: foldfitdict[targetpheno_old] = 0.0
						else: foldfitdict[f] = np.random.random_sample()


		#population PLASTIC fitness if samplenum = 1, or MAX fitness if samplenum > 1 
		popgen = {}
		probstarg0 = {}
		probstarg1 = {}
		count = 0
		cadapt = 0
		for ng,phs,ps in zip(genos,folds.values(),probs.values()):
			plasticfitness = 0 
			maxfitness = 0
			#maxfitness out of the sampled plastic phenotypes
			if plasticoption == False: 
				foldschoices = random.choices(list(phs),weights = list(ps), k=samplenum)
				maxfitness = max([foldfitdict[k] for k in foldschoices])
				popgen[(ng,count)] = maxfitness
			#plastic fitness (weighted)
			else:
				for index in range(0,len(phs)): plasticfitness += foldfitdict[phs[index]]*ps[index]
				popgen[(ng,count)] = plasticfitness#population data for plastic fitnesses

			#population data  for TPs probs
			probstarg0[(ng,count)] = prob0[count]
			probstarg1[(ng,count)] = prob1[count]
			count+=1
        
        
		#mean fitness of population at gen i       
		meanfitness[i] = sum(list(popgen.values()))/len(popgen)
		#data collection of target probabilities
		probstarg0dict[i] = probstarg0
		probstarg1dict[i] = probstarg1

		#roulette wheel selection
		newpop = []
		choices = []
		index = []
		fsum = sum(list(popgen.values()))
		probRWS = defaultdict(float)

		for g,ph in popgen.items():
			probRWS[g] = float(ph/fsum)
			choices.append(g[0])

		#sampling of new population 
		newpop = [random.choices(choices, weights=list(probRWS.values()), k=1)[0] for num in range(0,Npop)]
		#for j in range(0,Npop): newpop.append(random.choices(choices, weights=list(probRWS.values()), k=1)[0])  
		     
		#mutations for population i+1
		newpop = [mutation(seq, mu) for seq in newpop]
		#for j in range(0,Npop):
		#	seq = newpop[j]
		#	newpop[j] = mutation(seq,mu) #mutate member j with rate mu per base 
	probstarg0dict[i] = probstarg0
	probstarg1dict[i] = probstarg1
	#return meanrobustg, meanentropygenos, meanfitness,adaptedphenonum,targetgenosgen, probstarg0dict, probstarg1dict, listseqprobs, popgenos, meanevolvgnd #, indices #add foldfitdict to check fitnesses 
	return meanfitness,seqstargetlist,seqscommonlist,extreme0list,extreme1list,otherlist
if __name__ == "__main__":
	
	parameters = {}
	pair = int(sys.argv[1])
	fitlands = int(sys.argv[2])#1.random 2.hamming 3.random with inverse
	samplenum = int(sys.argv[3])
	plasticoption = int(sys.argv[4])
	gengap = int(sys.argv[5])
	generations = int(sys.argv[6])
	mu = float(sys.argv[7])
	minprob = int(sys.argv[8])
	minprob = minprob/100.
	maxprob = int(sys.argv[9])
	maxprob = maxprob/100.
	initoption = int(sys.argv[10]) #extremeoption -> 0: soft extremes, 1: true extremes, 2: random from seqscommon
	extremeoption = int(sys.argv[11]) #e.g. if we want to start simulation with target 0, the initial population should be at extreme of target 1
	Npop = int(sys.argv[12])
	run = int(sys.argv[13])
	
		
	path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"
	path1 = path + "fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+ "/"
	if not os.path.isdir(path): os.mkdir(path)
	if not os.path.isdir(path1): os.mkdir(path1)

	with open("/home/pg520/target_flipping/phenopairs.pkl","rb") as f: phenopairs = pickle.load(f)
	phenopair = phenopairs[pair]
	pheno0 = phenopair[0]
	pheno1 = phenopair[1]
	del phenopairs 
	

	with open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl","rb") as f: gpmap = pickle.load(f)

	seqscommon = np.load(path + "seqscommon.npy")
	#probs_common_0 = np.load(path + "probs_common_0.npy")
	#probs_common_1 = np.load(path + "probs_common_1.npy")

	#seqstarget = np.load(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy")
	seq_soft_extreme_0 = np.load(path + "seqssoftextreme0.npy")
	seq_soft_extreme_1 = np.load(path + "seqssoftextreme1.npy")
	seqstrueextreme0 = stringtoint(np.load(path + "seqstrueextreme0.npy"))
	seqstrueextreme1 = stringtoint(np.load(path + "seqstrueextreme1.npy"))
	seqscommon = stringtoint(np.load(path + "seqscommon.npy"))
	seqstarget = stringtoint(np.load(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy"))

	#seq_extreme_0_max = np.load(path + "seq_trueextreme_0_max.npy")
	#seq_extreme_1_max = np.load(path + "seq_trueextreme_1_max.npy")
	#probs_extreme_0_max = np.load(path + "/probs_trueextreme_0_max.npy")
	#probs_extreme_1_max = np.load(path + "/probs_trueextreme_1_max.npy")


	#with open("/rds/user/pg520/hpc-work/entropy.pkl","rb") as f: entropy = pickle.load(f)
	#with open("/rds/user/pg520/hpc-work/evgND.pkl","rb") as f: evgND = pickle.load(f)
	#with open("/rds/user/pg520/hpc-work/rhogND.pkl","rb") as f: rhogND = pickle.load(f)
	entropy = defaultdict(float)
	evgND = defaultdict(float)
	rhogND = defaultdict(float)
	
	#initial population
	#print("initoption:", initoption)
	#print("extremeoption:", extremeoption)

	initpop = []
	if initoption == 0:
		if extremeoption == 0: initpop = (1,seq_soft_extreme_0)
		if extremeoption == 1: initpop = (0,seq_soft_extreme_1)
	elif initoption == 1: 
		if extremeoption == 0: initpop = (1,seq_extreme_0_max)
	elif initoption == 1: 
		if extremeoption == 0: initpop = (1,seq_extreme_0_max)
		if extremeoption == 1: initpop = (0,seq_extreme_1_max)
	elif initoption == 2:
		p0val = 1
		p1val = 0
		while val == True:
			seq = random.choice(seqscommon)
			indexseq = list(seqscommon).index(seq)
			if seq not in seqstarget: val == False
		initpop = (random.choice([0,1]),seq)
	#print('targetoption and genotype start:', initpop)

	#meanrobust, meanentropy, meanfitness, adaptedphenonum, targetgenosgen, probstarg0dict, probstarg1dict, listseqprobs, popgenos, meanevolgnd = evodyn(rhogND=rhogND, evgND=evgND, entropy=entropy,targets=phenopair,initpop=initpop, genostargetset=seqstarget, seqscommon= seqscommon, samplenum = samplenum, plasticoption = plasticoption, gengap = gengap, L=12, gpmap = gpmap, mu = mu, generations = generations, Npop = Npop, fitlands = fitlands)
	
	meanfitness,seqstargetlist,seqscommonlist,extreme0list,extreme1list,otherlist = evodyn(rhogND=rhogND, evgND=evgND, entropy=entropy,targets=phenopair,initpop=initpop, genostargetset=seqstarget, seqscommon= seqscommon, extreme0 = seqstrueextreme0, extreme1 = seqstrueextreme1, samplenum = samplenum, plasticoption = plasticoption, gengap = gengap, L=12, gpmap = gpmap, mu = mu, generations = generations, Npop = Npop, fitlands = fitlands)
	
	del gpmap 
	dict = {'meanfitness': list(meanfitness.values()), 'seqscommon': seqscommonlist, 'seqstarget': seqstargetlist, 'extreme0': extreme0list, 'extreme1': extreme1list, 'other':otherlist} 
	df = pd.DataFrame(dict)
	del dict
	df.to_csv(path1+'seqdata_'+str(run)+'.csv')