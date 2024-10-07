import numpy as np
from collections import defaultdict
import random
from quantitiesND4 import *
from gpmapfunctions import *
from functions.genotype_populations import *
from fitnesslandscape import fitness_hamming
import pandas as pd
from evodyn_plots001 import stringtoint
from collections import Counter

stn = {'A': '0', 'U': '1', 'C': '2', 'G': '3'}

def evodyn(rhogND,evgND,entropy,targets,initpop, seqscommon,extreme0, extreme1, samplenum = 1,plasticoption = False, gengap = 10, L=12, gpmap = {}, mu = 0.05, generations = 10000, Npop = 100, fitlands = 2):
        popgenos = {}
        #initpop is a tupple with initial target pheno and initial population list
        targetpheno = targets[initpop[0]]
        prob0 = [0]*(Npop)
        prob1 = [0]*(Npop)
        popgenos = np.empty((generations+1, Npop),dtype=int)

        probRWS = defaultdict(float)
        meanfitness = defaultdict(float)
        switchnum = 0
        listseqprobs = {}
        meanentropygenos = {}
        meanevolvgnd = {}
        meanrobustg = {}
        probstarg0dict = defaultdict(dict)
        probstarg1dict = defaultdict(dict)
        foldfitdict = defaultdict(float)
        seqscommonCounter = Counter(seqscommon)
        #seqstargetCounter = Counter(seqstarget)
        seqsex0Counter = Counter(extreme0)
        seqsex1Counter = Counter(extreme1)
        meanprobs0 = np.zeros(generations)
        meanprobs1 = np.zeros(generations)
        meanprobs0seqscom = np.zeros(generations)
        meanprobs1seqscom = np.zeros(generations)
        #seqstargetlist = np.zeros(generations)
        seqscommonlist = np.zeros(generations)
        extreme0list = np.zeros(generations)
        extreme1list = np.zeros(generations)
        otherlist = np.zeros(generations)
        alphalist = np.zeros(generations)
        for i in range(1,generations+1):
                listseqprobs[i] = defaultdict(dict)
                meanevolvgnd[i] = 0
                meanentropygenos[i] = 0
                meanrobustg[i] = 0
                
        for i in range(1,generations+1): # i is generation number
                
                genos = np.array([initpop[1]]*Npop) if i == 1 else newpop
                popgenos[i] = [int(''.join([stn[n] for n in geno])) for geno in genos] #taken from i-1 mutated population
                #print(popgenos)
                unique, counts = np.unique(popgenos[i], return_counts=True)
                freq = dict(zip(unique, counts))
                #print(freq)
                countseqscommon = sum(freq[seq] for seq in freq if seqscommonCounter[seq] == 1)
                #countseqstarget = sum(freq[seq] for seq in freq if seqstargetCounter[seq] == 1)
                countseqsex1 = sum(freq[seq] for seq in freq if seqsex1Counter[seq] == 1)
                countseqsex0 = sum(freq[seq] for seq in freq if seqsex0Counter[seq] == 1)
                countother = sum(freq[seq] for seq in freq if seqscommonCounter[seq] != 1 and seqsex0Counter[seq] != 1 and seqsex1Counter[seq] != 1)

                #seqstargetlist[i-1] = countseqstarget / Npop
                seqscommonlist[i-1] = countseqscommon / Npop
                extreme0list[i-1] = countseqsex0 / Npop
                extreme1list[i-1] = countseqsex1 / Npop
                otherlist[i-1] = countother / Npop
                alphalist[i-1] = len(freq.keys())/Npop
                folds = {}
                probs = {}
                #data of population phenotype ensembles and TPs probs
                for j in range(0,Npop):
                        phvsprobseq = extractnormalisedprobs(gpmap[genos[j]],L)
                        listseqprobs[i][j] = phvsprobseq
                        for f,p in phvsprobseq.items():
                            if f==targets[0][0:L]: 
                                prob0[j] = p
                            if f==targets[1][0:L]: 
                                prob1[j] = p
                        folds[j] = list(phvsprobseq.keys())
                        probs[j] = list(phvsprobseq.values())

                #data collection of target genotypes and quantities
                #for g in genos:
                #       if g in genostargetset:targetgenos +=1
                #       meanentropygenos[i] += entropy[g]/len(genos)
                #       meanevolvgnd[i] += evgND[g]/len(genos)
                #       meanrobustg[i] += rhogND[g]/len(genos)
                   #meanevolvgnd[i]+=evolvabilitygND(gpmap,g,4,12)[g]/len(genos)

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
                probstarg0seqscom = []
                probstarg1seqscom = []
                probstarg1 = {}
                count = 0
                #cadapt = 0
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
                        if seqscommonCounter[popgenos[i][count]] == 1: 
                            probstarg0seqscom.append(prob0[count])
                            probstarg1seqscom.append(prob1[count])
                        probstarg1[(ng,count)] = prob1[count]
                        probstarg0[(ng,count)] = prob0[count]
                        count+=1
        
        
                #mean fitness of population at gen i       
                meanfitness[i] = sum(list(popgen.values()))/len(popgen)
                #data collection of target probabilities
                if i == generations:
                    probstarg0end = probstarg0
                    probstarg1end = probstarg1
                #probstarg0dict[i] = probstarg0
                #probstarg1dict[i] = probstarg1
                meanprobs0[i-1] = np.mean(list(probstarg0.values()))
                meanprobs1[i-1] = np.mean(list(probstarg1.values()))
                meanprobs0seqscom[i-1] =np.mean(probstarg0seqscom) 
                meanprobs1seqscom[i-1] =np.mean(probstarg1seqscom)
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
                     
                #mutations for population i+1
                newpop = [mutation(seq, mu) for seq in newpop]

        #return meanfitness,seqstargetlist,seqscommonlist,extreme0list,extreme1list,otherlist,probstarg0dict,probstarg1dict,meanprobs0,meanprobs1
        return meanfitness,seqscommonlist,extreme0list,extreme1list,otherlist,probstarg0end,probstarg1end,meanprobs0seqscom,meanprobs1seqscom,meanprobs0,meanprobs1,alphalist
if __name__ == "__main__":
        pair = int(sys.argv[1])
        fitlands = int(sys.argv[2])#1.random 2.hamming 3.random with inverse
        samplenum = int(sys.argv[3])
        plasticoption = int(sys.argv[4])
        gengap = float(sys.argv[5])
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
        samples = 100 
        path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"    
        pathx = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/regimedata/"
        path1 = pathx + "fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+ "/"
        if not os.path.isdir(pathx): os.mkdir(pathx)
        if not os.path.isdir(path): os.mkdir(path)
        if not os.path.isdir(path1): os.mkdir(path1)

        with open("/home/pg520/target_flipping/phenopairs.pkl","rb") as f: phenopairs = pickle.load(f)
        phenopair = phenopairs[pair]
        pheno0 = phenopair[0]
        pheno1 = phenopair[1]
        del phenopairs 
        

        with open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl","rb") as f: gpmap = pickle.load(f)

        seqscommon = np.load(path + "seqscommon.npy")

        seq_soft_extreme_0 = np.load(path + "seqssoftextreme0.npy")
        seq_soft_extreme_1 = np.load(path + "seqssoftextreme1.npy")
        seqstrueextreme0 = stringtoint(np.load(path + "seqstrueextreme0.npy"))
        seqstrueextreme1 = stringtoint(np.load(path + "seqstrueextreme1.npy"))
        seqscommon = stringtoint(np.load(path + "seqscommon.npy"))


        entropy = defaultdict(float)
        evgND = defaultdict(float)
        rhogND = defaultdict(float)
        

        initpop = []
        if initoption == 0:
                if extremeoption == 0: initpop = (1,seq_soft_extreme_0)
                if extremeoption == 1: initpop = (0,seq_soft_extreme_1)

        dfs = []
        column_mean = defaultdict(list)
        column_std = defaultdict(list)
        probstot0 = []
        probstot1 = []
        for sample in range(0,samples):
                meanfitness,seqscommonlist,extreme0list,extreme1list,otherlist,probstarg0end,probstarg1end,meanprobs0seqscom,meanprobs1seqscom,meanprobs0,meanprobs1,alphalist = evodyn(rhogND=rhogND, evgND=evgND, entropy=entropy,targets=phenopair,initpop=initpop, seqscommon= seqscommon, extreme0 = seqstrueextreme0, extreme1 = seqstrueextreme1, samplenum = samplenum, plasticoption = plasticoption, gengap = gengap, L=12, gpmap = gpmap, mu = mu, generations = generations, Npop = Npop, fitlands = fitlands)
                dictsampleprobs ={'meanfitness': list(meanfitness.values()),'seqscommon': seqscommonlist,'extreme0': extreme0list, 'extreme1': extreme1list, 'other':otherlist,'meanprobs0seqscom':meanprobs0seqscom, 'meanprobs1seqscom':meanprobs1seqscom, 'meanprobs0':meanprobs0, 'meanprobs1':meanprobs1,'alphalist': alphalist}
                df = pd.DataFrame(dictsampleprobs)
                dfs.append(df)
                for v0,v1 in zip(list(probstarg0end.values()), list(probstarg1end.values())):
                    probstot0.append(v0)
                    probstot1.append(v1)
                
    
        concatenated_df = pd.concat(dfs)
        dfmean = concatenated_df.groupby(concatenated_df.index).mean()
        dfstd = concatenated_df.groupby(concatenated_df.index).std()
        np.save(path1 + "probstarg0.npy", np.array(probstot0))
        np.save(path1 + "probstarg1.npy", np.array(probstot1))
        dfmean.to_csv(path1 + 'meanprobsalphadata_'+"fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+'.csv')
        dfstd.to_csv(path1 + 'stdprobsalphadata_'+"fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+'.csv')
