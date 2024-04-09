import itertools as it
from collections import Counter
import pickle
from collections import defaultdict
import time
import sys
import os
import math

def normalize(probs):
    prob_factor = 1 / sum(probs)
    return [prob_factor * p for p in probs]

def extractnormalisedprobs(pboltzlist,L):
    probsnorm = []
    for p in pboltzlist:
        probsnorm.append(float(p[L+2:]))
    prob = normalize(probsnorm)
    return {pboltzlist[pi][0:L]: prob[pi] for pi in range(0,len(pboltzlist))}

def mutationalneighbours(seq):
    mutations = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    return [seq[:j] + m + seq[j+1:] for j in range(0,len(seq)) for m in mutations[str(seq[j])]]

#%%%%%%%%%%%%%%%%%%%%%%%%%% ND quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



def robustnessgND(gpmap,seq,K,L):

         rhogndict = defaultdict(float)
         #first normalize probabilites:
         phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
         neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
         for phenotype, probg in phvsprobseq.items():
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                try:
                    probpgmut = phvsprobmut[phenotype] #if th is key does not exist it will jump into the KeyError bit, so we can have a single lookup to check if the Key exists and what value it maps to if it exists
                    rhogndict[seq] += probpgmut*probg/((K-1)*L) 
                except KeyError:
                    continue
         del neighbourvsphvsprob

         return rhogndict
def evolvabilitygND(gpmap,seq,K,L):
         evgndict = defaultdict(float)
         phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
         neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
         for phenotype, probg in phvsprobseq.items():
            probfold = defaultdict(lambda:1)
            evgndictp = 0
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut,probgmut in phvsprobmut.items():
                    if phenotype != phenomut:
                        probfold[phenomut] *=(1-probgmut)
                    else: continue
            for prob in probfold.values():
                evgndictp+=(1-prob)
            evgndict[seq]+=evgndictp*probg
            
         del neighbourvsphvsprob

         return evgndict
     
     
def NDfolds_setsizes(gpmap,K,L):
    
    folddict = defaultdict(list)
    folds = []
    NDsetsize = defaultdict(float)
    
    for seq in gpmap.keys():
        phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
        for phenotype,probg in phvsprobseq.items():
            folds.append(phenotype)
            NDsetsize[phenotype] += probg
            folddict[phenotype].append([seq,probg])
            
    return NDsetsize,folddict,folds



def robustnesspND(folddict,gpmap):
    rhopndict = defaultdict(float)
    
    for f,seqprobs in foldsdict.items():
        probgtot = 0
        for seq_p in seqprobs:
            seq=seq_p[0]
            probg =float(seq_p[1])
            probgtot += float(seq_p[1])
            neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
            
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut,probpgmut in phvsprobmut.items():
                    if f == phenomut:
                        rhopndict[f] += probpgmut*probg
            del neighbourvsphvsprob

        rhopndict[f]/=probgtot
        
    return rhopndict
    


def evolvabilitypND(folddict,gpmap): 

    evolp = defaultdict(float)

    for f,seqprobs in foldsdict.items():
        probgprime = defaultdict(lambda:1)
        for seq_p in seqprobs:
            seq=seq_p[0]
            probg =float(seq_p[1])
            neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
            
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut,probpgmut in phvsprobmut.items():
                    if f != phenomut:
                        probgprime[phenomut] *=(1-probpgmut*probg)
                    else: continue
            del neighbourvsphvsprob

    for prob in probgprime.values():
        evolp[f]+=(1-prob)
           
    return evolp

#%%%%%%%%%%%%%%%%%%%%%%%%%% averages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def sampleGP(dictsuboptRNA12):
    resolutiongp = {}
    for seq,subopt in dictsuboptRNA12.items():
        foldList = []
        probsList = []
        phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
        for phenotype,probg in phvsprobseq.values():
            foldList.append(phenotype)
            probsList.append(probg)
        resultgp[seq] = random.choices(foldList, weights = probsList, k=1)[0]
    return resolutiongp
"""
a_file = open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl", "rb")
dictsuboptRNA12 = pickle.load(a_file)
L=12
i = sys.argv[1]
evg = defaultdict(float)
f = open('/rds/user/pg520/hpc-work/seqsfiles/file'+i+'.txt',"rb")
lines = f.readlines()
for line in lines:
        seq = line[0:L]
        seq = seq.decode('utf-8')
        evg[seq] = evolvabilitygND(dictsuboptRNA12,seq,4,12)[seq]
a_file = open("evgND"+i+".pkl","wb")
pickle.dump(evg,a_file)
"""
