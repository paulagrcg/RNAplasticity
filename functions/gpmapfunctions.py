import numpy as np
import random

def normalize(probs):
    prob_factor = 1 / sum(probs)
    return [prob_factor * p for p in probs]

def mutationalneighbours(seq):
    mutations = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    return [seq[:j] + m + seq[j+1:] for j in range(0,len(seq)) for m in mutations[str(seq[j])]]

def extractnormalisedprobs(pboltzlist,L):
    probsnorm = []
    for p in pboltzlist:
        probsnorm.append(float(p[L+2:]))
    prob = normalize(probsnorm)
    return {pboltzlist[pi][0:L]: prob[pi] for pi in range(0,len(pboltzlist))}

def mutation(seq,mu):
    basemutation = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    for n in range(0,len(seq)):
        basemut = random.choices([1,0], weights=[mu,1-mu],k=1)[0]
        if basemut==1: #mutate a base
            mut = random.choices(basemutation[seq[n]],k=1)[0]
            p = seq[:n] + str(mut) + seq[n+1:]
            seq = p
    return seq
