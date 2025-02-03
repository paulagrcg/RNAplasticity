import numpy as np 
import subprocess
import os
import math
from functions.basefunctions import *
from collections import defaultdict
from functionalRNA import *
import sys
import time
import pickle


def plasticitydistance(fRNAprob1, fRNAprob2):
    distances = {}
    for x, y,key1,key2 in zip(fRNAprob1.values(), fRNAprob2.values(), fRNAprob1.keys(), fRNAprob2.keys()):
        if key1 == key2:
            distance = np.sqrt((x - 0.5)**2 + (y - 0.5)**2)/np.sqrt(2)
            distances[key1] = distance
    return distances


def mutatesite(seq, site):
    mutations = {'A': ['C', 'T', 'G'], 'C': ['A', 'T', 'G'], 'G': ['A', 'T', 'C'], 'T': ['A', 'G', 'C']}
    mutationchoices = [seq[:site] + m + seq[site+1:] for m in mutations[str(seq[site])]]

    return mutationchoices

def scan_sites(seq, samplesize):
    seqs = []
    folds1 = []
    folds2 = []
    probs1 = []
    probs2 = []
    samplesizecount = 0

    seqoutput,seqfolds = suboptfolding(seq)
    seqs.append(seq)

    fold1 = seqoutput[0][0]
    fold2 = seqoutput[1][0]
    prob2 = float(seqoutput[1][2])
    prob1 = float(seqoutput[0][2])

    #folds1.append(fold1)
    #folds2.append(fold2)
    #probs1.append(prob1)
    probs2.append(prob2)
    
    site = 0
    seqlen = len(seq)
    while site < seqlen:
            site_neutral = False
            mutationchoices = mutatesite(seq, site)  # random mutation at site

            while not site_neutral and mutationchoices:
                m = np.random.randint(0, len(mutationchoices))
                seqmut = mutationchoices.pop(m)
                mutoutput, mfolds = suboptfolding(seqmut)
                if fold1 in mfolds and fold2 in mfolds:  # both are still in the plastic repertoire
                    #print(f"Match found: {samplesizecount}")
                    seqs.append(seqmut)
                    site_neutral = True
                    fold1 = mutoutput[0][0]
                    fold2 = mutoutput[1][0]
                    prob2 = float(mutoutput[1][2])
                    #prob1 = float(mutoutput[0][2])

                    #folds1.append(fold1)
                    #folds2.append(fold2)
                    #probs1.append(prob1)
                    probs2.append(prob2)

                    samplesizecount += 1
                    if samplesizecount == samplesize:
                        break
                    site += 1
                    #print(f"Match found: {samplesizecount}")
                    if site == seqlen:
                        site = 0
                    seq = seqmut
            if samplesizecount == samplesize:
                break
            if not mutationchoices:
                #print(f"No more mutation choices at site {site}")
                site += 1
                if site == seqlen:
                    site = 0

    #print(f"Total sequences found: {samplesizecount}")
    #return seqs, folds1, folds2, probs1, probs2
    return seqs, probs2


if __name__ == "__main__":

    
    #with open('../data/fRNAhammingdistance.pkl','rb') as f:
    #    fRNAhammingdistance = pickle.load(f)
    with open('../data/fRNAprob2.pkl','rb') as f:
        fRNAprob2 =  pickle.load(f)
    with open('../data/fRNAprob1.pkl','rb') as f:
        fRNAprob1 =  pickle.load(f)
    #with open('../data/fRNAfolds.pkl','rb') as f:
    #    fRNAfolds =  pickle.load(f)


    with open(f"../data/selected_for_analysis.pkl","rb") as f:
        selected_for_analysis = pickle.load(f)

    samplesize = int(sys.argv[1])
    seqposition = int(sys.argv[2])
    #i = 0
    #start = time.time()
    #seqs, folds1, folds2, probs1, probs2 = scan_sites(selected_sequences[i], samplesize)
    #end = time.time()
    
  
    p1 = fRNAprob1[tuple(selected_for_analysis[seqposition])]
    p2 = fRNAprob2[tuple(selected_for_analysis[seqposition])]

    #seqs, folds1, folds2, probs1, probs2 = scan_sites(selected_names[seqposition][1], samplesize)
    #start = time.time()
    seqs, probs2 = scan_sites(tuple(selected_for_analysis[seqposition])[1], samplesize)
    #end = time.time()
    #print(f"Time taken for site scanning: {end-start}")

    #print(f"Time taken for site scanning: {end-start}")
    with open(f"../data/selected_for_analysis_{seqposition}_ssize{samplesize}.pkl","wb") as f:
        pickle.dump({'seqs': seqs, 'probs2': probs2},f)
    
    