import subprocess
import os
import math
from functions.basefunctions import *
from collections import defaultdict
import numpy as np
import time 
import pickle 

kbT_RNA = 0.6163207755


def suboptfolding(seq, energy_range=9.2448116325):
    # suboptimals and their energies
    rnasubopt_path = os.path.expanduser("~/ViennaRNA/bin/RNAsubopt")
    
    if not os.path.isfile(rnasubopt_path):
        raise FileNotFoundError(f"RNAfold not found at {rnasubopt_path}")
    
    rnasubopt = subprocess.Popen(
        [rnasubopt_path, "-e", str(energy_range)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={"PATH": os.environ["PATH"]}
    )
    
    stdout, stderr = rnasubopt.communicate(input=bytes(seq, 'UTF-8'))
    
    if rnasubopt.returncode != 0:
        raise RuntimeError(f"RNAfold failed with error: {stderr.decode()}")
    
    output = stdout.decode().strip().split('\n')    
    
    if not os.path.isfile(rnasubopt_path):
        raise FileNotFoundError(f"RNAfold not found at {rnasubopt_path}")
    
    #ensemble free energy

    rnasubopt_path = os.path.expanduser("~/ViennaRNA/bin/RNAfold")

    rnasubopt = subprocess.Popen(
        [rnasubopt_path, "-p"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={"PATH": os.environ["PATH"]}
    )
    
    stdout, stderr = rnasubopt.communicate(input=bytes(seq, 'UTF-8'))
    
    if rnasubopt.returncode != 0:
        raise RuntimeError(f"RNAfold failed with error: {stderr.decode()}")
    
    outputF = stdout.decode().strip().split('\n')
    Fens = float(outputF[2].replace('[','').replace(']','').strip().split()[1])

    # Extract the fold, energy and Boltzmann probabilities
    subopts = []
    folds = []
    probs = []
    energies = []
    suboptfolds = []
    if len(output) >= 2:
        for i in range(2, len(output)):
            foldandenergy = output[i].split()
            fold = (foldandenergy[0])
            folds.append(fold)
            energy = float(foldandenergy[1])
            energies.append(energy)
            prob = math.exp(Fens - energy) / kbT_RNA
            probs.append(prob)
            
    else:
        raise RuntimeError("Unexpected RNAfold output format")
    try:
        probs = normalize(probs)
    except ZeroDivisionError:
        return None
    
    for f,e,p in zip(folds,energies,probs):
        subopts.append((f,e,p))
        suboptfolds.append(f)

    subopts = sorted(subopts, key=lambda x: x[1])
    return subopts,suboptfolds



def read_fasta(file_path):
    sequences = {}
    with open(file_path) as f:
        seqname = ''
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seqname:
                    sequences[seqname] = seq
                seqname = line.lstrip('>')
                seq = ''
            else:
                seq += line
        if seqname:
            sequences[seqname] = seq
    return sequences  

def filtersequences(sequences,L=50):
    filtered_sequences = {}
    for seqname, seq in sequences.items():
        if len(seq) <= L:
            filtered_sequences[seqname] = seq
    return filtered_sequences


if __name__ == "__main__":
    sequence = "GGGAAAUCC"
    energy_range = 9.2448116325
    output = suboptfolding(sequence, energy_range)
    print(output)
    
    sequences = read_fasta('../data/sequence.fasta')
    sequences = filtersequences(sequences)

    hammingdistance = defaultdict(float)
    prob2= defaultdict(float)
    folds = defaultdict(list)
    prob1 = defaultdict(float)
    for seqname, seq in sequences.items():
        print(seqname)
        output, outfolds= suboptfolding(seq, energy_range)
        if output is None:
            continue
        folds[(seqname,seq)] = [output[0][0], output[1][0]]
        hammingdistance[(seqname,seq)] = hamming(output[0][0], output[1][0])
        prob2[(seqname,seq)] = float(output[1][2])
        prob1[(seqname,seq)] = float(output[0][2])
    
    with open("/home/pg520/phenodistance/data/fRNAhammingdistance.pkl","wb") as f:
        pickle.dump(hammingdistance,f)
    with open("/home/pg520/phenodistance/data/fRNAprob2.pkl","wb") as f:
        pickle.dump(prob2,f)
    with open("/home/pg520/phenodistance/data/fRNAprob1.pkl","wb") as f:
        pickle.dump(prob1,f)
    with open("/home/pg520/phenodistance/data/fRNAfolds.pkl","wb") as f:
        pickle.dump(folds,f)

    