import numpy as np 
import subprocess
import os
import math
from basefunctions import *
from collections import defaultdict
from functionalRNA import *
import sys
import pickle
import signal


kbT_RNA = 0.6163207755


def suboptfolding(seq, energy_range=9.2448116325):
    # suboptimals and their energies
    rnasubopt_path = os.path.expanduser("/usr/local/bin/RNAsubopt")
    
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

    rnasubopt_path = os.path.expanduser("/usr/local/bin/RNAfold")

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




# def save_checkpoint(seqposition, samplesize, site, seqs, probs1, probs2):
#     """Save the current state of seqs, probs1, and probs2 to a checkpoint file."""
#     checkpointlen = samplesize-len(seqs)
#     checkpoint_file = f"../data/checkpoint_{checkpointlen}_{site}_{seqposition}_ssize{samplesize}.pkl"
#     with open(checkpoint_file, "wb") as f:
#         pickle.dump({'seqs': seqs, 'probs1': probs1, 'probs2': probs2}, f)
#     print(f"Checkpoint saved to {checkpoint_file}")

# def handler(signum, frame):
#     """Handle termination signals by saving the current state and exiting."""
#     print("Received termination signal. Saving state...")
#     save_checkpoint()
#     sys.exit(0)

def plasticitydistance(fRNAprob1, fRNAprob2):
    distances = {}
    for x, y,key1,key2 in zip(fRNAprob1.values(), fRNAprob2.values(), fRNAprob1.keys(), fRNAprob2.keys()):
        if key1 == key2:
            distance = np.sqrt((x - 0.5)**2 + (y - 0.5)**2)/np.sqrt(2)
            distances[key1] = distance
    return distances

def hamming_plasticity_optimal_old(fRNAhammingdistance, fRNAfolds, distances):
    scores = {}
    for name in fRNAhammingdistance.keys(): 
        scores[name] = fRNAhammingdistance[name] - distances[name] 
    # Remove entries from scores where the sequence of dots is in fRNAfolds
    filtered_scores = {key: value for key, value in scores.items() if '.' * len(key[1]) not in fRNAfolds[key]}

    # Select the data points with the highest scores
    selected_data = sorted(filtered_scores.items(), key=lambda item: item[1], reverse=True)

    # Extract the Hamming values and their corresponding minimum distances
    selected_hamming_values = [fRNAhammingdistance[name] for name, score in selected_data]
    selected_min_distances = [distances[name] for name, score in selected_data]
    selected_names = [name for name, score in selected_data]

    return selected_hamming_values, selected_min_distances, selected_names

def hamming_plasticity_optimal_old(fRNAhammingdistance, fRNAfolds, distances):
    scores = {}
    for name in fRNAhammingdistance.keys(): 
        scores[name] = distances[name] 
    # Remove entries from scores where the sequence of dots is in fRNAfolds
    filtered_scores = {key: value for key, value in scores.items() if '.' * len(key[1]) not in fRNAfolds[key] and fRNAhammingdistance[key] >= 0.5}

    # Select the data points with the highest scores
    selected_data = sorted(filtered_scores.items(), key=lambda item: item[1])

    # Extract the Hamming values and their corresponding minimum distances
    selected_hamming_values = [fRNAhammingdistance[name] for name, score in selected_data]
    selected_min_distances = [distances[name] for name, score in selected_data]
    selected_names = [name for name, score in selected_data]

    return selected_hamming_values, selected_min_distances, selected_names

def hamming_plasticity_optimal(fRNAhammingdistance, fRNAfolds, distances):
    selected_names = []
    selected_hamming_values =[]
    selected_min_distances = []
    for name in fRNAhammingdistance.keys():
        if '.' * len(name[1]) not in fRNAfolds[name]:
            selected_names.append(name)
            selected_hamming_values.append(fRNAhammingdistance[name])
            selected_min_distances.append(distances[name])
    selected_names = np.array(selected_names)
    selected_hamming_values = np.array(selected_hamming_values)
    selected_min_distances = np.array(selected_min_distances)


    # Define the coordinates of the two points
    #point0 = (0.30688980382826875, 0.7777777777777778)
    #point1 = (0.20333026388667144, 0.16)
    point2 = (0.4702129756344944,0.7741935483870968)

    point1= (0.15382552893408447, 0.18181818181818182)
    point0 = (0.24487580716620733, 0.8095238095238095)
    
    m = (point1[1] - point0[1]) / (point1[0] - point0[0])
    b = point0[1] - m * point0[0]

    m1= (point2[1] - point0[1]) / (point2[0] - point0[0])
    b1 = point0[1] - m1 * point0[0]

    point0a = (0.30688980382826875,0.7777777777777778) #('FR481122||Fly small RNA', 'TCGAATCCGAAGATTGCA') 0.30688980382826875 0.7777777777777778
    point1a =  (0.20333026388667144,0.16)

    ma = (point1a[1] - point0a[1]) / (point1a[0] - point0a[0])
    ba = point1a[1] - ma * point1a[0]

    # Calculate the slope (m) and y-intercept (b) of the line
   

    # Define a function to determine if a point is to the left of the line
    def is_left_of_line(x, y, m, b):
        return y > m * x + b
    def is_hamming_greater_than_small(h,h0=0.15):
        return h > h0

    left_mask = is_left_of_line(selected_min_distances, selected_hamming_values, m, b) & is_hamming_greater_than_small(selected_hamming_values)
    left_mask_a = is_left_of_line(selected_min_distances, selected_hamming_values, ma, ba) & is_hamming_greater_than_small(selected_hamming_values)
    # Filter points that are to the left of the second line
    top_mask = is_left_of_line(selected_min_distances, selected_hamming_values, m1, b1)

    left_names = selected_names[left_mask]
    print(len(left_names))
    left_names_a = selected_names[left_mask_a]
    print(len(left_names_a))
    left_names_minus_a = []
    for i in left_names_a:
        if i in left_names: continue
        else: left_names_minus_a.append(i)

    print(len(left_names_minus_a))
    #left_min_distances = selected_min_distances[left_mask]
    #left_hamming_values = selected_hamming_values[left_mask]
    #combined_mask = left_mask & top_mask
    #topminuscombined_mask = top_mask & ~combined_mask
    #topminuscombined_min_distances = selected_min_distances[topminuscombined_mask]
    #topminuscombined_hamming_values = selected_hamming_values[topminuscombined_mask]
    #topminuscombined_names = selected_names[topminuscombined_mask]
    #joined_names = np.concatenate((topminuscombined_names, left_names))
    #joined_hamming_values = np.concatenate((topminuscombined_hamming_values, left_hamming_values))
    #joined_min_distances = np.concatenate((topminuscombined_min_distances, left_min_distances))

    #return left_names, topminuscombined_names,joined_names
    return left_names, left_names_a,left_names_minus_a

def mutatesite(seq, site):
    try:
        mutations = {'A': ['C', 'T', 'G'], 'C': ['A', 'T', 'G'], 'G': ['A', 'T', 'C'], 'T': ['A', 'G', 'C']}
        mutationchoices = [seq[:site] + m + seq[site+1:] for m in mutations[str(seq[site])]]
    except KeyError:
        print(f"KeyError: {seq[site]} not in mutations")
        print(seq)
        return []
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
    probs1.append(prob1)
    probs2.append(prob2)
    
    site = 0
    seqlen = len(seq)
    while site < seqlen:
            site_neutral = False
            mutationchoices = mutatesite(seq, site)  # random mutation at site
            if not mutationchoices:
                #print(f"No more mutation choices at site {site}")
                site += 1
                if site == seqlen:
                    break
                continue
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
                    prob1 = float(mutoutput[0][2])

                    #folds1.append(fold1)
                    #folds2.append(fold2)
                    probs1.append(prob1)
                    probs2.append(prob2)

                    samplesizecount += 1
                    #if samplesizecount % 100 == 0:  # Save progress every 10 samples
                        #save_checkpoint(seqposition, samplesize, site, seqs, probs1, probs2)
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
    return seqs, probs1, probs2


if __name__ == "__main__":

    import time
    start = time.time()
    #with open('../data/fRNAhammingdistance.pkl','rb') as f:
    #    fRNAhammingdistance = pickle.load(f)
    with open('../data/fRNAprob2.pkl','rb') as f:
        fRNAprob2 =  pickle.load(f)
    with open('../data/fRNAprob1.pkl','rb') as f:
        fRNAprob1 =  pickle.load(f)
    #with open('../data/fRNAfolds.pkl','rb') as f:
    #    fRNAfolds =  pickle.load(f)

    #distances = plasticitydistance(fRNAprob1, fRNAprob2)
    #start = time.time()
    #left_names, topminuscombined_names,joined_names = hamming_plasticity_optimal(fRNAhammingdistance, fRNAfolds, distances)
    #left_names, left_names_a,left_names_minus_a = hamming_plasticity_optimal(fRNAhammingdistance, fRNAfolds, distances)
    
    #selected_hamming_values, selected_min_distances, selected_names = hamming_plasticity_optimal_old(fRNAhammingdistance, fRNAfolds, distances)
    #site-scanning
    #end  = time.time()
    #print(f"Time taken for selecting: {end-start}")
    #with open(f"../data/selected_names_left.pkl","wb") as f:
    #    pickle.dump(left_names,f)
    #with open(f"../data/selected_names_left_a.pkl","wb") as f:
    #    pickle.dump(left_names_a,f)
    with open(f"not_in_seqanalysed.pkl","rb") as f:
        names = pickle.load(f)

    samplesize = 5000
    seqposition = 0
    #missing_numbers = np.load("../main/missing_numbers.npy")
    #seqposition = missing_numbers[int(sys.argv[2])]
    #print(seqposition)
    #i = 0
    #start = time.time()
    #seqs, folds1, folds2, probs1, probs2 = scan_sites(selected_sequences[i], samplesize)
    #end = time.time()
    
    #p1 and p2 for p-value calculation
    #p1 = fRNAprob1[selected_names[seqposition]]
    #p2 = fRNAprob2[selected_names[seqposition]]
    p1 = fRNAprob1[tuple(names[seqposition])]
    p2 = fRNAprob2[tuple(names[seqposition])]

    seqs, probs1, probs2 = scan_sites(tuple(names[seqposition])[1], samplesize)
    # try:
    #     seqs, probs1, probs2 = scan_sites(tuple(names[seqposition])[1], samplesize)
    # except Exception as e:
    #     print(f"An error occurred: {e}")
    #     #save_checkpoint()
    #     sys.exit(1)

    #end = time.time()
    #print(f"Time taken for site scanning: {end-start}")
    #print(f"Time taken for site scanning: {end-start}")
    with open(f"to_analyse_{tuple(names[seqposition])}_ssize{samplesize}.pkl","wb") as f:
        pickle.dump({'seqs': seqs, 'probs1': probs1, 'probs2': probs2},f)
    
    end = time.time()
    print(f"Time taken for site scanning: {end-start}")

