import numpy as np
import pickle
import operator
from collections import defaultdict
from itertools import combinations, product
import pandas as pd
import sys
def hamming_distance(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Strings must have equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))/12.

def hammingdistance_sets(dataset,seed):
    distance = []
    np.random.seed(seed)
    for tup in dataset:
        str1, str2 = tup
        distance.append(hamming_distance(str1, str2))

    distance_sort = np.unique(np.array(np.sort(distance)))
    print(distance_sort)

    max_distance = max(distance) 
    print('max:',max_distance)
 
    median_distance = np.median(distance)  
    print('med:',median_distance)

    min1_distance = distance_sort[1]
    print('min1', min1_distance)

    min_distance = min(distance)
    print('min:',min_distance)

    max_distance_inds = []  
    min_distance_inds = [] 
    min1_distance_inds = [] 
    median_distance_inds = []  
    for i, d in enumerate(distance):
        if d == max_distance:
            max_distance_inds.append(i)
        if d == min_distance:
            min_distance_inds.append(i)
        if d == median_distance:
            median_distance_inds.append(i)
        if d == min1_distance: 
            min1_distance_inds.append(i)
    max_distance_index = np.random.choice(max_distance_inds)
    min_distance_index = np.random.choice(min_distance_inds)
    min1_distance_index = np.random.choice(min1_distance_inds)
    median_distance_index = np.random.choice(median_distance_inds)
    return [dataset[max_distance_index], dataset[median_distance_index], dataset[min_distance_index], dataset[min1_distance_index]], [max_distance, median_distance, min_distance, min1_distance]

def phenosizescale(neutralsets):
    NSS_sorted = sorted(neutralsets.items(), key=lambda x: np.log10(x[1]), reverse=True)
    del NSS_sorted[0] #remove the deleterious set
    maxsetscale = np.log10(NSS_sorted[0][1])
    dictNSSscale = defaultdict(list)
    for i in range(3, int(maxsetscale)+1):
        for pheno, size in NSS_sorted:
            if np.log10(size) < i+1 and np.log10(size) >= i:
                dictNSSscale[i].append(pheno)
    
    complexcomplex = list(product(dictNSSscale[3], dictNSSscale[3]))
    complexsimple = list(product(dictNSSscale[3], dictNSSscale[5]))
    simplesimple = list(product(dictNSSscale[5], dictNSSscale[5]))
    return complexcomplex, complexsimple, simplesimple

if __name__ == "__main__":
    seed = int(sys.argv[1])
    with open("NDsetsize.pkl","rb") as f:
            neutralsets= pickle.load(f)
    
    complexcomplex, complexsimple, simplesimple = phenosizescale(neutralsets)
    cc, cc_dist = hammingdistance_sets(complexcomplex, seed)
    cs, cs_dist  = hammingdistance_sets(complexsimple, seed)
    ss, ss_dist = hammingdistance_sets(simplesimple, seed)
    cc = cc[:-2] + cc[-1:]
    ss = ss[:-2] + ss[-1:]
    cs = cs[:3]
    cc_dist = cc_dist[:-2] + cc_dist[-1:]
    ss_dist = ss_dist[:-2] + ss_dist[-1:]
    cs_dist = cs_dist[:3]
    np.save("cc_cs_ss_pairs_seed_"+str(seed)+".npy", np.array([cc,cs,ss]))
    np.save("cc_cs_ss_distances_seed_"+str(seed)+".npy", np.array([cc_dist,cs_dist,ss_dist]))

    print(np.load("cc_cs_ss_pairs_seed_"+str(seed)+".npy"))
    print(np.load("cc_cs_ss_distances_seed_"+str(seed)+".npy"))
