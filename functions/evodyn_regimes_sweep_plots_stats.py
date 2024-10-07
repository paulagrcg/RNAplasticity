import numpy as np 
from scipy.interpolate import interp1d
from scipy.integrate import simps
import numpy as np
import pandas as pd
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import sys
import pickle
import re
import os 

if __name__ == "__main__":

        seed = int(sys.argv[1])
        #complexity = int(sys.argv[2])
        #distance = int(sys.argv[3])
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
        samples = int(sys.argv[13])
        run = int(sys.argv[14])#out of 9 (0-8)
#       run = int((run/800.)*8)
        for seed in [1, 3, 5, 6, 7, 8, 9, 10, 11, 12]:
            for run in np.arange(9):

                path = "/rds/user/pg520/hpc-work/targetflipping/seed_"+str(seed)+"/"

                complexitylist = np.load(path + "cc_cs_ss_pairs_seed_"+str(seed)+".npy")
                sweep = ["cc_0","cc_1","cc_2","cs_0","cs_1","cs_2","ss_0","ss_1","ss_2"]
                flattened_complexitylist = [pair for sublist in complexitylist for pair in sublist]
                phenopair = flattened_complexitylist[run]
                pairphenopath = sweep[run]
                
                dir_name = path + pairphenopath
                
                pattern = r'.*gengap(.*?)_mu(.*?)_.*'
                regime_fitness = {}
                regime_dist= {}
                for filename in os.listdir(dir_name):
                        if filename.startswith('meanprobsalphadata_'):
                                match = re.search(pattern, filename)
                                if match:
                                        regime_dict = {}
                                        gengap, mu = match.groups()
                                        mu = float(mu)
                                        gengap = float(gengap)
                                        df = pd.read_csv(dir_name + '/' + filename)

                                        meanfitness = np.mean(df['meanfitness'][generations//2:generations])

                                        list_meanprobs0 = list(data['meanprobs0'])
                                        list_meanprobs1 = list(data['meanprobs1'])
                                        first_point = (list_meanprobs0[0], list_meanprobs1[0])
                                        last_point = (list_meanprobs0[-1], list_meanprobs1[-1])
                                        target_point = (0.5, 0.5)
                                        distance_last_to_target = np.sqrt((last_point[0] - target_point[0])**2 + (last_point[1] - target_point[1])**2)/(np.sqrt(2)/2)
                                        plasticitydist = np.round(distance_last_to_target,2)

                                        regime_fitness[mu,gengap] = meanfitness
                                        regime_dist[mu,gengap] = plasticitydist


                pathx = path + "regime_fitness_seed_"+str(seed)+"/"
                if not os.path.isdir(pathx): os.mkdir(pathx)
                with open(pathx + "regime_fitness_"+str(run)+".pkl","wb") as f: pickle.dump(regime_fitness,f)
                

                pathx = path + "regime_dist_seed_"+str(seed)+"/"
                if not os.path.isdir(pathx): os.mkdir(pathx)
                with open(pathx + "regime_dist_"+str(run)+".pkl","wb") as f: pickle.dump(regime_dist,f)

