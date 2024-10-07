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
        for run in range(0,9):

            path = "/rds/user/pg520/hpc-work/targetflipping/seed_"+str(seed)+"/"
            pathx = path + "regime_ratio_samplenum_run_"+str(run)+"/"
            if not os.path.isdir(pathx): os.mkdir(pathx)
            complexitylist = np.load(path + "cc_cs_ss_pairs_seed_"+str(seed)+".npy")
            sweep = ["cc_0","cc_1","cc_2","cs_0","cs_1","cs_2","ss_0","ss_1","ss_2"]
            flattened_complexitylist = [pair for sublist in complexitylist for pair in sublist]
            phenopair = flattened_complexitylist[run]
            pairphenopath = sweep[run]
            
            dir_name = path + pairphenopath+"/"
            column_names = ['meanfitness', 'seqscommon', 'seqstarget','extreme0', 'extreme1', 'other']
            df = pd.read_csv('/home/pg520/pairs.csv',float_precision='round_trip')

            mu_values = df['X']
            gengap_values = df['Y']
            #path = "/rds/user/pg520/hpc-work/targetflipping/pair_130/regimedata/"
            #        pathx = path + "regime_ratio/"
            #        if not os.path.isdir(pathx): os.mkdir(pathx)

            regime_alpha = {}
            regime_fitness = {}
            regime_dist= {}
            regime_ratio = {}

            for mu,gengap in zip(mu_values,gengap_values):
                                            filename = dir_name + "meanprobsalphadata_fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+ ".csv"
                                            regime_dict = {}
                                            mu = float(mu)
                                            gengap = float(gengap)
                                            df = pd.read_csv(filename)

                                            meanfitness = np.mean(df['meanfitness'][generations//2:generations])
                                            meanalpha = np.mean(df['alphalist'][generations//2:generations])

                                            list_meanprobs0 = list(df['meanprobs0'][generations//2:generations])
                                            list_meanprobs1 = list(df['meanprobs1'][generations//2:generations])
                                            list_dist = []
                                            for p0,p1 in zip(list_meanprobs0,list_meanprobs1):
                                                last_point = (p0, p1)
                                                target_point = (0.5, 0.5)
                                                distance_last_to_target = np.sqrt((last_point[0] - target_point[0])**2 + (last_point[1] - target_point[1])**2)/(np.sqrt(2)/2)
                                                list_dist.append(np.round(distance_last_to_target,2))

                                            plasticitydist = np.mean(list_dist)
                                            regime_alpha[mu,gengap] = meanalpha
                                            regime_fitness[mu,gengap] = meanfitness
                                            regime_dist[mu,gengap] = plasticitydist
                                            sum_seqscommon = np.sum(df['seqscommon'][5000:])
                                            sum_other =  np.sum(df['extreme0'][5000:]) + np.sum(df['extreme1'][5000:]) + np.sum(df['other'][5000:])
                                            regime_ratio[mu,gengap] = float(sum_seqscommon/sum_other)

            with open(pathx + "regime_ratio_samplenum_"+str(samplenum)+".pkl","wb") as f: pickle.dump(regime_ratio,f)

