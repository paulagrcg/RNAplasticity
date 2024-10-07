import numpy as np
import sys
import pandas as pd
from collections import defaultdict
import glob
import os 

if __name__ == "__main__":

        parameters = {}
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

        pathx = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/regimedata/"
        path1 = pathx + "fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+ "/"
        column_names = ['meanfitness', 'seqscommon', 'seqstarget','extreme0', 'extreme1', 'other']
        column_mean = defaultdict(list)
        for column_name in column_names:

                csv_files = glob.glob(path1 + 'seqdata_*.csv')

                dfs = []

                for file in csv_files:
                        df = pd.read_csv(file)
                        dfs.append(df[column_name])

                concatenated_df = pd.concat(dfs, axis=1)
                for i in range(0,generations):
                        column_mean[column_name].append(concatenated_df.iloc[i].mean())

        dfmean = pd.DataFrame(column_mean)
        pathxx = pathx + 'means/'
        if not os.path.isdir(pathxx): os.mkdir(pathxx)
        dfmean.to_csv(pathxx+'meandata_'+str(run)+'.csv')
