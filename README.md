# RNA plasticity 

-Data: organised by Figure

Data points per generation:

- seqscommon : plastic genotypes abundance 
- extreme0 : genotypes with only target 0
- extreme1: genotypes with only target 1 
- other: rest of gentoypes, noise. 
- meanfitness: mean fitness of the population
- meanprobs0: population mean of Boltzmann probability on target 0
- meanprobs1: population mean of Boltzmann probability on target 1

(ignore seqstarget, meanprobs0seqscom, meanprobs1seqscom, and alphalist)

- NDsetsize.pkl is the dictionary of secondary structures' NSS for the RNA12 ND GP map.
- pairs.csv is the mu (X) and delta t (Y) combinations for the parameter_sweep .py files to create the heatmaps
- phenopairs.pkl is an initial dictionary with targetpairs from which we select the pair 130 (key is 130) for the example in Fig.2.

-Functions: 

- evodyn001_* are .py files for Fig.2 example targetpair. 
- evodyn_regimes_sweep_* are .py files dealing with the targetpairs from different categories. (cc, cs, ss) (0,1,2)
- cc is complex-complex (both small NSS), cs is complex-simple (mixed NSS) and simple-simple (large NSS).
- 0 is small, 1 is medium and 2 is large Hamming distance.


gmapfunctions, phenotypesearchfinal, and targetpairs and data_genotype_targetpairs deal with categorising targetpairs according to NSS and Hamming distance as well as producing their genotype populations

evodyn_seqs.ipynb is a notebook used to produce plots and test data.

-Plots: directory with .png of plots
