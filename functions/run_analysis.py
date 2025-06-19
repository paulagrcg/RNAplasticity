import pickle
import os
import time
from sitescanning_manual import scan_sites

# Load data
with open('fRNAprob1.pkl', 'rb') as f:
    fRNAprob1 = pickle.load(f)
with open('fRNAprob2.pkl', 'rb') as f:
    fRNAprob2 = pickle.load(f)
with open('not_in_seqanalysed.pkl', 'rb') as f:
    names = pickle.load(f)

samplesize = 5000
total = len(names)

for seqposition in range(total):
    name_tuple = tuple(names[seqposition])
    out_file = f"to_analyse_{name_tuple}_ssize{samplesize}.pkl"

    if os.path.exists(out_file):
        print(f"Skipping {seqposition}: already exists.")
        continue

    try:
        print(f"Processing seqposition {seqposition}/{total - 1}: {name_tuple}")
        start = time.time()

        # Run scan
        seqs, probs1, probs2 = scan_sites(name_tuple[1], samplesize)

        with open(out_file, "wb") as f:
            pickle.dump({'seqs': seqs, 'probs1': probs1, 'probs2': probs2}, f)

        end = time.time()
        print(f"Done {seqposition} in {end - start:.2f} seconds.")

    except Exception as e:
        print(f"Failed at {seqposition} ({name_tuple}): {e}")

