import pickle
import os
import time
from multiprocessing import Pool, cpu_count
from sitescanning_manual import scan_sites

# Load shared data only once
with open('fRNAprob1.pkl', 'rb') as f:
    fRNAprob1 = pickle.load(f)
with open('fRNAprob2.pkl', 'rb') as f:
    fRNAprob2 = pickle.load(f)
with open('not_in_seqanalysed.pkl', 'rb') as f:
    names = pickle.load(f)

samplesize = 5000

def process_seq(seqposition):
    name_tuple = tuple(names[seqposition])
    out_file = f"to_analyse_{name_tuple}_ssize{samplesize}.pkl"

    if os.path.exists(out_file):
        return f"⏩ Skipped {seqposition}"

    try:
        start = time.time()
        seqs, probs1, probs2 = scan_sites(name_tuple[1], samplesize)

        with open(out_file, "wb") as f:
            pickle.dump({'seqs': seqs, 'probs1': probs1, 'probs2': probs2}, f)

        end = time.time()
        return f"✅ Done {seqposition} ({name_tuple}) in {end - start:.1f}s"

    except Exception as e:
        return f"❌ Failed {seqposition} ({name_tuple}): {e}"

if __name__ == "__main__":
    total = len(names)
    num_processes = max(1, cpu_count() - 1)
    print(f"Launching {num_processes} parallel workers...")

    # You can restrict to a subset if desired, e.g. range(0, 100)
    with Pool(num_processes) as pool:
        results = pool.map(process_seq, range(total))

    # Log output
    with open("parallel_log.txt", "w") as f:
        for r in results:
            print(r)
            f.write(r + "\n")
