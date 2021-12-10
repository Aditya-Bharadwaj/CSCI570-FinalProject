import tracemalloc
import sys
import time
from code import sequence_alignment_brute, process_input

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 <filename.py> <input.txt>")
        sys.exit()
    start_time = time.time()
    tracemalloc.start()
    X_orig, Y_orig = process_input(sys.argv[1])
    X_aligned, Y_aligned, opt_cost = sequence_alignment_brute(
        X_orig, Y_orig)
    with open("output.txt", "w") as f:
        f.write(X_aligned[: 50] + " " + X_aligned[-50:] + "\n")
        f.write(Y_aligned[: 50] + " " + Y_aligned[-50:] + "\n")
        f.write(str(opt_cost) + "\n")
        f.write(str(time.time() - start_time) + "\n")
        f.write(str(tracemalloc.get_traced_memory()[1]/1000) + "\n")
    tracemalloc.stop()
