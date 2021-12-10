import tracemalloc
import sys
import time
from code import divide_and_conquer_alignment, process_input

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 <filename.py> <input.txt>")
        sys.exit()
    start_time = time.time()
    tracemalloc.start()
    X_orig, Y_orig = process_input(sys.argv[1])
    X_aligned, Y_aligned, opt_cost = divide_and_conquer_alignment(
        X_orig, Y_orig)
    print(X_aligned[: 50], X_aligned[-50:])
    print(Y_aligned[: 50], Y_aligned[-50:])
    print(opt_cost)
    print(time.time() - start_time)
    print(tracemalloc.get_traced_memory()[1]/1000)
    with open("output.txt", "w") as f:
        f.write(X_aligned[: 50] + " " + X_aligned[-50:] + "\n")
        f.write(Y_aligned[: 50] + " " + Y_aligned[-50:] + "\n")
        f.write(str(opt_cost) + "\n")
        f.write(str(time.time() - start_time) + "\n")
        f.write(str(tracemalloc.get_traced_memory()[1]/1000) + "\n")
    tracemalloc.stop()
