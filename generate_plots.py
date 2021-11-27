import random
import matplotlib.pyplot as plt
import numpy as np
import code
import time
import tracemalloc


def generate_randstring(fn, input_size=1000):
    """
    Generates two random equal-length strings and performs sequence alignment 
    """
    LETTERS = ['A', 'G', 'T', 'C']
    x_string = ''
    y_string = ''
    cpu_times = list()
    problem_size = list()
    for i in range(0, input_size):
        x_string += random.choice(LETTERS)
        y_string += random.choice(LETTERS)
        start_time = time.time()
        tracemalloc.start()
        fn(x_string, y_string)
        cpu_times.append(time.time() - start_time)
        problem_size.append(i)
        # print(i)
    return (cpu_times, problem_size)


if __name__ == '__main__':
    input_size = 16
    cpu_datapoints, problem_sizes = generate_randstring(
        code.sequence_alignment_brute, input_size)
    xaxis = np.array(problem_sizes)
    yaxis = np.array(cpu_datapoints)

    plt.xlabel("Problem Size")
    plt.ylabel("CPU Time (in seconds)")
    plt.title("Comparison of Sequence Alignment Implementations")
    plt.plot(xaxis, yaxis, label="Brute Force")
    plt.legend()
    plt.savefig('cpu_problemsize.svg', format='svg')
    plt.savefig('cpu_problemsize.png', format='png')
