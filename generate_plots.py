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
    problem_comparison = dict()
    x_string = ''
    y_string = ''
    for i in range(0, input_size):
        # print(i)
        problem_comparison[i] = dict()
        x_string += random.choice(LETTERS)
        y_string += random.choice(LETTERS)
        start_time = time.time()
        tracemalloc.start()
        fn(x_string, y_string)
        problem_comparison[i]['cpu_time'] = time.time() - start_time
        problem_comparison[i]['memory_util'] = tracemalloc.get_traced_memory(
        )
        tracemalloc.stop()
    return problem_comparison


def build_plots(input_size=100):
    dnc_problem_sizes = generate_randstring(
        code.divide_and_conquer_alignment, input_size)
    brute_problem_sizes = generate_randstring(
        code.sequence_alignment_brute, input_size)
    xaxis = [i * 2 for i in range(input_size)]
    yaxis_cpu_brute = np.array([brute_problem_sizes[key]
                                ['cpu_time'] for key in range(input_size)])
    yaxis_cpu_dnc = np.array([dnc_problem_sizes[key]
                              ['cpu_time'] for key in range(input_size)])
    yaxis_mem_brute = np.array([brute_problem_sizes[key]
                                ['memory_util'][1]/1000 for key in range(input_size)])
    yaxis_mem_dnc = np.array([dnc_problem_sizes[key]
                              ['memory_util'][1]/1000 for key in range(input_size)])
    # Generate Problem Size vs CPU
    plt.xlabel("Problem Size")
    plt.ylabel("CPU Time (in seconds)")
    plt.title("Comparison of Sequence Alignment Implementations")
    plt.plot(xaxis, yaxis_cpu_brute, label="Brute Force")
    plt.plot(xaxis, yaxis_cpu_dnc, label="D&C")
    plt.legend()
    # plt.savefig('cpu_problemsize.svg', format='svg')
    plt.savefig('cpu_problemsize.png', format='png')

    # Clear figure for new plot
    plt.clf()

    # Generate Problem Size vs Memory
    plt.xlabel("Problem Size")
    plt.ylabel("Memory Utilization (in kb)")
    plt.title("Comparison of Sequence Alignment Implementations")
    plt.plot(xaxis, yaxis_mem_brute, label="Brute Force")
    plt.plot(xaxis, yaxis_mem_dnc, label="D&C")
    plt.legend()
    # plt.savefig('mem_problemsize.svg', format='svg')
    plt.savefig('mem_problemsize.png', format='png')


if __name__ == '__main__':
    input_size = 200
    build_plots(input_size)
