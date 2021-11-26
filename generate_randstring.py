import random
import matplotlib.pyplot as plt
import numpy as np
import code
import time
import tracemalloc
def generate_randstring(input_size = 10000):
    LETTERS = ['A', 'G', 'T', 'C']
    x_string = ''
    y_string = ''
    for i in range(0, input_size):
        x_index = random.randint(0, 3)
        y_index = random.randint(0, 3)
        x_string += LETTERS[x_index]
        y_string += LETTERS[y_index]
    return (x_string, y_string)

def generate_datapoints(fn, input_size=10000):
    cpu_times = list()
    X, Y = generate_randstring(input_size)
    # mem_usage = list()
    for i in range(input_size):
        for j in range(input_size):
            start_time = time.time()
            tracemalloc.start()
            fn(X[:i], Y[:j])
            cpu_times.append(time.time() - start_time)
    return cpu_times    

if __name__ == '__main__':
    input_size = 10000
    cpu_datapoints = generate_datapoints(code.sequence_alignment, 100)
    print(cpu_datapoints)
    

    