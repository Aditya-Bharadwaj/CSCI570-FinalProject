# CSCI570 - Final Project

Implementation of Sequence Alignment Algorithms

## Description

---

This project presents two implementations of 2-sequence alignment - a common dynamic programming problem.

The first implementation is a straightforward dynamic programming approach that computes and stores the values of every subproblem in-memory.

The second implementation is a more memory efficient version that integrates both Divide-and-Conquer and Dynamic Programming paradigms without sacrificing on time complexity.

This project was part of CSCI 570 - Analysis of Algorithms course at the University of Southern California.

## Authors

---

Aditya Bharadwaj

Hiranmaya Gundu

Karan Agrawal

## Running the code

To run the shell scripts run the following commands:

Efficient version:

```sh
./5889672313_3586943483_9346073553_efficient.sh
```

Basic version:

```sh
./5889672313_3586943483_9346073553_basic.sh
```

Both the shell scripts also take optional input parameter:

```sh
./5889672313_3586943483_9346073553_efficient.sh new_input.txt
```

To generate the plots, run:

```sh
python3 generate_plots.py
```
