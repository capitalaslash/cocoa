#! /usr/bin/env python3.11

import os
import sys

import numpy as np
import matplotlib.pyplot as plt


def main(filename):

    dirname = os.path.dirname(filename)

    # start = 0.0
    # end = 1.0
    # nElems = 10
    time = 0.0
    finalTime = 1.0
    dt = 0.1
    outFile = "fd1d"
    with open(filename, "r") as f:
        for line in f:
            tokens = line.split(" ")
            if tokens[0] == "start_time:":
                time = float(tokens[1])
            elif tokens[0] == "final_time:":
                finalTime = float(tokens[1])
            elif tokens[0] == "dt:":
                dt = float(tokens[1])
            elif tokens[0] == "out_file:":
                outFile = tokens[1]
            # elif tokens[0] == "start:":
            #     start = float(tokens[1])
            # elif tokens[0] == "end:":
            #     end = float(tokens[1])
            # elif tokens[0] == "n_elems:":
            #     nElems = int(tokens[1])

    it = 0
    iters = [0]
    while time < finalTime - 1.0e-6:
        it += 1
        iters.append(it)
        time += dt

    print(f"iter: {iters}\n")

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))

    for it in iters:
        data = np.loadtxt(f"{dirname}/{outFile}.{it}.dat")
        ax.plot(data[:, 0], data[:, 1], label=f"{it}")

    ax.legend()
    ax.grid()

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} configFile")
        exit(1)
    print(f"configFile: {sys.argv[1]}")

    main(sys.argv[1])
