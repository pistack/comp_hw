import argparse
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot result')
    parser.add_argument('mode', 
                        choices=['plot', 'compare'],
                        help='mode: 1. plot: plot one result, ' + \
                        '2. compare: compare with reference file')
    parser.add_argument('filename', 
                        help='filename to plot')
    parser.add_argument('--ref',
                        help='filename of reference result required by ' + \
                              'compare mode')

    arg = parser.parse_args()
    A = np.genfromtxt(arg.filename)

    if arg.mode == 'plot':
        plt.plot(A[:, 0], A[:, 1], label='result')
        plt.legend()
        plt.show()
    else:
        ref = np.genfromtxt(arg.ref)
        plt.plot(A[:, 0], A[:, 1], marker='o', markersize=5,
                 linestyle='none', label='result')
        plt.plot(ref[:, 0], ref[:, 1], label='ref')
        plt.legend()
        plt.show()
