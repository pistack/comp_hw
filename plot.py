import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    A = np.genfromtxt(sys.argv[1])
    plt.plot(A[:, 0], A[:, 1], marker='o', markersize=5,
             linestyle='none')
    plt.show()
