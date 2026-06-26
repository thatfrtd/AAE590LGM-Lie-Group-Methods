import numpy as np
import matplotlib.pyplot as plt

def plot_convergence(solution):
    '''
    Plot convergence of objective, delta, and defect
    '''

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    a = np.array([solution["info"][i]["J"] for i in range(solution["converged_i"])])

    ax1.plot(np.arange(0, solution["converged_i"]), [solution["info"][i]["J"] for i in range(solution["converged_i"])])
    ax1.set_xlabel("Iteration")
    #ax1.set_ylabel("Objective")
    ax1.set_title("Objective vs Iteration")
    ax1.grid(True)

    ax2.plot(np.arange(0, solution["converged_i"]), solution["delta"][0:solution["converged_i"]])
    ax2.set_xlabel("Iteration")
    #ax2.set_ylabel("Delta")
    ax2.set_yscale("log")
    ax2.set_title("Stopping Criteria vs Iteration")
    ax2.grid(True)

    ax3.plot(np.arange(0, solution["converged_i"] + 1), solution["defect"][0:(solution["converged_i"] + 1)])
    ax3.set_xlabel("Iteration")
    #ax3.set_ylabel("Defect")
    ax3.set_yscale("log")
    ax3.set_title("Defect vs Iteration")
    ax3.grid(True)

    plt.show()