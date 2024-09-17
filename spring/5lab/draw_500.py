#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def draw(filenames, labels):

    plt.rcParams["legend.markerscale"] = 2.0
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['PT Sans']
    plt.rcParams['font.size'] = '12'
    plt.rcParams["legend.loc"] = "upper left"
    plt.rcParams["legend.fontsize"] = "8"
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(10*cm, 7*cm))
    ax = fig.add_subplot(111)
    ax.set_title("500000 элементов")
    ax.set(xlabel="Количество $p$ процессов", ylabel="Ускорение")
    ax.label_outer()
    ax.grid()
    ax.xaxis.set_tick_params(direction='in', which='both')
    ax.yaxis.set_tick_params(direction='in', which='both')

    for (fname, datalabel) in zip(filenames, labels):
        data = np.loadtxt(fname)
        x = data[:, 0]
        y = data[:, 1]

        ax.plot(x, y, '-o', markersize=1, linewidth=0.5, label=datalabel)

    ax.plot(range(2, 9), range(2, 9), '-', c="blue", linewidth=0.5, label="Линейное ускорение")

    plt.tight_layout()
    ax.legend()
    fig.savefig('chart_500.png', dpi=300)

if __name__ == "__main__":
    filenames = ["500000and200.dat", "500000and400.dat", "500000and600.dat", "500000and800.dat", "500000and1000.dat", "500000and1200.dat"]
    labels = ["threshold = 200", "threshold = 400", "threshold = 600", "threshold = 800", "threshold = 1000", "threshold = 1200"]
    draw(filenames, labels)