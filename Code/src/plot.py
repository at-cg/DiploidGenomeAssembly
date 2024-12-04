import click
import math
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

H = 61317360
G = 2 * H
maxl = 25000
EPS = 0.01

# # 1000
# len_double = 16810
# mat_triple = 9317; mat_triple_F = 1
# pat_triple = 2784; pat_triple_F = 1
# statA = 16810; statA_F = 1
# statB = 25831; statB_F = 1
# statC = 18130; statC_F = 4

# # 20000
len_double = 16810
mat_triple = 9317; mat_triple_F = 1
pat_triple = 7371; pat_triple_F = 1
statA = 16810; statA_F = 1
statB = 277321; statB_F = 1
statC = 285826; statC_F = 4

def interleaved(path, L):
    maxv1 = -1; maxv1_F = 1
    maxv2 = -1; maxv2_F = 1
    with open(path, 'r') as file:
        for line in file:
            tokens = line.strip().split()
            l1 = int(tokens[0]); l2 = int(tokens[1])
            if(l1 < L - 1 and l2 < L - 1):# Type 2
                val = (l1 + l2) / 2
                if(val > maxv2):
                    maxv2 = val; maxv2_F = 1
                elif(val == maxv2):
                    maxv2_F += 1
            elif(l1 < L - 1 and l2 >= L - 1):# Type 1
                val = l1
                if(val > maxv1):
                    maxv1 = val; maxv1_F = 1
                elif(val == maxv1):
                    maxv1_F += 1
    return  max(G / (2 * (L - maxv1 - 1)) * math.log(maxv1_F / (2 * EPS)), 
                G / (4 * (L - maxv2 - 1)) * math.log(maxv2_F / (2 * EPS)))

def triple(x, xf, L):
    return G / (3 * (L - x - 1)) * math.log(xf / (2 * EPS))

def verify(path, N, L):
    val = 0
    with open(path, 'r') as file:
        for line in file:
            g = int(line.strip().split()[0])
            for j in range(max(0, g - L + 1), min(g, L - 1) + 1):
                val += 2 * (1 / G * (g - j + 1) / G * (2 * math.exp(-(N - 2)/ G * (L - (g - j) - 1))))
            
            if(L < g + 2):
                val += 2 * math.exp(-N / G * (L - (g - L + 1) - 1))
    
    return val <= EPS

@click.command()
@click.option('--gap_path', '-gpath', help='path to the file storing the length of gap b/w two adjacent loci',
              type=click.Path(exists=True))
@click.option('--mat_interleaved_path', '-mpath', help='path to the maternal interleaved file',
              type=click.Path(exists=True))
@click.option('--pat_interleaved_path', '-ppath', help='path to the paternal interleaved file',
              type=click.Path(exists=True))
@click.option('--special_interleaved_path', '-spath', help='path to the special interleaved file',
              type=click.Path(exists=True))
@click.option('--epsilon', '-ep', help='substitution error rate', default=0.01,
              type=float)
@click.option('--output_path', '-opath', help='path to the file that will save the plot',
              type=click.Path(exists=False))
def main(gap_path, mat_interleaved_path, pat_interleaved_path, special_interleaved_path, epsilon, output_path):
    global EPS
    GPATH = gap_path
    MPATH = mat_interleaved_path
    PPATH = pat_interleaved_path
    SPATH = special_interleaved_path
    EPS = epsilon
    OUTPATH = output_path

    minit = 9319
    xi = []
    yi = []

    ming = 16812
    xg = []
    yg = []

    mind = 16811
    xd = []
    yd = []

    for i in range(minit, maxl):
        L = i

        # Information-theoretic
        xi.append(L)
        val_LW = G / L * math.log(G / (L * EPS))
        val = val_LW # Lander-Waterman

        val = max(val, interleaved(MPATH, L)) # mat-interleaved
        val = max(val, interleaved(PPATH, L)) # pat-interleaved
        val = max(val, interleaved(SPATH, L)) # special-interleaved

        val = max(val, triple(mat_triple, mat_triple_F, L)) # mat-triple
        val = max(val, triple(pat_triple, pat_triple_F, L)) # pat-triple

        yi.append(val / val_LW)

        # Greedy
        if(i >= ming):
            xg.append(L)
            valg = val # Information-theoretic 

            valg = max(valg, G / (2 * (L - statA - 1)) * math.log(statA_F / EPS)) # Well-bridging
            valg = max(valg, G / (5 * (L - statB / 5 - 1)) * math.log(statB_F / EPS))
            valg = max(valg, G / (8 * (L - statC / 4 - 1)) * math.log(statC_F / EPS))

            # assert verify(GPATH, valg, L) # Checking for pairwise overlap

            yg.append(valg / val_LW)

        # DBG
        if(i >= mind + 1):
            xd.append(L)
            vald = max(val / val_LW, 1 / (1 - (len_double + 1) / L))
            yd.append(vald)
        
        
    plt.plot(xi, yi, label='Information-Theoretic Lower Bound')
    plt.plot(xg, yg, label='Greedy Algorithm')
    plt.plot(xd, yd, label='de Bruijn Graph based Algorithm')
    plt.vlines(x = len_double, ymin = 0, ymax = 20, color = 'k', linestyle = 'dashed', linewidth = 1)
    plt.text(len_double - 500, 10, 'Maximum Length of a Double Repeat', color = 'k', rotation = 90, verticalalignment = 'center')
    plt.axhline(y = 1, color = 'k', linestyle = 'dashed', linewidth = 1)
    plt.xlabel('Read Length')
    plt.ylabel('Minimum Normalized Coverage')
    plt.ylim(0, 20)
    plt.title('Feasibility Plot for Diploid Human T2T Chromosome 19')
    plt.legend()
    plt.savefig(OUTPATH, dpi=1000)

if __name__ == '__main__':
    main()

