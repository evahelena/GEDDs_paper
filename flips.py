# Flips

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 23:43:50 2015

@author: lepraunch
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import math
import sys

fileIN = sys.argv[2]

with open(fileIN) as f:
    data = f.read().splitlines()
    
size = len(data)
lowess = sm.nonparametric.lowess

chrnum = str(sys.argv[1])

def randomizeArr(arr):
    size = len(arr)
    newArr = []
    for i in range(size):
        newArr.append(arr[i])
    for i in range(size):
        ran = int(np.random.uniform(0, size))
        newArr.append(newArr[ran])
        newArr.pop(ran)
        size -= 1
    return newArr

def fillCurve(arr):
    size = len(arr)
    fillaPos = []
    fillaNeg = []
    for i in range(size):
        if(arr[i] > 0):
            fillaPos.append(arr[i])
            fillaNeg.append(0)
        else:
            fillaPos.append(0)
            fillaNeg.append(arr[i])
    return [fillaPos, fillaNeg]

def numOfRegions(arr):
    first = arr[0]
    count = 1
    for i in range(1, len(arr)):
        if(math.copysign(1, first) != math.copysign(1, arr[i])):
            count += 1
            first = -1*first
    return count

def energy(arr):
    size = len(arr)
    E = arr[0]*arr[1]
    E += arr[size-1]*arr[size-2]
    for i in range(1,size-1):
        E += arr[i]*arr[i-1]+arr[i]*arr[i+1]
    return E

def plotFDR(bins):
    state5 = False
    state1 = False

    prob1 = float(bins[0][0])
    prob5 = float(bins[0][0])

    for i in range(1, len(bins[0])):

        if (prob1/boots <= 0.001 and state1 == False):
            prob1+=float(bins[0][i])
        elif(prob1/boots > 0.001 and state1 == False):
            state1 = True
            plt.plot([bins[1][i-1], bins[1][i-1]],[0,float(max(bins[0]))/7], 'r',linewidth=2, label='3sigma 0.1%')

        if (prob5/boots <= 0.023 and state5 == False):
            prob5+=float(bins[0][i])
        elif(prob5/boots > 0.023 and state5 == False):
            state5 = True
            plt.plot([bins[1][i-1], bins[1][i-1]],[0,float(max(bins[0]))/3], 'y',linewidth=2, label='2sigma 2.3%')

        if (prob1/boots <= 0.999 and state1 == True):
            prob1+=float(bins[0][i])
        elif(prob1/boots > 0.999 and state1 == True):
            plt.plot([bins[1][i], bins[1][i]],[0,float(max(bins[0]))/7], 'r',linewidth=2)
            break

        if (prob5/boots <= 0.977 and state5 == True):
            prob5+=float(bins[0][i])
        elif(prob5/boots > 0.977 and state5 == True):
            plt.plot([bins[1][i], bins[1][i]],[0,float(max(bins[0]))/3], 'y',linewidth=2)
            state5 = 3

for k in range(1):
    count = 0
    x = []
    y = []
    startLibrary = {}
    for i in range(size):
        values = data[i].split("\t")
        if (values[19] == chrnum and values[3] != '#VALUE!' and values[3] != 'NA'):
            startLibrary[int(values[20])] = values

    Starts = startLibrary.keys()
    sortedStarts = Starts
    sortedStarts.sort()

    firstStart = int(startLibrary[sortedStarts[0]][20])

    for i in range(len(sortedStarts)):
        start = sortedStarts[i]
        end = int(startLibrary[sortedStarts[i]][21])
        length = end - start
        logFold = float(startLibrary[sortedStarts[i]][3])
        fold = float(startLibrary[sortedStarts[i]][4])

        if(logFold < 0):
            pass
        else:
            pass
        x.append(i)
        y.append(logFold)
        firstStart = firstStart + length

    f = 15
    z = lowess(y, x, frac=1./f, is_sorted= True)
    a = []
    b = []
    for i in range(len(z)):
        a.append(z[i][0])
        b.append(z[i][1])
    initb = b
    #
    # DO BOOTSTRAPING
    #

    Nregions = []
    boots = 100000
    print("number of regions in ordered chromosome "+chrnum+":  " + str(numOfRegions(initb)))
    for i in range(boots):
        if(i%(boots/100) == 0):
            print(str(100*i/boots)+'%')
        nA = randomizeArr(y)
        Nregions.append(numOfRegions(nA))

    fig = plt.figure()
    bins = plt.hist(Nregions, 40)
    plotFDR(bins)
    plt.axvline(numOfRegions(y), linewidth=2, color='g', label = 'original value')
    plt.legend()


    nameOffile = fileIN+'_chr'+chrnum+'_frac'+str(f)+'-'+str(numOfRegions(y))+'_'+str(boots)+'_boots_flips.png'
    fig.savefig(nameOffile, dpi=300)



