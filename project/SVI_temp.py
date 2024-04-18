import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

def OptSVI1test(i, Vol, w, p, K, S, TT):
    """
    Introduction
        Given the implied volatilities in market, weights, exercise prices,
        initial asset price, Expiry time,
        this function calculate objective function for SVI calibration.

    Input
        Vol (a list)                  : implied volatilities in market
        w (a list)                    : weights
        p (a list)                    : SVI parameters
        K (a list)                    : exercise prices
        S (a float)                   : initial asset price
        TT (a list)                   : Expiry time

    Output
        e (a float)                   : objective function
    """
    
    lK = len(K)
    ImpVol = np.zeros(lK)  # implied volatility
    k = np.zeros(lK)  # logarithmic forward moneyness
    e = 0
    for j in range(lK):
        ImpVol[j] = Vol[i][j]
        k[j] = np.log(K[j] / S)
        SVIVol = np.sqrt(max(p[0] + p[1] * (p[2] * (k[j] - p[3]) + np.sqrt(max((k[j] - p[3]) ** 2 + p[4], 0))), 0)) / np.sqrt(TT[i] / 365)  # max is to avoid error
        e += w[i][j] * (SVIVol / ImpVol[j] - 1) ** 2
    return e


