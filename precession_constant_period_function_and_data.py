# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 14:09:58 2018

@author: cana5
"""
# Import the following programs:
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define the following to convert values to seconds:
seconds_in_year = 86400 * 365
seconds_in_day = 86400

# This is the independent variable (time):
#MJD = np.linspace(50000, 55000, 1000)
#MJD_seconds = MJD * seconds_in_day

# DATA SECTION START
file_name = 'data/1828-11_100vf_fig.dat'
cols = (0, 3, 4)
names = ["MJD", "F1", "F1_err"]
df = pd.read_csv(file_name, sep=' ', usecols=cols, header=None, names=names,
                 dtype=None, skipinitialspace=True, index_col=False)

MJD_data = df.MJD.values
nudot_data = df.F1.values * 1e-15
MJD_seconds = MJD_data * seconds_in_day
# DATA SECTION END

# Define the terms in the function:
def term(MJD, tau_age, P, n, t_ref, theta, chi, tau_p, psi_initial):
    
    # Define the function in parts:
    psi = (2 * np.pi * (MJD - t_ref) / tau_p) + psi_initial
    mean = 1 / (tau_age * P)
    
    a = -1
    b = (n * (MJD - t_ref)) / tau_age
    c = 2 * theta * (np.cos(chi) / np.sin(chi)) * np.sin(psi)
    d = (- theta**2 / 2) * np.cos(2 * psi)
    
    # Define the function as a whole:
    return mean * (a + b + c + d)

# Assign values for the appropriate terms:
tau_age = 213827.91 * seconds_in_year
P = 0.405
n = 16.08
t_ref = 49621 * seconds_in_day

theta = 0.05155138281293769
chi = 1.5541505500913704
tau_p = 41906142.83814837
psi_initial = 3.8694873697852117

# Define 'nudot' in order to make it easier to plot the function:
nudot = term(MJD_seconds, tau_age, P, n, t_ref, theta, chi, tau_p, psi_initial)

# Plot the function: 
fig = plt.figure(figsize=(10,6))
plt.plot(MJD_seconds, nudot, label='DATA FITTED PRECESSION FUNCTION')
plt.scatter(MJD_seconds, nudot_data, color='gray', marker='.', 
            label='PULSAR DATA')
plt.title('PRECESSION FUNCTION FOR PULSAR B1828-11')
plt.legend()
plt.xlabel('TIME (s)')
plt.ylabel('$\dot\\nu$ (Hz/s)')
txt = ("Graph showing the data fitted precession function, and the pulsar's "
       "data points.")
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', 
            fontsize=10)
plt.tight_layout()
plt.savefig('precession_constant_period_function_and_data.pdf')