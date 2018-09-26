# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:57:31 2018

@author: cana5
"""
# Import the following programs:
import numpy as np
import matplotlib.pyplot as plt

# Define the following to convert values to seconds:
seconds_in_year = 86400 * 365
seconds_in_day = 86400

# This is the independent variable (time):
MJD = np.linspace(50000, 55000, 1000)
MJD_seconds = MJD * seconds_in_day

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

theta = 0.049
chi = 1.5517
tau_p = 485.56 * seconds_in_day
psi_initial = 3.8709

# Define 'nudot' in order to make it easier to plot the function:
nudot = term(MJD_seconds, tau_age, P, n, t_ref, theta, chi, tau_p, psi_initial)

# Plot the function:
fig = plt.figure(figsize=(10,6))
plt.plot(MJD_seconds, nudot, label='PRECESSION FUNCTION')
plt.legend()
plt.title('PRECESSION FUNCTION FOR PULSAR B1828-11')
plt.xlabel('TIME (s)')
plt.ylabel('$\dot\\nu$ (Hz/s)')
txt = ("Graph showing the precession function, where the pulsar's spin-down "
       "goes to three distinct states, and then repeats.")
plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', 
            fontsize=10)
plt.tight_layout()
plt.savefig('precession_constant_period_function.pdf')