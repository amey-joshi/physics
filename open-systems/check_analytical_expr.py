#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:57:17 2020

@author: ajoshi
"""

import numpy as np
import matplotlib.pyplot as plt

max_iter = 2500
omega_1 = 1
omega_2 = 1
tau = .3
T_1 = 9
T_2 = 1
omega = 1
lambda_ = 0.1 # lambda is a Python keyword!!

phi = tau * (omega_1 + omega_2)/2
theta = tau * np.sqrt(omega**2 * lambda_**2 + (omega_1 - omega_2)**2/4)

if omega_1 == omega_2:
  omega_minus = np.pi/2
else:
  omega_minus = np.arctan((2 * omega * lambda_)/(omega_1 - omega_2))

b1 = np.zeros(max_iter)
b2 = np.zeros(max_iter)

def get_beta_inv(v, omega):
  return 0.5 * omega / np.arctanh(np.reciprocal(2 * v))

v1 = 0.5/np.tanh(omega_1/(2 * T_1))
v2 = 0.5/np.tanh(omega_2/(2 * T_2))
alpha = 1
b1[0] = T_1
b2[0] = T_2

for n in range(1, max_iter):
  curr_v1 = v1 + (v2 - v1) * np.sin(theta)**2 * np.sin(omega_minus)**2
  curr_v2 = v2 - (v2 - v1) * np.sin(theta)**2 * np.sin(omega_minus)**2

  b1[n] = get_beta_inv(curr_v1, omega_1)
  b2[n] = get_beta_inv(curr_v2, omega_2)

  v1 = curr_v1
  v2 = curr_v2

params = r'$\omega_1=$ {:0.1f}, $\omega_2=$ {:0.1f}, $\lambda=$ {:0.1f}, $\tau=$ {:0.1f}'. \
      format(omega_1, omega_2, lambda_, tau)
t = np.zeros(max_iter)
for n in range(1, max_iter):
  t[n] = tau * n
plt.plot(t, b1, color = 'black', label = 'Oscillator 1')
plt.plot(t, b2, color = 'red', label = 'Oscillator 2')
plt.ylabel(r'$kT/\hslash\omega$')
plt.xlabel(r'$t$')
plt.legend()
plt.title(params)
plt.show()
