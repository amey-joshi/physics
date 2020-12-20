#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 14:00:02 2020

@author: ajoshi
"""

import numpy as np
# import scipy.constants as sc
import matplotlib.pyplot as plt
# import pandas as pd

omega_1 = 1
omega_2 = 1
tau = .3
T_1 = 9
T_2 = 1
omega = 1
lambda_ = 0.1 # lambda is a Python keyword!!

if omega_1 == omega_2:
  omega_minus = np.pi/2
else:
  omega_minus = np.arctan((2 * omega * lambda_)/(omega_1 - omega_2))

def get_S(phi, theta, omega_minus):
  S2 = np.matrix([[np.cos(phi), np.sin(phi), 0, 0],
                  [-np.sin(phi), np.cos(phi), 0, 0],
                  [0, 0, np.cos(phi), np.sin(phi)],
                  [0, 0, -np.sin(phi), np.cos(phi)]])
# =============================================================================
#   S3 = np.matrix([[np.cos(theta), 0, np.sin(theta) * np.cos(omega_minus), np.sin(theta) * np.sin(omega_minus)],
#                   [0, np.cos(theta), np.sin(theta) * np.sin(omega_minus), -np.sin(theta) * np.cos(omega_minus)],
#                   [-np.sin(theta) * np.cos(omega_minus), -np.sin(theta) * np.sin(omega_minus), np.cos(theta), 0],
#                   [-np.sin(theta) * np.sin(omega_minus), np.sin(theta) * np.cos(omega_minus), 0, np.cos(theta)]])
# =============================================================================

  S3 = np.matrix([[np.cos(theta), np.sin(theta) * np.cos(omega_minus), 0, np.sin(theta) * np.sin(omega_minus)],
                  [-np.sin(theta) * np.cos(omega_minus), np.cos(theta), -np.sin(theta) * np.sin(omega_minus), 0],
                  [0, np.sin(theta) * np.sin(omega_minus), np.cos(theta), -np.sin(theta) * np.cos(omega_minus)],
                  [-np.sin(theta) * np.sin(omega_minus), 0, np.sin(theta) * np.cos(omega_minus), np.cos(theta)]])

  S = S2 @ S3

  return S

def get_beta(omega, T):
  # return sc.hbar * omega/(2 * sc.Boltzmann * T)
  return omega/(2 * T)

def get_v(beta):
  return 1/(2 * np.tanh(beta))

def build_V(v1, v2):
    V = np.matrix([[v1, 0, 0, 0],
                   [0, v1, 0, 0],
                   [0, 0, v2, 0],
                   [0, 0, 0, v2]])
    return V

def get_temperature(u, omega):
  x1 = np.arctanh(np.reciprocal(2 * u))
  # T = sc.hbar * omega/(2 * sc.Boltzmann) * np.reciprocal(x1)
  T = omega/2 * np.reciprocal(x1)

  return T

def run_simulation(max_iter, tau):
  beta_1 = get_beta(omega_1, T_1)
  beta_2 = get_beta(omega_2, T_2)
  v1 = get_v(beta_1)
  v2 = get_v(beta_2)
  u00 = np.zeros(max_iter)
  u22 = np.zeros(max_iter)

  V = build_V(v1, v2)
  u00[0] = V[0, 0]
  u22[0] = V[2, 2]

  phi = tau * (omega_1 + omega_2)/2
  theta = tau * np.sqrt(omega**2 * lambda_**2 + (omega_1 - omega_2)**2/4)

  for n in range(1, max_iter):
    S = get_S(phi, theta, omega_minus)
    U = S @ V @ np.transpose(S)
    u00[n] = U[0, 0]
    u22[n] = U[2, 2]
    # If both systems are considered to evolve. They do in Chimonidou-Sudarshan
    # scheme. They don't in the original Jayseetha Rau's scheme.
    V = build_V(u00[n], u22[n])

  return (u00, u22)

def plot_results(omega_1, omega_2, lambda_, T_1, T_2, u00, u22, max_iter, tau, rname):
  params = r'$\omega_1=$ {:0.1f}, $\omega_2=$ {:0.1f}, $\lambda=$ {:0.1f}, $\tau=$ {:0.1f}'. \
      format(omega_1, omega_2, lambda_, tau)
  T1 = get_temperature(u00, omega_1)
  T2 = get_temperature(u22, omega_2)
  t = np.zeros(max_iter)
  for n in range(1, max_iter):
    t[n] = tau * n

  plt.plot(t, T1, color = 'black', label = rf'$T_1={T_1}$')
  plt.plot(t, T2, color = 'red', label = rf'$T_2={T_2}$')
  plt.xlim(0, max(t))
  plt.legend()
  plt.ylabel(r'$kT/\hslash\omega$')
  plt.xlabel(r'$t$')
  plt.title(params)
  # plt.show()

  fname = f'fig_{rname}.png'
  plt.savefig(fname)
  plt.close()

def main():
  max_iter = 2500
  u00, u22 = run_simulation(max_iter, tau)
  rname = 'x1'
  plot_results(omega_1, omega_2, lambda_, T_1, T_2, u00, u22, max_iter, tau, rname)

if __name__ == '__main__':
  main()
