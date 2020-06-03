#!/usr/bin/env python
# coding: utf-8

# Ramsey Fringes

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import rb_constants as rbc

g31 = rbc.g31
g32 = rbc.g32
gf = rbc.gf
gs = rbc.gs

first_step=0.5
max_step=10

# eoms

def get_w1(I1):
    return rbc.w1(I1)


def get_w2(I2):
    return rbc.w2(I2)


def dP_dt(t, P, w1, w2, d1, d2):
    p11, p22, p12, p13, p23 = P
    p33 = 1-p11-p22
    p21 = np.conj(p12)
    p31 = np.conj(p13)
    p32 = np.conj(p23)

    Pt = np.zeros(5, dtype='complex')

    Pt[0] = 1j * w1/2 * (p31 - p13) + g31*(p33) + gs*(p22 - p11)
    Pt[1] = 1j * w2/2 * (p32 - p23) + g32*(p33) - gs*(p22 - p11)
    Pt[2] = 1j * ((d1+d2)*p12 + w1/2*p32 - w2/2*p13) - gs*p12
    Pt[3] = 1j * (d1*p13 + w1/2*(p33-p11) - w2/2*(p12)) - gf*p13
    Pt[4] = 1j * (-1*d2*p23 - w1/2*p21 + w2/2*(p33-p22)) - gf*p23
    return Pt


def cw(t, P, w1, w2, d1, d2, fs=first_step, ms=max_step):
    if (fs > t/3):
        fs = t/3
    tspan=(0,t)
    on=(w1,w2)
    detuning=(d1,d2)
    return solve_ivp(dP_dt, tspan, P, args=on+detuning, first_step=fs, max_step=ms)


def free(t, P, d1, d2, fs=first_step, ms=max_step):
    tspan=(0,t)
    off=(0,0)
    detuning=(d1,d2)
    return solve_ivp(dP_dt, tspan, P, args=off+detuning, first_step=fs, max_step=ms)


def ramsey(tpulse, tfree, tmeas, P, w1, w2, d1, d2, fs=first_step, ms=max_step):

    # pulse
    result = cw(tpulse, P, w1, w2, d1, d2, fs=fs, ms=ms)

    p_curr = np.zeros(5, dtype='complex')
    for j in range(5):
        p_curr[j] = result.y[j, -1]

    # free evolution
    result2 = free(tfree, p_curr, d1, d2, fs=fs, ms=ms)
    p_curr2 = np.zeros(5, dtype='complex')
    for j in range(5):
        p_curr2[j] = result2.y[j, -1]

    result3 = cw(tmeas, p_curr2, w1, w2, d1, d2, fs=fs, ms=ms)

    return result3


def cw_resonance(t, P, I1, I2, d1_max, n=100, fs=first_step, ms=max_step):
    w1=rbc.w1(I1)
    w2=rbc.w2(I2)
    d2 = 0
    d1_lin = np.linspace(-1*d1_max, d1_max, num=n)
    p33_lin=np.ones(n, dtype='complex')
    for i in range(n):
        d1 = d1_lin[i]
        result = cw(t, P, w1, w2, d1, d2, fs=fs, ms=ms)
        p11f = result.y[0,-1]
        p22f = result.y[1,-1]
        p33f = 1-p11f-p22f
        p33_lin[i] = p33f
    return (d1_lin, p33_lin)


def ramsey_resonance(tpulse, tfree, tmeas, P, I1, I2, d1_max, n=100, fs=first_step, ms=max_step):
    w1=rbc.w1(I1)
    w2=rbc.w2(I2)
    d2 = 0
    d1_lin = np.linspace(-1*d1_max, d1_max, num=n)
    p33_lin=np.ones(n, dtype='complex')
    for i in range(n):
        d1 = d1_lin[i]
        result = ramsey(tpulse, tfree, tmeas, P, w1, w2, d1, d2, fs=fs, ms=ms)
        p11f = result.y[0,-1]
        p22f = result.y[1,-1]
        p33f = 1-p11f-p22f
        p33_lin[i] = p33f
    return (d1_lin, p33_lin)





