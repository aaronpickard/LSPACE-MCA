# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 21:43:25 2020

@author: Aaron Pickard
L'SPACE Mission Concept Academy Spring 2020
Team Ice-Nine
Mission Planning Equations
mission_math_functions.py

This file provides basic math functions to team participants in order to facilitate basic mission analysis tasks.

The functions are not grouped into classes (OOP) because this file is designed to be used in command line by STEM students who are not software developers.
The developer believes it will be easier to explain how to use the tool if the intended user is 
importing the file into IDLE and call functions than if they have to create and manipulate objects.

Presently, the file contains functions to:
1 Calculate Keplerian orbital parameters

WIP
1 Assess the size of the propulsive vehicle

Test functions not provided in this file (presently at all, but I'm hoping to change that eventually)
"""
import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Keplerian parameters
# SMAD pg 137
def get_semimajor_axis(r_a, r_p):
    a = (r_a + r_p)/2
    return a

def get_semimajor_axis(u, e):
    a = -u/(2*e)
    return a

def get_eccentricity(r_p, a):
    e = 1 - (r_p/a)
    return e

def get_eccentricity(r_a, a):
    e = (r_a/a) - 1
    return e

def get_inclination( h_z, h):
    i = math.acos(h_z/h)
    return i
    
def get_RAAN( n_x, n):
    omega = math.acos(n_x/n)
    return omega
    
def get_arg_of_pe(n, e):
    num = n * e
    den = abs(n) * abs(e)
    w = math.acos(num/den)
    return w
    
def get_true_anomaly(e, r):
    num = e * r
    den = abs(e) * abs(r)
    v = math.acos(num/den)
    return v
    
def get_r_p(a, e):
    r = a*(1-e)
    return r
    
def get_r_a(a, e):
    r = a*(1+e)
    return r
    
def get_period(a, u):
    p = 2 * math.pi * ((a**3)/u)**(0.5)
    return p
    
def get_frequency(a, u):
    w_0 = (u/(a**3))**(0.5)
    return w_0
 
# http://propagation.ece.gatech.edu/ECE6390/project/Fall2007/Crself.escentComm/teaself.m8 = 0/trajself.sim = 0.htmlself.
    
def get_prop_mass(dV, isp, mi, mp):
    # dV = g0*isp*ln((mi+mp)/mi)
    # mi is inert mass (structural mass + payload mass)
    # mp = propellant mass
    # prop = (e^(a/(b*c))-1)*(x+y)
    # prop = (e^(dV/(g0*isp))-1)*(inertmf+payloadmf)
    # mission dV = 2380m/s
    # lunar g0 = 1.625
    # total mission mass is 10 kg - if  mp + mi > 10 kg, report error
    # g0 = 1.625 # on the moon
    mprop = (math.exp(dV/(g0*isp))-1)*(mi+mp)
    return mprop
        
def get_inert_mass():
    pass
    
def get_payload_mass():
    pass
    
def get_inertmf(mi, mp,  mprop):
    self.imf = (mi + mp)/ (mi + mp + mprop)
    return imf
    
def get_propmf():
    pass
    
def get_max_mass_ratio(imf):
    mr_max = 1 / imf
    return self.mr_max
    
def get_mass_ratio(self):
    pass
    
def iterative_m_inert():
    pass
    
def iterative_m_propellant():
    pass
    
def iterate_master(max_counter, min_change):
    it_counter = 0
    it_mi_last = 0 # last calculated value for inert mass
    it_mpr_last = 0 # last calculated value for prop mass
    it_list = []
    while it_counter < max_counter:
        it_mi = iterative_m_inert()
        it_mpr = iterative_m_propellant()
        it_list = [it_mi, it_mpr]
        if abs(it_mi - it_mi_last) < min_change:
            return it_list
        if abs(it_mpr - it_mpr_last) <  min_change:
            return it_list
        it_mi_last = it_mi
        it_mpr_last = it_mpr
        it_counter+=1
    return it_list
         
