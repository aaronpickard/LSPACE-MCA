# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 21:43:25 2020

@author: Aaron Pickard
L'SPACE Mission Concept Academy Spring 2020
Team Ice-Nine
Mission Planning Equations
mission_math_functions.py

This file provides basic functions to team participants in order to facilitate basic mission analysis tasks.

The functions are not grouped into classes 
because this file is designed to be used in command line by STEM students who are not software developers.
The developer believes it will be easier to explain how to use the tool if the intended user is 
importing the file into IDLE and call functions than if they have to create and manipulate objects.

Presently, the file contains functions to:
1 Calculate Keplerian orbital parameters
2 Assess the size of the propulsive vehicle

WIP

TODO
1 Convert Keplerian orbital elements to Cartesian orbital elements
2 Convert Cartesian orbital elements to Keplerian orbital elements
3 Test completed functions (will happen in different Python file)
"""
import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

g0 = 9.81

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
    
# Source for vehicle sizing equations is notes from Summer 2018 
# lecture series @ NASA Langley Systems Analysis & Concepts Directorate
def m_prop_tsiolkovsky(dV, isp, m_inert, m_payload):
    # dV = g0*isp*ln((mi+mp)/mi)
    # mi is inert mass (structural mass + payload mass)
    # mp = propellant mass
    # prop = (e^(a/(b*c))-1)*(x+y)
    # prop = (e^(dV/(g0*isp))-1)*(inertmf+payloadmf)
    # mission dV = 2380m/s
    # lunar g0 = 1.625
    # total mission mass is 10 kg - if  mp + mi > 10 kg, report error
    # g0 = 1.625 # on the moon
    m_prop = (math.exp(dV/(g0*isp))-1)*(m_inert+m_payload)
    return m_prop
        
def m_ratio(dV, isp):
    m_ratio = math.exp(dv/g0*isp))
    return m_ratio

def m_inert_from_imf(m_payload, m_ratio, imf):
    # payload must be non-zero
    num = m_payload*(m_ratio-1)*(imf)
    den = 1 - m_ratio*imf
    m_inert = num/denom
    return m_inert
    
def m_propellant_from_imf(m_payload, m_ratio, imf): 
    # payload must be non-zero
    num = m_payload*(m_ratio-1)*(1-imf)
    den = 1 - m_ratio*imf
    m_prop = num/denom
    return m_prop
     
def imf(m_inert, m_prop):
    imf = m_inert/(m_inert+m_prop)
    return m_inert

def pmf(m_inert,  m_prop):
    pmf = m_prop/(m_inert+m_prop)
    return pmf
    
def m_ratio_max(imf):
    m_ratio_max = 1 / imf
    return m_ratio_max
     
def m_inert_iterative(imf, m_inert, m_prop):
    m_inert = imf*(m_inert+m_prop)
    return m_inert
    
def m_prop_iterative(m_ratio, m_inert, m_payload):
    # to be used if we know what the intended inert and payload masses are
    # probably to verify that we have enough propellant after allocating system masses
    m_prop = (m_ratio-1)*(m_inert+m_payload)
    return m_prop

def m_prop_iterative(m_ratio, m_non_reactive):
    # to be used if we haven't allocated system masses yet
    # earlier in the process
    m_prop = (m_ratio-1)*(m_non_reactive)
    return m_prop

def iterate_master(max_counter, min_change, m_inert, m_prop, m_ratio):
    # this function takes in as input variables to represent:
    # the max number of iterations in a loop, the minimum changes in the inert and propellant mass,
    # the initial inert mass, initial propellant mass, and the specified mass ratio for the design
    # this function returns as output a Python list of [inert mass, propellant mass, inert mass fraction]
    it_counter = 0
    it_mi_last = 0 # last calculated value for inert mass
    it_mpr_last = 0 # last calculated value for prop mass
    it_mi = m_inert
    it_mpr = m_prop
    it_mr = m_ratio
    it_imf = imf(it_mi, it_mpr)
    it_list = []
    while it_counter < max_counter:
        it_mi = m_inert_iterative(it_imf, it_mi, it_mpr)
        it_mpr = m_prop_iterative(m_ratio, it_mmi)
        it_imf = imf(it_mi, it_mpr)
        it_index = [it_mi, it_mpr, it_imf]
        if abs(it_mi - it_mi_last) < min_change:
            return it_list
        if abs(it_mpr - it_mpr_last) <  min_change:
            return it_list
        it_list = it_index
        it_mi_last = it_mi
        it_mpr_last = it_mpr
        it_counter+=1
    return it_list
         
