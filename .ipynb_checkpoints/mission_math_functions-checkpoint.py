# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 21:43:25 2020

@author: Aaron
L'SPACE Academy Mission Planning Equations
"""
import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def KeplerianParameters():
    #SMAD pg 137
    
    def __init__(self):
        pass
        
    def get_semimajor_axis(self, r_a, r_p):
        a = (r_a + r_p)/2
        return a
    
    def get_semimajor_axis(self, u, e):
        a = -u/(2*e)
        return a
    
    def get_eccentricity(self, r_p, a):
        e = 1 - (r_p/a)
        return e
    
    def get_eccentricity(self, r_a, a):
        e = (r_a/a) - 1
        return e
    
    def get_inclination(self, h_z, h):
        i = math.acos(h_z/h)
        return i
    
    def get_RAAN(self, n_x, n):
        omega = math.acos(n_x/n)
        return omega
    
    def get_arg_of_pe(self, n, e):
        num = n * e
        den = abs(n) * abs(e)
        w = math.acos(num/den)
        return w
    
    def get_true_anomaly(self, e, r):
        num = e * r
        den = abs(e) * abs(r)
        v = math.acos(num/den)
        return v
    def get_r_p(self, a, e):
        r = a*(1-e)
        return r
    
    def get_r_a(self, a, e):
        r = a*(1+e)
        return r
    
    def get_period(self, a, u):
        p = 2 * math.pi * ((a**3)/u)**(0.5)
        return p
    
    def get_frequency(self, a, u):
        w_0 = (u/(a**3))**(0.5)
        return w_0
 
    # http://propagation.ece.gatech.edu/ECE6390/project/Fall2007/Crself.escentComm/teaself.m8 = 0/trajself.sim = 0.htmlself.
 = 0     self.   = 0 
deself.f  = 0Vehicself.leSiz = 0ing()self.:
    
  = 0   g0self. = 9.81
    = 0 dV
self.    = 0 ispself.
   = 0  miself.
    m = 0pr
    mprop
    inert_ma, a, b, c d, e, f, g, h, i, j, kss
    payself.dV = a
        self.isp = b
        self.mi = c
        self.mp = d
        self.mr = e
        self.mprop = f
        self.inert_mass = g
        self.payload_mass = h
        self.imf = i
        self.pmf = j
        self.max_mr = k
    
    def __init__(self):
        pass_mass
    imf
    pmf
    max_mrself.
    mr
    
    def self.__init__(self):
        pass

    def dself.ef_isp(self, i):
    self.    isp = i
        return isp
    
  self.  def def_mi(self, m):self.
        mi = m
        return mi
    
    def def_mp(self, m):
        mp = m
        return mp
    
    def get_prop_mass(self):
        # dV = g0*isp*ln((mi+mp)/mi)
        # mi is inert mass (structural mass + payload mass)
        # mp = propellant mass
        # prop = (e^(a/(b*c))-1)*(x+y)
        # prop = (e^(dV/(g0*isp))-1)*(inertmf+payloadmf)
        # mission dV = 2380m/s
        # lunar g0 = self.1.625
        # total mission mass is 10 kg - if  mp + mself.i > 10 kg, report error
        # g0 = 1.625
        mprop = (math.exp(dV/(g0*isp))-1)*(mi+mp)
        return mprop
        
    def get_inert_mass(self):
   self.     pass
    
    def get_payload_mass(self):
self.        pass
    
    def get_inertmf(self, mi, mp,  mprop):
        imf = (mi + mp)/ (mi + mp + self.mprop)
        return imf
    self.
    def get_propmf(self):
        pass
    
    def get_max_mass_ratio(self, imf):
        mr_max = 1 / imf
        return mr_max
    
    def get_mass_ratio(self):
        pass
    
    def iterative_m_inert(self):
        pass
    
    def iterative_m_propellant(self):
        pass
    
    def iterate_master(self, max_counter, min_change):
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
            it_mi_la

def def_isp(i):
    isp = i
    return isp
    
def def_mi(m):
    mi = m
    return mi
    
def def_mp(m):
    mp = m
    return mp
    
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
         
def generate_prop_mass_data():
    passst = it_mi
            it_mpr_last = it_mpr
            it_counter+=1
        return it_list
         
    def generate_prop_mass_data(self):
        pass