#from popstar import synthetic, evolution, atmospheres, reddening, ifmr
#from popstar.imf import imf, multiplicity
#import os, sys, math
import numpy as np
from gcwork import orbits
from astropy.table import Table, Column, MaskedColumn
import matplotlib.pyplot as plt
from astropy.constants import G
from astropy import units as u
from copy import deepcopy
from random import choice

def a_to_P(mass, a):
    """
    Goes from semimajor axis in AU to period in years
    
    mass - float or array-like
    primary object mass in Msun
    
    a - float or array-like
    semimajor axis in AU
    """
    
    G_units = G.to("AU3/(M_sun*year2)").value
    period = (a**3*4*(np.pi**2)/G_units/mass)**(1/2)
    return period

def add_positions(ss):
    """
    Adds x and y positions randomly in a box of length and width 40000 AU for each system.
    
    ss - star system table without positions
    """
    ss_temp = deepcopy(ss)
    
    ss_temp.add_column( Column(np.zeros(len(ss_temp), dtype=float), name='x', description='AU') )
    ss_temp.add_column( Column(np.zeros(len(ss_temp), dtype=float), name='y', description='AU') )
    
    sign_x = np.array([choice([-1,1]) for i in range(len(ss_temp))])
    sign_y = np.array([choice([-1,1]) for i in range(len(ss_temp))])
    ss_temp['x'] = sign_x*20000*np.random.rand(len(ss_temp))
    ss_temp['y'] = sign_y*20000*np.random.rand(len(ss_temp))
    
    return ss_temp

def add_mult_positions(companions, ss_pos, logAge):
    """
    Adds x and y positions of multiple companions by transforming keplerian parameters to xyz in AU
    using Siyao's code and random initial times. Then adding them to the random posiiton of the primary object.
    
    ss_pos - star system table with positions added with add_positions()
    
    companion - companion table without positions
    
    logAge - float or int
    log of age of cluster with age in years
    """
    companions_temp = deepcopy(companions)
    
    companions_temp.add_column( Column(np.zeros(len(companions_temp), dtype=float), name='x', description='AU') )
    companions_temp.add_column( Column(np.zeros(len(companions_temp), dtype=float), name='y', description='AU') )
        
    orb = orbits.Orbit()
    for i in companions_temp:
        orb.w = i['omega'] #degrees
        orb.o = i['Omega'] #degrees
        orb.i = i['i'] #degrees
        orb.e = i['e'] #between 0 and 1
        orb.p = a_to_P(ss_pos[i['system_idx']]['mass'],10**i['log_a']) #year
        orb.t0 = (10**logAge)*np.random.rand() #year
        
        (r, v, a) = orb.kep2xyz(np.array([10**logAge]),mass=ss_pos[i['system_idx']]['mass'], dist=1.)
               
        x = r[0][0]
        y = r[0][1]
        
        #putting positions relative to primary object
        i['x'] = ss_pos[i['system_idx']]['x'] + x #AU
        i['y'] = ss_pos[i['system_idx']]['y'] + y #AU
    
    
    return companions_temp

def plot_projected_cluster(ss_pos, companions_pos):
    """
    Plots projected cluster with lines between companions and primary stars
    
    Takes companion and star system tables modified with add_mult_positions() and add_positions() respectively
    """
    plt.figure(figsize=(10,10))
    plt.plot(ss_pos['x'], ss_pos['y'],linestyle='none',marker='.' )
    plt.plot(companions_pos['x'], companions_pos['y'],linestyle='none',marker='.' )
    
    #makes lines between companion and primary star
    for i in companions_pos:
        plt.plot([i['x'], ss_pos[i['system_idx']]['x']],[i['y'], ss_pos[i['system_idx']]['y']],color='grey',linewidth=1)
        
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.show()
    
    return

def plot_companion_orbit(ss, companions_pos, logAge, t0 = None, system = None):
    """
    Plots the orbit of one system assuming the primary object is at (0,0). By default random companion and initial time.
    
    ss - star system table (does not matter if it has positions or not)
    
    companion_pos - companion table with positions added with add_mult_positions()
    
    logAge - float or int
    log of age of cluster with age in years
    
    t0 - float or int
    initial time of creation of the system in years (by default random)
    
    system - int
    index of desired companion in companion_pos table (by default random)
    """
    
    if system == None:
        system = np.random.randint(len(companions_pos))
    if t0 == None:
        t0 = (10**logAge)*np.random.rand() #year

    companion = companions_pos[system]
    
    orb = orbits.Orbit()
    orb.w = companion['omega'] #degrees
    orb.o = companion['Omega'] #degrees
    orb.i = companion['i'] #degrees
    orb.e = companion['e'] #between 0 and 1
    orb.p = a_to_P(ss[companion['system_idx']]['mass'],10**companion['log_a']) #year
    orb.t0 = t0 #year
        
    (r, v, a) = orb.kep2xyz(np.linspace(1, 10**logAge), mass=ss[companion['system_idx']]['mass'], dist=1.)
        
    for i in r:
        plt.plot(i[0],i[1], marker='.', color = 'blue')
        if i[0] == r[-1][0]:
            plt.plot(i[0],i[1], marker='*', color = 'gold', markersize=15)
            
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.show()
    
    return