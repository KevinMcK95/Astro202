# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 13:52:46 2018

Author: Kevin McKinnon
Title: Astro 202 Assignment 1
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

#constants (in cgs)
h = np.float128(6.62607004 * 10**-27)
c = np.float128(2.99792458 * 10**10)
k = np.float128(1.38064852 * 10**-16)
sig = (2 * np.pi**5 * k**4)/float(15 * c**2 * h**3)

def freqPlanck(nu,T):
    #Planck's law with frequency
    coeff = (2 * h * nu**3) / (c**2)
    exponent = (h * nu) / (k * T)
    return coeff * (np.exp(exponent)-1)**-1


def wavPlanck(lam,T):
    #Planck's law with wavelength
    coeff = (2 * h * c**2 / lam**5)
    exponent = (h * c) / (lam * k * T)
    return coeff * (np.exp(exponent)-1)**-1


def planckPlotter(T,form,figNum=0):
    #plots the Planck function in the given form at the given temperature
    #Providing a figNum is useful when trying to put multiple T on the same plot

    if figNum == 0:
        plt.figure()
    else:
        plt.figure(num=figNum)
    plt.title('Planck Function')
    plt.ylabel('Brightness [erg$\cdot$s$^{-1}$$\cdot$cm$^{-2}$$\cdot$sr$^{-1}$$\cdot$Hz$^{-1}$]')
    
    if form == 'wavelength':
        top = 10*10**-3
        x = np.arange(0.000000000001,top,top/10000.0) #lambda values from 0 to inifinity
        y = wavPlanck(x,T)
        plt.xlabel('Wavelength [cm]')
    elif form == 'frequency':
        x = np.logspace(5,15,num=1000)
        y = freqPlanck(x,T)
        plt.xlabel('Frequency [Hz]')
        plt.loglog()
    else:
        print "ERROR: The form must be either 'frequency' or 'wavelength'"
        return
    plt.plot(x,y,label='T = %d'%T)
    if figNum != 0:
        plt.legend(loc='best')
    plt.tight_layout()
    return
        

def planckIntegrator(T):
    #choose a step size and integrate using a trapezoidal method of numerical integration
    #stop integrating at a given step size when the change in the integrated value is less than 0.001%
    #chose a new step size 1/10th the previous and repeat
    #return an integral answer when the step size changes produce a change of less than 10^-10%
    #or it reaches 10000 times the lambda_max
    
    step = np.float128(0.01)

    integrals = [0] #initially, integral value is 0
    intChange = 100 #initiator for the integral change as a function of step size
    
    while (intChange >= 0.1):
        integral = np.float128(0)
        pastIntegral = np.float128(0)
        lam = np.float128(10**-15)  #initial lambda value
        percentChange = 100         #initial value
        bVals = np.float128(0),np.float128(0)
        passedPeak = False #calue to check if peak value has been passed yet
        peakVal = np.inf
        while (not passedPeak) or ((percentChange > 10**-10) and (lam < 10000 * peakVal)): 
        #need to be passed the peak wavelength and also the percent change to be less than 10^-10% to stop
            bVals = wavPlanck(lam,T),wavPlanck(lam+step,T)
            if (bVals[1] - bVals[0] < 0) and not passedPeak:
                peakVal = lam
                passedPeak = True
            if passedPeak and bVals[1] == 0:
                break
            integral += (bVals[1]+bVals[0])*0.5 * step
            if pastIntegral != 0:
                percentChange = abs(integral - pastIntegral)/float(pastIntegral) * 100
            pastIntegral = integral
            lam += step
        intChange = abs(integral - integrals[-1])/float(integrals[-1]) * 100
        integrals.append(integral)
        step = step * 0.1 #reduce step size by a factor of 10
    return integrals[-1]



def derivs(x,func,side=10**-10):
    #returns the func value and derivative at x
    fs = func(x-side),func(x),func(x+side) #B values
    fPrime = (fs[2] - fs[0])/float(2*side)
    return fs[1],fPrime


def wavPeak():
    #use Newton's method formula on dB/d lam = 0 
    #Checks to see the percent change from the previous root answer
    #and stops after the changes are less than 0.001%

    oldX = 10 #first guess
    side = 10**-10
        
    percentChange = 100 #initial value
        
    def planckDeriv(x):
        #this is the simplified result of dB/d lam = 0 and changing variables to x = (h * c)/ (lam * k * T)
        return 5*(np.exp(x)-1) - x*np.exp(x)

    while percentChange > 0.001:
        fatx, fprimeatx = derivs(oldX,planckDeriv,side=side)
        newX = oldX - (fatx)/float(fprimeatx)
        if oldX != 0:
            percentChange = abs(newX - oldX)/float(oldX) * 100
        oldX = newX
    return h*c / (oldX * k) #return the value that T * lam_max is equal to
    
def freqPeak():
    #use Newton's method formula on dB/d nu = 0 
    #Checks to see the percent change from the previous root answer
    #and stops after the changes are less than 0.001%
   
    oldX = 5 #first guess, want it to be away from 0 otherwise we might get that solution
    side = 10**-10
        
    percentChange = 100 #initial value
    def planckDeriv(x):
        #this is the simplified result of dB/d nu = 0 and changing variables to x = (h * nu)/ (k * T)
        return 3*(np.exp(x)-1) - x*np.exp(x)
        
    while percentChange > 0.001:
        fatx, fprimeatx = derivs(oldX,planckDeriv,side=side)
        newX = oldX - (fatx)/float(fprimeatx)
        if oldX != 0:
            percentChange = abs(newX - oldX)/float(oldX) * 100
        oldX = newX
    return oldX * k / h #return the value that nu_max / T is equal to


def tempDerivPeak():
    #use Newton's method formula on d(dB/dT)/d nu = 0 
    #Checks to see the percent change from the previous root answer
    #and stops after the changes are less than 0.001%

    oldX = 5 #first guess, want it to be away from 0 otherwise we might get that solution
    side = 10**-10
        
    percentChange = 100 #initial value
    
    def planckDoubleDeriv(x):
        #this is the simplified result of d(dB/dT)/d nu = 0 and changing variables to x = (h * nu)/ (k * T)
        return (2*x**3+8*x**2)*(np.exp(x)-1)-4*x**3*np.exp(x)
        
    while percentChange > 0.001:
        fatx, fprimeatx = derivs(oldX,planckDoubleDeriv,side=side)
        newX = oldX - (fatx)/float(fprimeatx)
        if oldX != 0:
            percentChange = abs(newX - oldX)/float(oldX) * 100
        oldX = newX
    return oldX * k / h #return the value that nu_max / T is equal to
   
   
def runAll():
    #produces results for T=500, 800, and 1000 Kelvin
    plt.close('all')
    for temp in [500,800,1000]:
        planckPlotter(temp,'wavelength',figNum=1)
        planckPlotter(temp,'frequency',figNum=2)
        
        planckIntegration = planckIntegrator(temp)
        print '\nThe Planck function integral at T=%d is %.5f'%(temp,planckIntegration)
        trueInt = sig * temp**4/np.pi
        diff = trueInt/planckIntegration * 100
        print '(The true integral value is %.5f, which is %f percent of the numerical solution)\n'%(trueInt,diff)
    
    print 'lambda_max * T = %0.5f for the Planck function\n'%(wavPeak())
    print 'nu_max / T = %0.5f for the Planck function\n'%(freqPeak())
    print 'nu_max / T = %0.5f for the temperature derivative of the Planck function\n'%(tempDerivPeak())
    return
    
#runAll()
