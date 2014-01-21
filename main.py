import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

pi = np.pi
sqrt = np.sqrt
cos = np.cos
sin = np.sin

def coulomb_henkel_h1_asym(l,rho):
    return (1j)**(-l)*np.exp(1j*rho)
def coulomb_henkel_h2_asym(l,rho):
    return (1j)**(l)*np.exp(-1j*rho)
def coulomb_henkel_h1_deriv_asym(l,rho,k):
    return k*(1j)**(-l+1)*np.exp(1j*rho)
def coulomb_henkel_h2_deriv_asym(l,rho,k):
    return -k*(1j)**(l+1)*np.exp(-1j*rho)

class single_channel():
    def __init__(self):
        self.E=1.0
        self.l=0
    def func(self,state, t):
        V=-61.1/(1.0+np.exp((t-1.2*np.power(10,1.0/3.0))/0.65))
        x, xdot = state
        return [xdot, (self.l*(self.l+1)/(t*t)+0.0478450*(V-self.E))*x]
    def return_k(self):
        k=np.sqrt(.047845*self.E)
        return k
    def set_E(self,E):
        self.E=E
    def set_l(self,l):
        self.l=l

if __name__=='__main__':
    plots = []
    fig = plt.figure(figsize=(6,7))
    Be11 = single_channel()
  
    angular_momenta = np.linspace(0,2,3)    
    energies = np.linspace(0.0001,4,400)
    for angular_momentum in angular_momenta:
        print "L = ",angular_momentum
        deltas = []    
        cross_section_el = []
        Be11.set_l(int(angular_momentum))
        for energy in energies:
            #SciPy Runge-Kutta
            print energy
            Be11.set_E(energy)
            domain = np.linspace(0.0000001, 500, 10000)
            state_init = [0.0, 1.0]
            state = integrate.odeint(Be11.func, state_init, domain)
            x, xdot = state.T

            #Calculate S_L, \delta_L, and \sigma_L
            a=8000 #radius index corresponding to radius=400
            radius = domain[a]
            k = Be11.return_k()
            rho = k*radius
            r_matrix_el = (1/radius)*(x[a]/xdot[a])    
            s_matrix_el = (coulomb_henkel_h2_asym(Be11.l,rho)-radius*r_matrix_el*coulomb_henkel_h2_deriv_asym(Be11.l,rho,k))/\
                (coulomb_henkel_h1_asym(Be11.l,rho)-radius*r_matrix_el*coulomb_henkel_h1_deriv_asym(Be11.l,rho,k))
            delta = (1/(2j))*np.log(s_matrix_el)
            delta = np.real(delta)
            crsxn = 4*np.pi/np.power(k,2)*np.power(sin(delta),2)
            deltas.append(delta)
            cross_section_el.append(crsxn)            
            #cross_section_el = deltas  #use this to plot deltas on the fly
        #MatPlotLib routines
        if Be11.l==0:            
            ll = plt.subplot(3,1,Be11.l+1)
            grid()
            ll.plot(energies, cross_section_el) 
            ll.set_yscale('log')
            #ll.yaxis.set_ticks(plt.LogLocator())            
            plt.minorticks_on()
            plots.append(ll)
            plt.ylabel(r"$\sigma_{el,L=%s}(E)$"% (Be11.l, ))
            plt.xlabel(r"E (MeV)")
            for tic in ll.xaxis.get_major_ticks():
                tic.tick1On = tic.tick2On = False
                tic.label1On = tic.label2On = False                    
        else:
            ll = plt.subplot(3,1,Be11.l+1,sharex=plots[Be11.l-1],sharey=plots[Be11.l-1])
            grid()
            ll.plot(energies, cross_section_el) 
            ll.set_yscale('log')
            #ll.yaxis.set_ticks(plt.LogLocator())            
            plt.minorticks_on()
            plt.ylabel(r"$\sigma_{el,L=%s}(E)$"% (Be11.l, ))            
            if Be11.l != angular_momenta[len(angular_momenta)-1]:
                plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
                for tic in ll.xaxis.get_major_ticks():
                    tic.tick1On = tic.tick2On = False
                    tic.label1On = tic.label2On = False        
            plt.xlabel(r"E (MeV)")
            ll.xaxis.set_major_locator(MaxNLocator(nbins=5,prune='lower'))
            ll.yaxis.set_major_locator(MaxNLocator(prune='upper'))
            plt.subplots_adjust(hspace=0)
            plots.append(ll)
            
    plt.savefig('/user/sullivan/public_html/research/matplotlib_output.pdf')




