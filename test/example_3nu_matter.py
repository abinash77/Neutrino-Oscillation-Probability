# -*- coding: utf-8 -*-
r"""Run the matter 3nu example shown in README.md.

Runs the three-neutrino example of oscillations in matter
"""

from __future__ import print_function

import sys,os
sys.path.append('../src')

import numpy as np
import ROOT
from ROOT import gROOT, TCanvas, TFile,TH1D, TH2D,TH2F, TTree, TList, gStyle, TH2, Double, TObject

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

#c1 = TCanvas( 'c1', 'Probabilities', 700, 500 )

Emin = 0.1 #GeV
Emax = 100.1
binwidthE = 0.1
nbinsE = int((Emax-Emin)/binwidthE)
Lmin = 1 #km
Lmax = 13001
binwidthL = 1
nbinsL = int((Lmax-Lmin)/1)
print( "nbinsL=%6.5f \t nbinsE=%6.5f\n" %(nbinsL, nbinsE))
h2d_Pee = TH2D('h2d_Pee', 'Oscillation Probabilities', nbinsE, Emin, Emax, nbinsL, Lmin, Lmax) 
h2d_Pem = TH2D('h2d_Pem', 'Oscillation Probabilities', nbinsE, Emin, Emax, nbinsL, Lmin, Lmax)
h2d_Pmm = TH2D('h2d_Pmm', 'Oscillation Probabilities', nbinsE, Emin, Emax, nbinsL, Lmin, Lmax)
h2d_Pme = TH2D('h2d_Pme', 'Oscillation Probabilities', nbinsE, Emin, Emax, nbinsL, Lmin, Lmax)

for baseline in np.arange(3001, 4000 , binwidthL) :    #baseline units are in km
	for energy in np.arange(Emin, Emax, binwidthE):     # Neutrino energy [GeV]
	
		h_vacuum_energy_indep = \
    		hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(S12_NO_BF,
                                                                S23_NO_BF,
                                                                S13_NO_BF,
                                                                DCP_NO_BF,
                                                                D21_NO_BF,
                                                                D31_NO_BF,
							        							compute_matrix_multiplication=True)

# Units of VCC_EARTH_CRUST: [eV]
		if (int(baseline)%100 == 0 and int(energy)%10 == 0):
			print("baseline =%6.5f and energy = %6.5f\n" %(baseline, energy))
		h_matter = hamiltonians3nu.hamiltonian_3nu_matter(h_vacuum_energy_indep,
                                                    energy*1e9,
                                                    VCC_EARTH_CRUST)

		Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    		oscprob3nu.probabilities_3nu(h_matter, baseline*CONV_KM_TO_INV_EV)
		
		nubin = h2d_Pee.FindBin(energy,baseline)		
		h2d_Pee.SetBinContent(nubin,Pee)
		h2d_Pem.SetBinContent(nubin,Pem)
		h2d_Pme.SetBinContent(nubin,Pme)
		h2d_Pmm.SetBinContent(nubin,Pmm)
	#Energy loop ends here
#baseline loop ends here
f = TFile("probability_table_4.root","RECREATE")		
h2d_Pee.Write()
h2d_Pem.Write()
h2d_Pme.Write()
h2d_Pmm.Write()
f.Close()


