# this code is to be used in DWSIM Python Unit Operation

from DWSIM.Thermodynamics import *
import math
import System

# note: the format is [0, 1, 2, 3] = [H2, O2, H2O, N2]
# the position of each species is important!!!

# assigning input and output streams into variables
# format: inlet = i OR outlet = o, energy = e OR material = m, s = stream, number
enerin = ies1
enerout = oes1
inlet1 = ims1 # fuel
inlet2 = ims2 # air    
outlet1 = oms1 # fuel
outlet2 = oms2 # air

hydrogen = Flowsheet.SelectedCompounds['Hydrogen']
oxygen = Flowsheet.SelectedCompounds['Oxygen']
water = Flowsheet.SelectedCompounds['Water']
nitrogen = Flowsheet.SelectedCompounds['Nitrogen']

# constants
faraday = 96485.3365 # faraday constant in As/mol
R = 8.314 # gas constant in J/molK
MH2 = 2
MO2 = 32
MH2O = 18
MN2 = 28
MR = [2, 32, 18, 28]

# getting inlet properties
composition_f_in = inlet1.GetOverallComposition() # getting the mole fractions
molflow_f_in = inlet1.GetMolarFlow()
mass_f_in = inlet1.GetMassFlow()
T_f_in = inlet1.GetTemperature()
enthalpy_f_in = inlet1.GetMassEnthalpy() # in kJ/kg 
composition_a_in = inlet2.GetOverallComposition() # getting the mole fractions
molflow_a_in = inlet2.GetMolarFlow()
mass_a_in = inlet2.GetMassFlow()
T_a_in = inlet2.GetTemperature()
enthalpy_a_in = inlet2.GetMassEnthalpy() # in kJ/kg 
nocomps = len(composition_f_in)

# temperature used 
T_in = T_f_in

# creating empty arrays for calculation results
compounds_f_in = [] # mole
compounds_a_in = []
massfrac_f_in = [] # mass fractions
massfrac_a_in = []
massfrac_f_out = []
massfrac_a_out = []
composition_f_out = [] # mole fractions
comp_f_out_arr = []
composition_a_out = []
comp_a_out_arr = []

# converting mole to mass fraction
mass_i_list_f = []
mass_i_list_a = []
for i in range(nocomps):
    mass_i_f = composition_f_in[i] * MR[i]
    mass_i_list_f.append(mass_i_f)
    mass_i_a = composition_a_in[i] * MR[i]
    mass_i_list_a.append(mass_i_a)
massesnorm_f = sum(mass_i_list_f)
massesnorm_a = sum(mass_i_list_a)
for i in range(nocomps):
    massfrac_f = mass_i_list_f[i]/massesnorm_f
    massfrac_f_in.append(massfrac_f)
    massfrac_a = mass_i_list_a[i]/massesnorm_a
    massfrac_a_in.append(massfrac_a)

el_tr = (current/(4*faraday))
# according to 2H2 + O2 -> 2H2O 
v = [-2, -1, 2, 0] # O2 based
# [H2, O2, H2O, N2]
el_tr_list = [] # mol/s
for c in range(len(v)):
    tr_calc = el_tr*v[c] # electron transfer for each species (mol/s)
    el_tr_list.append(tr_calc)

# Mass Balances
# Total mass balance for F streams
mass_f_out = mass_f_in + el_tr*MO2
# Total mass balance for A streams 
mass_a_out = mass_a_in - el_tr*MO2 

# Species mass balances
massfrac_f_out = []
massfrac_a_out = []
for i in range(nocomps):
    if massfrac_f_in[i] != 0:
        msfrout_f = ((massfrac_f_in[i]*mass_f_in) + (el_tr_list[i]*MR[i]))/mass_f_out
        massfrac_f_out.append(msfrout_f)
    else:
        msfrout_f = 0
        massfrac_f_out.append(msfrout_f)
    if massfrac_a_in[i] != 0:
        msfrout_a = ((massfrac_a_in[i]*mass_a_in) + (el_tr_list[i]*MR[i]))/mass_a_out
        massfrac_a_out.append(msfrout_a)
    else:
        msfrout_a = 0
        massfrac_a_out.append(msfrout_a)

# Calculating mid mass fractions
midmassfr_f = []
midmassfr_a = []
for i in range(nocomps):
    mmff = (massfrac_f_in[i] + massfrac_f_out[i])/2
    midmassfr_f.append(mmff)
    mmfa = (massfrac_a_in[i] + massfrac_a_out[i])/2
    midmassfr_a.append(mmfa)
# Calculating mid mole fractions
midmolelst_f = []
midmolelst_a = []
midmolefr_f = []
midmolefr_a = []
for i in range(nocomps):
    mmlff = midmassfr_f[i]/MR[i]
    midmolelst_f.append(mmlff)
    mmlfa = midmassfr_a[i]/MR[i]
    midmolelst_a.append(mmlfa)
molesnorm_f = sum(midmolelst_f)
molesnorm_a = sum(midmolelst_a)
for i in range(nocomps):
    midmolef = midmolelst_f[i]/molesnorm_f
    midmolelst_f.append(midmolef)
    midmolea = midmolelst_a[i]/molesnorm_a
    midmolelst_a.append(midmolea)

# eqs for power input, based on fuel cell using mid mole fractions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
voltage = 1.23 + ((R*T_in)/(2*faraday))*math.log((midmolelst_f[0]*(math.sqrt(midmolelst_a[1])))/midmolelst_f[2])
power = (current*voltage*cells)/1000
enerin.EnergyFlow = power
power_used = power

# Convert to mole fractions
molelst_f = []
molelst_a = []
for i in range(nocomps):
    mmlff = massfrac_f_out[i]/MR[i]
    molelst_f.append(mmlff)
    mmlfa = massfrac_a_out[i]/MR[i]
    molelst_a.append(mmlfa)
moles_f = sum(molelst_f)
moles_a = sum(molelst_a)
for i in range(nocomps):
    molef = molelst_f[i]/moles_f
    composition_f_out.append(molef)
    molea = molelst_a[i]/moles_a
    composition_a_out.append(molea)
# Calculating mole flow out
molflow_f_out = mass_f_out/(composition_f_out[0]*MH2+composition_f_out[2]*MH2O)
molflow_a_out = mass_a_out/(composition_a_out[1]*MO2+composition_a_out[3]*MN2)
# Convert list of mole fractions into array
comp_f_out_arr = System.Array[float](composition_f_out)
comp_a_out_arr = System.Array[float](composition_a_out)

# Setting fuel output stream 
outlet1.SetTemperature(T_f_in)
outlet1.SetOverallComposition(comp_f_out_arr)
outlet1.SetMolarFlow(molflow_f_out)
outlet1.SetMassFlow(mass_f_out)
# Setting air output stream
outlet2.SetTemperature(T_a_in)
outlet2.SetOverallComposition(comp_a_out_arr)
outlet2.SetMolarFlow(molflow_a_out)
outlet2.SetMassFlow(mass_a_out)

# calculating energy balance

heat_out = -48.7 * el_tr_list[2] # kJ/s, where -48.7 kJ/mol is TdS part of H of formation of H2O

enthalpy_f_out = (enthalpy_f_in*mass_f_in)/mass_f_out + heat_out
enthalpy_a_out = (enthalpy_a_in*mass_a_in)/mass_a_out 

heatloss = mass_f_in*enthalpy_f_in + mass_a_in*enthalpy_a_in - mass_a_out*enthalpy_a_out - mass_f_out*enthalpy_f_out + power - power_used

enerout.EnergyFlow = heatloss
