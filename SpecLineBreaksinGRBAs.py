import math

nu_sa
nu_m
nu_c


k = 0, 2
# radial dependence of the ambient medium
# density is proportional to r^(-k)
# k = 0, constant density, k = 2, 1/r^2

beta_1 =
beta_2 = 


p = 2.5
# represents the spectral index of the electron distribution
# number of electrons with energy E is proportional to E^(-p)
# range between 2 - 4, 2.5 preferred value.

z = 1 
# cosmological red shift


epsilon_e = 
# fraction of the total energy density in electrons
# between 10^-6 and 0.4

epsilon_e_bar = (epsilon_e*(p-2)) / (p-1)


epsilon_B = 
# fraction total energy density in the magnetic field 
# between 10^-9 and 0.3

epsilon_p = 1 - (epsilon_e + epsilon_B)


n_0 = 
# ambient density in units of particles/cm^3
# height of the density function for k = 0
# between 10^-4 to 10^3

E_52 = 
# explosion energy of the GRB in units of 10^52 ergs
# between 10^-4 and 10^2

A_star = 
# height of the density function for k = 2

t_day = 
# time since explosion in units of days in the observer frame(?)
# 0.01 and 100

d_L28 = 
# luminosity distance in units of 10^28 cm
# ned wright cosmology calculator to convert z to d_L28



#Krupa's Section
#When b = 1 beta_1 = 2 and beta_2 = 1/3
var_1_1 = (pow((p-1), 3/5)) / (pow((3p+2), 3/5))
nu_1 = 1.24 * var_1_1 * (10**9) * pow((1+z), -1) * pow(epsilon_e_bar, -1) * pow(epsilon_B, 1/5) * pow(n_0, 3/5) * pow(E_52, 1/5)
# b = 1, for k = 0

nu_2 =
# b = 2, for k = 0

nu_3 = 
# b = 3, for k = 0


#Damien's Section
nu_7 = 
# b = 7, for k = 0

nu_9 = 
# b = 9, for k = 0

nu_10 = 
# b = 10, for k = 0


#Niru's section
nu_11 = 
# b = 11, for k = 0

nu_1_ext = 
# b = 1, ext for k = 0

nu_7_ext = 
# b = 7, ext for k = 0

#Lily's Section
def F(): 
# function (nu) = flux at desired frequency, with passed argument of F_nu_b, nu_b, s, beta_1, beta_2
# implementation of eqn 1

def F_tilde():
# implementation of eqn 4
# function (nu) = flux at desired frequency, with passed argument of 

def (): 
# implementation of eqn 5	

def ():
# implementation of eqn 9
