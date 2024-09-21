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
#b = 1, for k = 0
#nu_b = nu_sa

var_1_1 = (pow((p - 1), 3/5)) / (pow((3*p + 2), 3/5))
var_1_2 = pow((1 + z), -1)
var_1_3 = pow(epsilon_e_bar, -1)
var_1_4 = pow(epsilon_B, 1/5)
var_1_5 = pow(n_0, 3/5)
var_1_6 = pow(E_52, 1/5)

nu_1 = 1.24 * var_1_1 * (10**9) * var_1_2 * var_1_3 * var_1_4 * var_1_5 * var_1_6

#When b = 2, for k = 0
#beta_1 = 1/3, beta_2 = (1-p)/2
#nu_b = nu_m

var_2_1 = (p - 0.67) * (10**15)
var_2_2 = pow(1 + z), 1/2)
var_2_3 = pow(E_52, 1/2)
var_2_4 = pow(epsilon_e_bar, 2)
var_2_5 = pow(epsilon_B, 1/2)
var_2_6 = pow(t_day, -3/2)

nu_2 = 3.71 * var_2_1 * var_2_2 * var_2_3 * var_2_4 * var_2_5 * var_2_6

#When b = 3, for k = 0
#beta_1 = (1 - p)/2, beta_2 = -p/2
#nu_b = nu_c

var_3_1 = (p - 0.46) * (10**13)
var_3_2 = pow(e, -1.16p) # is the e here the function e or something else
# if it is exponential function: math.exp(-1.16*p)
var_3_3 = pow((1 + z), -1/2)
var_3_4 = pow(epsilon_B, -3/2)
var_3_5 = pow(n_0, 1)
var_3_6 = pow(E_52, -1/2)
var_3_7 = pow(t_day, -1/2)
              
nu_3 = 6.37 * var_3_1 * var_3_2 * var_3_3 * var_3_4 * var_3_5 * var_3_6 * var_3_7

#Damien's Section


#I think this could be put into some sort of if statement / case system this could be more elegant
var_7_1 = ( ((3p-1)**(8/5)) / ((3p+2)**(8/5)) )
var_7_2 = (1+z)**(-13/10)
var_7_3 = (epsilon_e_bar**(-8/5))
var_7_4 = (epsilon_B**(-2/5))
var_7_5 = (n_0**(3/10))
var_7_6 = (E_52**(-1/10))
var_7_7 = (t_day**(3/10))

nu_7 = 1.12 * var_7_1 * (10**8) * var_7_2 * var_7_3 * var_7_4 * var_7_5 * var_7_6 * var_7_7
# b = 7, for k = 0
# beta_1 = 2
# beta_2 = 11/8
# nu_b = nu_ac

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
def F(nu_b_ext,nu,nu_b,s,beta_1,beta_2):
    return nu_b_ext*(((nu/nu_b)**(-s*beta_1))+((nu/nu_b)**(-s*beta_2))**(-1/s))  #KP - Lily I think you are missing a parenthasis
# implementation of eqn 1
# function = flux at desired frequency, with passed argument of nu_b_ext, nu_b, s, beta_1, beta_2

def F_tilde(nu,nu_b,s,beta_1,beta_2):
    return (1+(nu/nu_b)**(s*(beta_1-beta_2)))**(-1/s)
# implementation of eqn 4, with passed argument of nu_b, s, beta_1, beta_2

def F5(nu,s,beta_1,beta_2):
    return F(nu_1_ext,nu,nu_1,s,beta_1,beta_2)*F_tilde(nu,nu_2,s,beta_1,beta_2)*F_tilde(nu,nu_3,s,beta_1,beta_2)
# implementation of eqn 5, with passed argument of nu, s, beta_1, beta_2

def F9(nu,s,beta_1,beta_2):
    return F(nu_7_ext,nu,nu_7,s,beta_1,beta_2)*F_tilde(nu,nu_10,s,beta_1,beta_2)*F_tilde(nu,nu_11,s,beta_1,beta_2)*F_tilde(nu,nu_9,s,beta_1,beta_2)
# implementation of eqn 9, with passed argument of nu, s, beta_1, beta_2
