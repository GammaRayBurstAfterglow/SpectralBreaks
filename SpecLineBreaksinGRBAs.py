import math

###################################################
#
# Variable definitions and assigning constants
#
###################################################

# k = 0,2
#for now, only worried about k=0, but useful for future expansion

# b = 1,2,3,7,9,10,11
# possible values of b. Method of selection is not considered, could be input() if desired

# nu =
#range(10**(6), 2.418*(10**26))
# range between 10 MHz and 1 TeV
# convert energies into frequencies
#KP - Converted TeV to Hz. Someone please check it just to make sure.

p = 2.23
# represents the spectral index of the electron distribution
# number of electrons with energy E is proportional to E^(-p)
# range between 2 - 4, 2.5 preferred value.

z = 1 
# cosmological red shift

epsilon_e = 0.1
#range(10**(-6), 0.4)
# fraction of the total energy density in electrons

epsilon_e_bar = (epsilon_e*(p-2)) / (p-1)

epsilon_B = 0.01
#range(10**(-9) 0.3)
# fraction total energy density in the magnetic field 

epsilon_p = 1 - (epsilon_e + epsilon_B)

n_0 = 1
#range(10**(-4), 10**(3))
# ambient density in units of particles/cm^3
# height of the density function for k = 0

E_52 = 1
#range(10**(-4), 10*(2))
# explosion energy of the GRB in units of 10^52 ergs

# A_star = 
# height of the density function for k = 2
# Only needed for k = 2, ignoring for now

t_days = [1.157*(10**-4), 5.64*(10**-2), 5.585*(10**-1), 3.495, 13.831] 
#Convereted from seconds to days
# 3.083E+03 Convert from seconds to days
#range(0.01, 100)
# time since explosion in units of days in the observer frame(?)

d_L28 =  2.095 # (10**28)cm
# 6787.5 Mpc convert to 10^28
# luminosity distance in units of 10^28 cm
# ned wright cosmology calculator to convert z to d_L28
# Implement directly at some point if z != 1
# https://www.astro.ucla.edu/~wright/CosmoCalc.html

###############################################################
#
# Corresponding flux densities
#
# F and F_tilde will be implemented directly into BreakCase():
#
##############################################################
'''
def F(nu_b_ext, nu, nu_b, s, beta_1, beta_2):
	return nu_b_ext*(((nu/nu_b)**(-s*beta_1))+((nu/nu_b)**(-s*beta_2))**(-1/s)) 
	# implementation of eqn 1
	# function = flux at desired frequency, with passed argument of nu_b_ext, nu_b, s, beta_1, beta_2

def F_tilde(nu, nu_b, s, beta_1, beta_2):
	return (1+(nu/nu_b)**(s*(beta_1-beta_2)))**(-1/s)
	# implementation of eqn 4, with passed argument of nu_b, s, beta_1, beta_2
'''


# Needs to be rewritten to align the new BreakCase():
def F5(nu,s,beta_1,beta_2):
	return F(nu_1_ext,nu,nu_1,s,beta_1,beta_2)*F_tilde(nu,nu_2,s,beta_1,beta_2)*F_tilde(nu,nu_3,s,beta_1,beta_2)
	# implementation of eqn 5, with passed argument of nu, s, beta_1, beta_2

def F9(nu,s,beta_1,beta_2):
	return F(nu_7_ext,nu,nu_7,s,beta_1,beta_2)*F_tilde(nu,nu_10,s,beta_1,beta_2)*F_tilde(nu,nu_11,s,beta_1,beta_2)*F_tilde(nu,nu_9,s,beta_1,beta_2)
	# implementation of eqn 9, with passed argument of nu, s, beta_1, beta_2

##############################################################
#
# Main function (Refer to line 328 of GS2002_test.f90)
#
# At the end of the fucntion, return F5, F9?
#
##############################################################


for i in loop():



def BreakCase(b, k, beta_1, beta_2, s, nu_b, nu):
	match(b, k, beta_1, beta_2, s, nu_b, nu):
	# Values for b, k, beta_1, beta_2, s, and nu_b are given by Granot and Sari



		# Break 1 for k = 0
		case(b = 1, k = 0, beta_1 = 2, beta_2 = 1/3, s = 1.64 nu_b = nu_sa):

			# Dummy variables for calculating nu_b | nu_b = nu_
			# Dummy variation notation will follow as such for all break cases:
			# var_i = var_x_i for b = x and i being the iteration, x is understood by the case of BreakCase():
			var1 = (pow((p - 1), 3/5)) / (pow((3*p + 2), 3/5))
			var2 = pow((1 + z), -1)
			var3 = pow(epsilon_e_bar, -1)
			var4 = pow(epsilon_B, 1/5)
			var5 = pow(n_0, 3/5)
			var6 = pow(E_52, 1/5)

			nu_sa = 1.24 * var1 * (10**9) * var2 * var3 * var4 * var5 * var6
			

			# The chunk below solves for F_nu equation 1
			var_1_ext_1 = (pow((p - 1), 6/5)) / ((3*p - 1)*(pow((3*p + 2), 1/5)))
			var_1_ext_2 = pow((1 + z), 1/2)
			var_1_ext_3 = pow(epsilon_e_bar, -1)
			var_1_ext_4 = pow(epsilon_B, 2/5)
			var_1_ext_5 = pow(n_0, 7/10)
			var_1_ext_6 = pow(E_52, 9/10)
			var_1_ext_7 = pow(t_days, 1/2)
			var_1_ext_8 = pow(d_L28, -2)

			nu_1_ext = 0.647 * var_1_ext_1 * var_1_ext_2 * var_1_ext_3 * var_1_ext_4 * var_1_ext_5 * var_1_ext_6 * var_1_ext_7 * var_1_ext_8

			F_nu = nu_1_ext * (pow((nu/nu_b), -s*beta_1) + pow((nu/nu_b), -s*beta_2))**(-1/s)

			return F_nu, beta_1, beta_2 



		# Break 2 for k = 0
		case(b = 2, k = 0, beta_1 = 1/3, beta_2 = ((1-p)/2), s = (1.84-(0.40*p)), nu_b = nu_m):
			
			# Dummy variables for calculating nu_2
			var1 = (p - 0.67) * (10**15)
			var2 = pow((1 + z), 1/2)
			var3 = pow(E_52, 1/2)
			var4 = pow(epsilon_e_bar, 2)
			var5 = pow(epsilon_B, 1/2)
			var6 = pow(t_days, -3/2)

			nu_2 = 3.73 * var1 * var2 * var3 * var4 * var5 * var6

			F_tilde_2 = (1 + pow((nu/nu_sa), s*(beta_1 - beta_2)))**(-1/s)

			return F_tilde_2



		# Break 3 for k = 0
		case(b = 3, k = 0, beta_1 = ((1-p)/2), beta_2 = (-p/2), s = (1.15-(0.06*p)), nu_b = nu_c):

			# Dummy variables for calculating nu_3
			var1 = (p - 0.46) * (10**13)
			var2 = math.exp(-1.16*p)
			var3 = pow((1 + z), -1/2)
			var4 = pow(epsilon_B, -3/2)
			var5 = pow(n_0, -1)
			var6 = pow(E_52, -1/2)
			var7 = pow(t_days, -1/2)
					  
			nu_3 = 6.37 * var1 * var2 * var3 * var4 * var5 * var6 * var7

			F_tilde_3 = (1 + pow((nu/nu_sa), s*(beta_1 - beta_2)))**(-1/s)

			return F_tilde_3



		# DK - Still missing s, and nu_b, need to refer to Granot and Sari for these later
		# Break 7 for k = 0
		case(b = 7, k = 0, beta_1 = 2, beta_2 = 11/8, s = , nu_b = ):

			#dummy variables for calculating nu_7
			var1 = ( ((3p-1)**(8/5)) / ((3p+2)**(8/5)) )
			var2 = (1+z)**(-13/10)
			var3 = (epsilon_e_bar**(-8/5))
			var4 = (epsilon_B**(-2/5))
			var5 = (n_0**(3/10))
			var6 = (E_52**(-1/10))
			var7 = (t_day**(3/10))


			# Function given by nu_b of Granot and Sari
			nu_7 = 1.12 * var1 * (10**8) * var2 * var3 * var4 * var5 * var6 * var7



		# DK - Still missing s, need to refer to Granot and Sari for this
		# Break 9 for k = 0
		case(b = 9, k = 0, beta_1 = -1/2, beta_2 = -p/2, s = , nu_b = nu_m):

			#dummy variables calculating nu_9
			var1 = (p - 0.74)
			var2 = (1 + z)**(1/2)
			var3 = epsilon_e_bar**2
			var4 = epsilon_B**(1/2)
			var5 = E_52**(1/2)
			var6 = t_days**(-3/2)

			# Function given by nu_b of Granot and Sari
			nu_9 = 3.94 * var1 * 10**15 * var2 * var3 * var4 * var5 * var6



		# DK - Still missing s, need to refer to Granot and Sari for this
		# Break 10 for k = 0
		case(b = 10, k = 0, beta_1 = 11/8, beta_2 = 1/3, s = , nu_b = nu_sa):

			# Dummy variables to calculate nu_10
			var1 = (1+z)**(-1/2)
			var2 = epsilon_B**(6/5)
			var3 = n_0**(-1)
			var4 = E_52**(7/10)
			var5 = t_days**(-1/2)

			# Function given by nu_b of Granot and Sari
			nu_10 = 1.32 * 10**10 * var1 * var2 * var3 * var4 * var5 



		# DK - Still missing s, and nu_b, need to refer to Granot and Sari for these later
		# Break 11 for k = 0
		case(b = 11, k = 0, beta_1 = 1/3, beta_2 = -1/2, s = , nu_b =):

			#dummy variables to shorten nu_11
			var1 = pow((1+z), -1/2)
			var2 = pow(epsilon_B, -3/2)
			var3 = pow(n_0, -1)
			var4 = pow(E_52, -1/2)
			var5 = pow(t_days, -1/2)

			# Function given by nu_b of Granot and Sari
			nu_11 = 5.86 * 10**12 * var1 * var2 * var3 * var4 * var_11_5



		# Else case, I.E. invalid option for b, k, breaks code and outputs error message
		case _:
			print("Error: Invalid option for b, k")
			break


#KP - Tried to implement equation 5 here
def F1(F_nu, F_tilde_2, F_tilde_3):
	return F_nu * F_tilde_2 * F_tilde_3



# External Flux density with incorrect (Should be a flux density instead of nu_b_ext)

#nu_1_ext: 
# b = 1, ext for k = 0

#dummy variables to shorten nu_1_ext
#var_1_ext_1 = (pow((p - 1), 6/5)) / ((3*p - 1)*(pow((3*p + 2), 1/5)))
#var_1_ext_2 = pow((1 + z), 1/2)
#var_1_ext_3 = pow(epsilon_e_bar, -1)
#var_1_ext_4 = pow(epsilon_B, 2/5)
#var_1_ext_5 = pow(n_0, 7/10)
#var_1_ext_6 = pow(E_52, 9/10)
#var_1_ext_7 = pow(t_days, 1/2)
#var_1_ext_8 = pow(d_L28, -2)

#nu_1_ext = 0.647 * var_1_ext_1 * var_1_ext_2 * var_1_ext_3 * var_1_ext_4 * var_1_ext_5 * var_1_ext_6 * var_1_ext_7 * var_1_ext_8

#nu_7_ext: 
# b = 7, ext for k = 0

#dummy variables to shorten nu_7_ext
#var_7_ext_1 = ( (pow(3p-1, 11/5)) / (pow(3p+2, 11/5)) )
#var_7_ext_2 = pow(1+z, -1/10)
#var_7_ext_3 = pow(epsilon_e_bar, -4/5))
#var_7_ext_4= pow(epsilon_B, -4/5)
#var_7_ext_5 = pow(n_0, 1/10)
#var_7_ext_6 = pow(E_52, 3/10)
#var_7_ext_7 = pow(t_days, 11/10)
#var_7_ext_8 = pow(d_L28, -2)

#nu_7_ext = 5.27 * pow(10, -3) * var_7_ext_1 * var_7_ext_2 * var_7_ext_3 * var_7_ext_4 * var_7_ext_5 * var_7_ext_6 * var_7_ext_7 * var_7_ext_8
