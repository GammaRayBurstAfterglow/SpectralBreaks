import math

k = 0,2
#for now, only worried about k=0, but useful for future expansion

b = 1,2,3,7,9,10,11
# possible values of b. Method of selection is not considered, could be input() if desired

p = 2.5
# represents the spectral index of the electron distribution
# number of electrons with energy E is proportional to E^(-p)
# range between 2 - 4, 2.5 preferred value.

z = 1 
# cosmological red shift

epsilon_e = range(10**(-6), 0.4)
# fraction of the total energy density in electrons

epsilon_e_bar = (epsilon_e*(p-2)) / (p-1)

epsilon_B = range(10**(-9) 0.3)
# fraction total energy density in the magnetic field 

epsilon_p = 1 - (epsilon_e + epsilon_B)

n_0 = range(10**(-4), 10**(3))
# ambient density in units of particles/cm^3
# height of the density function for k = 0

E_52 = range(10**(-4), 10*(2))
# explosion energy of the GRB in units of 10^52 ergs

A_star = 
# height of the density function for k = 2

t_days = range(0.01, 100)
# time since explosion in units of days in the observer frame(?)

d_L28 = 
# luminosity distance in units of 10^28 cm
# ned wright cosmology calculator to convert z to d_L28

# Function for defining the break frequencies
def BreakCase(b,k):
	match(b,k):
		case(b = 1, k = 0):
			beta_1 = 2
			beta_2 = 1/3
			b = 1
			nu_b = nu_sa

			var1 = (pow((p - 1), 3/5)) / (pow((3*p + 2), 3/5))
			var2 = pow((1 + z), -1)
			var3 = pow(epsilon_e_bar, -1)
			var4 = pow(epsilon_B, 1/5)
			var5 = pow(n_0, 3/5)
			var6 = pow(E_52, 1/5)

			nu_1 = 1.24 * var1 * (10**9) * var2 * var3 * var4 * var5 * var6

		case(b = 2, k = 0):
			beta_1 = 1/3
			beta_2 = (1-p)/2
			nu_b = nu_m

			var1 = (p - 0.67) * (10**15)
			var2 = pow((1 + z), 1/2)
			var3 = pow(E_52, 1/2)
			var4 = pow(epsilon_e_bar, 2)
			var5 = pow(epsilon_B, 1/2)
			var6 = pow(t_days, -3/2)

			nu_2 = 3.73 * var1 * var2 * var3 * var4 * var5 * var6

		case(b = 3, k = 0):
			beta_1 = (1 - p)/2
			beta_2 = -p/2
			nu_b = nu_c

			var1 = (p - 0.46) * (10**13)
			var2 = pow(e, -1.16*p) # is the e here the function e or something else
			# if it is exponential function: math.exp(-1.16*p)
			var3 = pow((1 + z), -1/2)
			var4 = pow(epsilon_B, -3/2)
			var5 = pow(n_0, -1)
			var6 = pow(E_52, -1/2)
			var7 = pow(t_days, -1/2)
		              
			nu_3 = 6.37 * var1 * var2 * var3 * var4 * var5 * var6 * var7

		case(b = 7, k = 0):
			#dummy variables for writing nu_7, E.G. var_1 = var_7_1
			var1 = ( ((3p-1)**(8/5)) / ((3p+2)**(8/5)) )
			var2 = (1+z)**(-13/10)
			var3 = (epsilon_e_bar**(-8/5))
			var4 = (epsilon_B**(-2/5))
			var5 = (n_0**(3/10))
			var6 = (E_52**(-1/10))
			var7 = (t_day**(3/10))

			#Given by table 2 for b = 7
			beta_1 = 2
			beta_2 = 11/8

			#function given by nu_b of granot and sari
			nu_7 = 1.12 * var1 * (10**8) * var2 * var3 * var4 * var5 * var6 * var7

			#nu_b for the given case of b = 7, represents nu_ac
			return nu_7, beta_1, beta_2, s

		case(b = 9, k = 0):
			#dummy variables to shorten nu_9
			var1 = (p - 0.74)
			var2 = (1 + z)**(1/2)
			var3 = epsilon_e_bar**2
			var4 = epsilon_B**(1/2)
			var5 = E_52**(1/2)
			var6 = t_days**(-3/2)

			nu_9 = 3.94 * var1 * 10**15 * var2 * var3 * var4 * var5 * var6
			beta_1 = -1/2
			beta_2 = -p/2
			nu_b = nu_m

		case(b = 10, k = 0):
			#dummy variables to shorten nu_10
			var1 = (1+z)**(-1/2)
			var2 = epsilon_B**(6/5)
			var3 = n_0**(-1)
			var4 = E_52**(7/10)
			var5 = t_days**(-1/2)

			nu_10 = 1.32 * 10**10 * var1 * var2 * var3 * var4 * var5 

			beta_1 = 11/8
			beta_2 = 1/3
			nu_b = nu_sa

		case(b = 11, k = 0):
			beta_1 = 1/3
			beta_2 = -1/2

			#dummy variables to shorten nu_11
			var1 = pow((1+z), -1/2)
			var2 = pow(epsilon_B, -3/2)
			var3 = pow(n_0, -1)
			var4 = pow(E_52, -1/2)
			var5 = pow(t_days, -1/2)

			nu_11 = 5.86 * 10**12 * var1 * var2 * var3 * var4 * var5

		case _:
			print("Error: Invalid option for b, k")
			break
			#invalid option for b, k


# External flux density

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

#KP- I tried to add all the functions needed in one case break. I don't know if it's right or not. I am going to add it below as a comment:
'''
		case(beta_1 = 2, beta_2 = 1/3, s = 1.64 nu_b = nu_sa, b = 1, k = 0):

			#this case gives out the nu_b, F_tilde, F_ext, and F_nu for break 1 of k = 0

			var1 = (pow((p - 1), 3/5)) / (pow((3*p + 2), 3/5))
			var2 = pow((1 + z), -1)
			var3 = pow(epsilon_e_bar, -1)
			var4 = pow(epsilon_B, 1/5)
			var5 = pow(n_0, 3/5)
			var6 = pow(E_52, 1/5)

			nu_sa = 1.24 * var1 * (10**9) * var2 * var3 * var4 * var5 * var6
			F_tilde = (1 + pow((nu/nu_sa), s*(beta_1 - beta_2)))**(-1/s)
			
			#KP - I tried and included the ext flux density equation in this break case
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
'''

# Corresponding flux densities

def F(nu_b_ext, nu, nu_b, s, beta_1, beta_2):
    return nu_b_ext*(((nu/nu_b)**(-s*beta_1))+((nu/nu_b)**(-s*beta_2))**(-1/s))  
    #KP - Lily I think you are missing a parenthasis
	# implementation of eqn 1
	# function = flux at desired frequency, with passed argument of nu_b_ext, nu_b, s, beta_1, beta_2

def F_tilde(nu, nu_b, s, beta_1, beta_2):
    return (1+(nu/nu_b)**(s*(beta_1-beta_2)))**(-1/s)
	# implementation of eqn 4, with passed argument of nu_b, s, beta_1, beta_2

def F5(nu,s,beta_1,beta_2):
    return F(nu_1_ext,nu,nu_1,s,beta_1,beta_2)*F_tilde(nu,nu_2,s,beta_1,beta_2)*F_tilde(nu,nu_3,s,beta_1,beta_2)
	# implementation of eqn 5, with passed argument of nu, s, beta_1, beta_2

def F9(nu,s,beta_1,beta_2):
    return F(nu_7_ext,nu,nu_7,s,beta_1,beta_2)*F_tilde(nu,nu_10,s,beta_1,beta_2)*F_tilde(nu,nu_11,s,beta_1,beta_2)*F_tilde(nu,nu_9,s,beta_1,beta_2)
	# implementation of eqn 9, with passed argument of nu, s, beta_1, beta_2
