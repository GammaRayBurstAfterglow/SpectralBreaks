k = 0,2
#for now, only worried about k=0, but useful for future expansion

b = 1,2,3,7,9,10,11
# possible values of b. Method of selection is not considered, could be input() if desired

match(b,k):
	case(b = 1, k = 0):

	case(b = 2, k = 0):

	case(b = 3, k = 0):

	case(b = 7, k = 0):
		#dummy variables for writing nu_7
		var_7_1 = ( ((3p-1)**(8/5)) / ((3p+2)**(8/5)) )
		var_7_2 = (1+z)**(-13/10)
		var_7_3 = (epsilon_e_bar**(-8/5))
		var_7_4 = (epsilon_B**(-2/5))
		var_7_5 = (n_0**(3/10))
		var_7_6 = (E_52**(-1/10))
		var_7_7 = (t_day**(3/10))

		#Given by table 2 for b = 7
		beta_1 = 2
		beta_2 = 11/8

		#I think there is an error here maybe. nu_b = nu_ac = nu_7
		nu_b = nu_ac

		#function given by nu_b of granot and sari
		nu_7 = 1.12 * var_7_1 * (10**8) * var_7_2 * var_7_3 * var_7_4 * var_7_5 * var_7_6 * var_7_7

		#nu_b for the given case of b = 7
		nu_b = nu_7

		return nu_b

	case(b = 9, k = 0):

	case(b = 10, k = 0):

	case(b = 11, k = 0):

	case _:
		print("Error: Invalid option for b, k")
		break
		#invalid option for b, k
