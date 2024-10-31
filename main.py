if __name__ == '__main__':

	# Loop through the specified days
	for t_days in t_days_values:
		print(f'\n\nResults for t_days={t_days}:\n')


		# These are each an array for the value t_days
		F_nu_1 = BreakCase(1, t_days, nu)
		F_tilde_2 = BreakCase(2, t_days, nu)
		F_tilde_3 = BreakCase(3, t_days, nu)

		# Multiplying each array f_nu_1, f_tilde_2, f_tilde_3 for F5

		F5 = F_nu_1 * F_tilde_2 * F_tilde_3


		# Multiplying each array f_nu_7, f_tilde_9, f_tilde_10, f_tilde_11 for F9

		F_nu_7 = BreakCase(7, t_days, nu)
		F_tilde_9 = BreakCase(9, t_days, nu)
		F_tilde_10 = BreakCase(10, t_days, nu)
		F_tilde_11 = BreakCase(11, t_days, nu)
		
		F9 = F_nu_7 * F_tilde_9 * F_tilde_10 * F_tilde_11


		matrix = np.column_stack((F5, F9, nu))



		np.savetxt("output.csv", matrix, delimiter=",", header="F5,F9,nu", comments='', fmt='%.50f')

