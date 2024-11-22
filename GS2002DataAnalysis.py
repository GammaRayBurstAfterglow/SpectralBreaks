import linecache
import numpy as np
import matplotlib.pyplot as plt

# Random iterator numbers selected to compare previous results from GS2002_text.f90 to grba.py
itm_numbers = [0, 27, 35, 45, 51]

# Fetches the line numbers based on each itm

# This function is very specific to the fort.810 output from GS2002_test.f90
def getLineNum(itm):

	line_start = 1 + itm * 644
	line_end = line_start + 641

	return line_start, line_end


# Opens fort.810 and names it HisData
with open("./WarrenCode/fort.810", "r") as HisData:

	# Seperates the data by line
	lines=HisData.readlines() 

	# Setting an empty array to define log10(Jy)
	y = []

	# Setting an empty array to define log10(Hz)
	x = []


	# iterates through the itm numbers
	for i in itm_numbers:
		# finds the start and end lines for each itm number
		line_start, line_end = getLineNum(i)

		# Creates new array in y, x for each itm_number
		y.append([])
		x.append([])

		for j in range(line_start, line_end):
	
	
			# y,x [-1] adds all of the data to the newly created arrays in the last for loop
			# Split is going through each section of the line, with the delimiter being empty space
	
			x[-1].append(float(lines[j].split()[3]))
			# Split 3 is finding log10(nuFnu) in fort.810

			y[-1].append(float(lines[j].split()[5]))
			# Split 5 is finding log10(Jy) in fort.810
	

	

			# The final output of y and x is a list of lists with each being a list of values for each itm_number

	HisData.close()


# plotting stuff
fig, ax = plt.subplots()
ax.plot(x[0], y[0], 'k-', linewidth=2, label="iteration 0")
ax.plot(x[1], y[1], 'r-', linewidth=2, label="iteration 27")
ax.plot(x[2], y[2], 'b-', linewidth=2, label="iteration 35")
ax.plot(x[3], y[3], 'g-', linewidth=2, label="iteration 45")
ax.plot(x[4], y[4], 'm-', linewidth=2, label="iteration 51")
plt.xscale('linear')
ax.set_yscale('linear')
ax.legend(loc="upper right")
ax.set_xlabel("Log(nu) [Hz]")
ax.set_ylabel("Log(F_nu) [Jy]")
plt.savefig("./output/GS2002.png")





GS2002_Data = np.empty((641, 6))

GS2002_Data[:,0] = x[0]
GS2002_Data[:,1] = y[0]
GS2002_Data[:,2] = y[1]
GS2002_Data[:,3] = y[2]
GS2002_Data[:,4] = y[3]
GS2002_Data[:,5] = y[4]




with open("./output/GS2002_Data.csv", "a") as f:
    np.savetxt(f, GS2002_Data, delimiter=",", header="Log(nu) [Hz], Log(F_nu) [Jy] i_tm=1, Log(F_nu) [Jy] i_tm=2, Log(F_nu) [Jy] i_tm=3, Log(F_nu) [Jy] i_tm=4, Log(F_nu) [Jy] i_tm=5", fmt='%.5f')
