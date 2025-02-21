# This look 2 directories deep from the starting directory and strips out all of the elapsed time per job (Expected to be ~20s)
# It takes this data, sorts it, creates a graph with it
# It also searches for non-empty error files and outputs the data

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# file.enswith("*_out_*.dat") wasn't working right, so I imported this weird module
# I was missing the first * in "*_out_*.dat" This took 90 minutes to figure out :)
from fnmatch import fnmatch

# Function to search through the training data and extract the time elapsed per indivdual job (~280,000?)
def process_files(directory):
    elapsed_times = []
    for root, dirs, files in os.walk(directory):
        #print(f"Checking directory: {root}")
        #print(f"Files found: {files}")

        # Ensure to only check two levels deep to save computing time
        if root[len(directory):].count(os.sep) == 2:
            for file in files:

                #print(f"Processing file: {file}")

                # Only searching the output files
                if fnmatch(file, "*_out_*.dat"):
                
                    file_path = os.path.join(root, file)
                    #print(f"Processing file: {file_path}")
                    with open(file_path, 'r') as f:
                        for line in f:
                            # Each output files has "Emission calculated!" in the same line as the time elapsed
                            if "Emission calculated!" in line:
                                time_value = line.split('time = ')[-1].split(' sec')[0]
                                elapsed_times.append(time_value)
                                #print(time_value)
    return elapsed_times


# This is looking for non-empty error files
# This function can probably be implemented into process_files(), most of the start is the same
def errors(directory):
    error_files = []
    for root, dirs, files in os.walk(directory):
        # Ensure to only check two levels deep to save computing time
        if root[len(directory):].count(os.sep) == 2:
            for file in files:
                if fnmatch(file, "*.err*"):
                    file_path = os.path.join(root, file)
                    with open(file_path, 'r') as f:
                        error_output = f.read()
                        # only appends the file if it's not empty
                        if error_output.strip(): 
                            error_files.append((file, error_output))
    return error_files




def main():
    ########################
    #   Processing  Data   #
    ########################

    # Not sure what this is on the lab pc, just putting .
    # I.E. make sure to run this script in the right directory
    base_directory = os.path.expanduser('.')

    elapsed_times = process_files(base_directory)

    # Save all of the elapsed times to output.txt
    with open('output.txt', 'w') as f:
        for time in elapsed_times:
            f.write(f"{time}\n")

    # Sort the file content
    with open('output.txt', 'r') as f:
        sorted_times = sorted(f.readlines())

    # Overwrite the file with the sorted times
    with open('output.txt', 'w') as f:
        f.writelines(sorted_times)


    ########################
    #    Error Analysis    #
    ########################

    # Looking for error files
    error_output = errors(base_directory)

    # outputs all of the nonempty error files and the content of their files into errors.txt
    with open("errors.txt", 'w') as f:
        for file, content in error_output:
            f.write(f"File: {file}")
            f.write(f"Content: {content}\n")


    ########################
    # Using Processed Data #
    ########################

    # Still need to implement find -iname '*.err*' and output it to a text file

    # Loads the outputted file from TimeElapsed.sh
    data = np.loadtxt('output.txt')

    # Creates a histogram from all of the data in output.txt
    plt.hist(data, edgecolor='black')
    plt.xlabel("Time (s)")
    plt.ylabel("Frequency")

    # saves the graph in the current working directory
    plt.savefig('histogram.png')

    # Loads the elapsed times into a dataframe for head() and tail()
    df = pd.read_csv("output.txt")

    # Saves the outliar values to head.txt and tail.txt
    np.savetxt(r'head.txt', df.head(3).values, '%f')
    np.savetxt(r'tail.txt', df.tail(3).values, '%f')


if __name__ == "__main__":
    main()
