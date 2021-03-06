#-------------extract_gauge.py-------------------------------------------------#
#
# Purpose: This file takes a .mat file for a gauge field and turns it into 
#          something actually useable by GPUE / c++
#
#------------------------------------------------------------------------------#

import scipy.io

hbar = 1.05457148e-34

# function to extract gauge field from .mat file
def extract_field(filename, varname):
    mat = scipy.io.loadmat(filename)
    return mat[varname]

# function to output file in 2d format
def write_2d(var, outfile, fudge_factor):
    file = open(outfile,'w')
    for i in range(0,len(var)):
        for j in range(0,len(var)):
            file.write(str(var[i][j] * fudge_factor) + '\n')
    file.close()

# function to output 2d gauge field in 3d for GPU simulation
def write_3d(var, outfile, fudge_factor):
    file = open(outfile,'w')
    for i in range(0,len(var)):
        for j in range(0,len(var)):
            for k in range(0,len(var)):
                file.write(str(var[k][j] * fudge_factor) + '\n')
    file.close()

def write_constant(outfile, res):
    file = open(outfile,'w')
    for i in range(0,res*res*res):
        file.write(str(0) + '\n')
    file.close()

var = extract_field("Avec_broad_256.mat","avec")
#var = extract_field("Avec_128.mat","avec")

# These are the necessary parameters for when the BEC is 10x too large
#write_3d(var, "gauge_3d",  0.001)
#write_2d(var, "gauge_2d",  0.001)
write_3d(var, "gauge_3d",  1.0/3100.)
write_2d(var, "gauge_2d",  1.0/3100.)
#write_constant("const_256", 256)

# This is for the case where the BEC is the appropriate size
#write_3d(var, "gauge_3d",  1.0/2000000.0)
#write_2d(var, "gauge_2d",  1.0/2000000.0)
