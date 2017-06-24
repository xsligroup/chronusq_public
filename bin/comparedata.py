#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2022 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#


# This script compares the data from two CQ binfiles
import os,sys
import h5py
import numpy as np


fName1 = sys.argv[1] # Name of first file to compare
fName2 = sys.argv[2] # Name of second file to compare

dSetName = sys.argv[3] # Full path of H5 data set in both files
                       # e.g. /RESP/RESIDUE/EIGENVALUES

# Open H5 bins files
f1 = h5py.File(fName1)
f2 = h5py.File(fName2)


print ""
print "Comparing " + dSetName + " in"
print "  * File 1 = " + fName1 
print "  * File 2 = " + fName2 
print ""
print ""



# Check if files have dataset
f1HasData = dSetName in f1
f2HasData = dSetName in f2

if not f1HasData:
    print fName1 + " does not contain dataset " + dSetName
    sys.exit(1)

if not f2HasData:
    print fName1 + " does not contain dataset " + dSetName
    sys.exit(1)



# Load the data

f1Data = np.array(f1[dSetName])
f2Data = np.array(f2[dSetName])

if not (len(f1Data) == len(f2Data)):
    print "Lengths of dataset " + dSetName + " do not match in the two files"
    sys.exit(1)




diff = f1Data - f2Data
maxDiff = np.max(diff)
maxAbsDiff = np.max(np.abs(diff))

print "Max Diff     = ", maxDiff
print "Max Abs Diff = ", maxAbsDiff

print ""
