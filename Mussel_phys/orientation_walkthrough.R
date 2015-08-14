# orientation_walkthrough.R
# 
# Author: Luke Miller Aug 3, 2015
###############################################################################
setwd('./Mussel_phys/')
# Load accel/magnetometer into memory
source("D:/R/Mussel_phys/orientation_functions.R") 


# Step 1, generate a set of offset and scale corrections for the accelerometer
# and magnetometer. Ideally you have a set of calibration data with data from
# a variety of orientations all around a sphere, which will reveal whether 
# certain axes display an offset, or a non-linearity, in response to gravity 
# or magnetic fields. 

# Example calibration data for the Pololu board
fname = "D:/R/Mussel_phys/test_rotations/pololu_calib.csv"
dfo = read.csv(fname, header = FALSE)

# Extract raw acceleration and magnetometer vectors
accVec = dfo[,c(1,2,3)]
magVec = dfo[,c(4,5,6)]
# Default settings for the Pololu board
accRange = 2	# +/- 2g
magRange = 4	# +/- 4gauss

# Use the processRawLSM303D function to convert the raw LSM303D counts into
# real-world units, also with the axes inverted as needed to conform to the
# AN4248 North-East-Down navigation reference frame.
procData = processRawLSM303D(accVec,magVec,accRange,magRange)
# Extract the accVec and magVec data
accVecG = procData$accVec
magVecT = procData$magVec

# Local norms. Gravity will always be ~1000 milli-g, but magnetic field total
# norm value depends on your location on earth and the time. See the site
# http://www.ngdc.noaa.gov/geomag-web/#igrfwmm for your local conditions
accNorm = 1000 # 1000 milli-g local acceleration (gravity) norm
magNorm = 47894.9	# nanoTesla local magnetic field total norm

##########################################
# Use calibrate function to produce offset and scaling factors
# These may also be referred to as bias or hard iron (for offset), or 
# 'soft iron' (for scaling)
accCal = calibrate(accVecG[,1],accVecG[,2],accVecG[,3], accNorm)
magCal = calibrate(magVecT[,1],magVecT[,2],magVecT[,3], magNorm)

# Extract the bias and scale
accOffset = accCal$bias		# This will be a single value
accScale = accCal$softiron	# This will be a 3x3 matrix
magOffset = magCal$bias		# This will be a single value
magScale = magCal$softiron	# This will be a 3x3 matrix
# With the appropriate offset and scaling values, new input data from the sensor
# can be calibrated.

# If we're still working with the calibration data, the calibrate function
# also returned the calibrated accVec and magVec values. 
accVecCal = accCal$Calibrated
magVecCal = magCal$Calibrated


## Alternatively, for testing you can load the bigger pololu board calibration
# routine I did on the gimbal jib. This takes a while to generate the
# values, so it's quicker to just reload the values from a RData file. 
# Load calibration offset and scaling factors produced by my 'calibrate' 
# function and the RMagneto.dll
#load('D:/R/Mussel_phys/test_rotations/pololu_calib_20150805.RData')
#accOffset = accCal$bias	# vector of 3 values
#accScale = accCal$softiron	# 3x3 matrix
#magOffset = magCal$bias
#magScale = magCal$softiron


################################################################################
# Part 2 - determining a rotation matrix between sensor frame and animal frame

# In Part 1 we used a set of slow rotations to determine how the accel/mag axes
# were offset from zero and non-linear, so that we can bring any reading back
# to a value that would be read with an ideal sensor in the same orientation.

# We can now use readings from the calibration file when the animal was in a 
# known orientation (horizontal, facing north) to estimate the rotation matrix
# needed to rotate the sensor frame into line with the animal frame. We use
# horizontal, north, because this allows us to estimate rotations to the simple
# case where accelerometer X+Y axes should read 0 ideally, and where the 
# magnetometer X axis should be maximal and the Y-axis is at zero when 
# horizontal

# Here are some simple "known" orientations of the sensor when the animal frame
# is assumed horizontal and north-facing. This could be replaced with a set
# of averaged accel/mag readings (with the calibration corrections already 
# applied) from the calibration file when the animal was held horizontal and
# anterior/nose pointing north.
accVecTest = c(0, 0, 17400)		# LSM303 +X axis pointing north, horizontal
magVecTest = c(1370, 0, -1792) 	# LSM303 +X axis pointing north, horizontal
#
#accVecTest = c(9078, 91, 14718)		# LSM303 +X axis pointing north, pitch ~+30deg
#magVecTest = c(-137, -67, -2369) 	# LSM303 +X axis pointing north, pitch ~+30deg
#
#accVecTest = c(0, 0, 17400)		# LSM303 +X axis pointing west, horizontal
#magVecTest = c(0, -1995, -1778) # LSM303 +X axis pointing west, horizontal

#accVecTest = c(145, 17089, 927)	# LSM303 pointing west, rolled +90 to right
#magVecTest = c(-10, -2327, 2256) # I think LSM303 pointing west, rolled +90 to right

#accVecTest = c(-66, 63, 17324)	# LSM303 pointing south, horizontal
#magVecTest = c(-2612, 102, -1499) # LSM303 pointing south, horizontal

#accVecTest = c(-58, 7822, 15559)	# LSM303 pointing south, rolled +30 degrees
#magVecTest = c(-2561, -815, -1494) # LSM303 pointing south, rolled +30 degrees

# Convert to n x 3 matrix
accVecTest = matrix(accVecTest,  nrow = 1, ncol = 3, byrow = TRUE)
magVecTest = matrix(magVecTest, nrow = 1, ncol = 3, byrow = TRUE)
# Process these new input matrices to use real world units and to arrange
# the axes so they conform to the AN4248 North-East-Down frame
procData = processRawLSM303D(accVecTest,magVecTest,accRange,magRange)
# Extract the accVec and magVec data
accVecG = procData$accVec
magVecT = procData$magVec

# Apply the calibration corrections
accVecCal = accVecG - accOffset
magVecCal = magVecT - magOffset
## Apply the scaling corrections
accVecCal  =  accScale %*% t(accVecCal)
magVecCal =  magScale %*% t(magVecCal)


# Call calcWrotation function to estimate the W rotation matrix that is used
# to bring data from the sensor frame to the animal frame. The input data should
# be for an orientation where the animal (not necessarily the sensor) was 
# oriented horizontal, with its anterior/rostrum pointing to magnetic north. 
W = calcWrotation(accVecCal, magVecCal)

################################################################################
################################################################################
# Part 3
# With the offset and scale factors in hand, along with the W rotation matrix,
# we are now in a position to take in new data with the animal in any position
# and derive the animal's orientation in the earth frame (i.e. relative to 
# horizontal and magnetic north). 
# Raw data from the sensor should be read in, possibly filtered, and then the
# calibration factors (offset and scale) should be applied. After that, the 
# rotation matrix W should be applied to the values to rotate the sensor values
# into the animal frame. 
# From there, for each time point, you can calculate pitch and roll, then use 
# those values to re-orient the magnetometer data into the horizontal plane 
# (aka gimbal the magnetometer), and derive the heading relative to magnetic
# north. 

# Read in some test data
#df = read.csv("D:/R/Mussel_phys/test_rotations/forward_pitch_north.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/LSM303D_rotation_test3.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/downward_pitch_360.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/upward_pitch_360.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/roll_360.csv", header = FALSE)
df = read.csv("D:/R/Mussel_phys/test_rotations/rotate_heading_360.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/pitch90_south_test.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/pitch60_wobble.csv", header = FALSE)
#df = read.csv("D:/R/Mussel_phys/test_rotations/pitch60_wobble2.csv", header = FALSE)

# First run the raw data through the processRawLSM303D function to make sure
# the columns are arranged properly and values are converted to real world
# units
procData = processRawLSM303D(df[,1:3],df[,4:6],accRange = 2, magRange = 4)

accVec = procData$accVec
magVec = procData$magVec

# Apply the offset+scale factors derived earlier from the calibration data
accVecCal = accVec
magVecCal = magVec
for (i in 1:nrow(accVec)){
	# Subtract bias value for each axis
	accVecCal[i,1] = accVec[i,1] - accOffset[1]
	accVecCal[i,2] = accVec[i,2] - accOffset[2]
	accVecCal[i,3] = accVec[i,3] - accOffset[3]
	# Pre-multiply the bias-corrected values by the softiron matrix
	accVecCal[i,1:3] =  accScale %*% accVecCal[i,1:3]
	# Subtract bias value for each axis
	magVecCal[i,1] = magVec[i,1] - magOffset[1]
	magVecCal[i,2] = magVec[i,2] - magOffset[2]
	magVecCal[i,3] = magVec[i,3] - magOffset[3]
	# Pre-multiply the bias-corrected values by the softiron matrix
	magVecCal[i,1:3] =  magScale %*% magVecCal[i,1:3]
}

# Estimate heading, pitch, and roll at each time point by supplying the sensor
# data and the W rotation matrix to rotate values from sensor frame to animal
# frame. Results are in radians
#myhpr = hpr(accVecCal,magVecCal)
myhpr = hpr(accVecCal,magVecCal, W)
myhpr = myhpr * 180/pi
plothpr(myhpr)
# Plot un-rotated data (though still calibrated) as tiny dots
myhpr2 = hpr(accVecCal,magVecCal)
myhpr2 = myhpr2 * 180/pi
points(myhpr2[,3], pch = '.', col = 'black')
points(myhpr2[,2], pch = '.', col = 'green')
points(myhpr2[,1], pch = '.', col = 'red')


