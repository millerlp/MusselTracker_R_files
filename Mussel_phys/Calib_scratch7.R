# Calib_scratch7.R
# 
# Author: Luke Miller Aug 5, 2015
###############################################################################


setwd("D:/R/Mussel_phys")
source("D:/R/Mussel_phys/orientation_functions.R")
# Load calibration offset and scaling factors produced by my 'calibrate' 
# function and the RMagneto.dll
load('pololu_calib_20150805.RData')

accOffset = accCal$bias	# vector of 3 values
accScale = accCal$softiron	# 3x3 matrix
magOffset = magCal$bias
magScale = magCal$softiron

# Here is a simple "known" orientation of the sensor when the animal frame
# is assumed horizontal and north-facing. This could be replaced with a set
# of averaged accel/mag readings (with the calibration correction already 
# applied) from the calibration file when the animal was held horizontal and
# anterior/nose pointing north.
#accVecTest = c(0, 0, 17400)		# LSM303 +X axis pointing north, horizontal
#magVecTest = c(1370, 0, -1792) 	# LSM303 +X axis pointing north, horizontal
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

# LSM pointing slightly west of north, rolled -90 degrees to port, and rotated 
# so that it is
# facing down about -20 degrees from horizon, similar to what it might be on a
# mussel.
accVecTest = c(-6210, -15349, 1398)	# LSM303 pointing south, rolled +30 degrees
magVecTest = c(2171, 833, 462) # LSM303 pointing south, rolled +30 degrees


# Convert to n x 3 matrix
accVecTest = matrix(accVecTest,  nrow = 1, ncol = 3, byrow = TRUE)
magVecTest = matrix(magVecTest, nrow = 1, ncol = 3, byrow = TRUE)
# Process these new input matrices to use real world units and to arrange
# the axes so they conform to the AN4248 North-East-Down frame
procData = processRawLSM303D(accVecTest,magVecTest,accRange=2,magRange=4)
# Extract the accVec and magVec data
accVecG = procData$accVec
magVecT = procData$magVec

# Apply the calibration corrections
accVecCal = accVecG - accOffset
magVecCal = magVecT - magOffset
## Apply the scaling corrections
accVecCal  =  accScale %*% t(accVecCal)
magVecCal =  magScale %*% t(magVecCal)


########################################################
# Call calcWrotation function to estimate the W rotation matrix that is used
# to bring data from the sensor frame to the animal frame. The input data should
# be for an orientation where the animal (not necessarily the sensor) was 
# oriented horizontal, with its anterior/rostrum pointing to magnetic north. 
W = calcWrotation(accVecCal, magVecCal)


# This assumes that the animal/phone is currently horizontal, but the sensor is
# offset by some roll angle, r
# Eq 13
r = atan2(accVecCal[2],accVecCal[3])	# result in radians

# Calculate pitch angle as the inverse tangent of accelerometer x axis divided
# by the accel y + z axes and roll angle (in radians)
# Eq 15 uses the atan function so that the result is restricted in the range
# from -90 to +90. 
p = atan(-accVecCal[1]/ (accVecCal[2]*sin(r) + accVecCal[3]*cos(r)))
# Restrict pitch to range -90 to 90 degrees
if (p > 90 * pi/180) p = (180 * pi/180) - p
if (p < -90 * pi/180) p = (-180 * pi/180) - p
Rrot = RollMatrix(r)
Prot = PitchMatrix(p)
# Bring the magVec into the horizontal plane by applying the roll and pitch
# inverse rotation matrices
#magVecH = t(Rrot) %*% t(Prot) %*% magVec # multiply by inverses of Rrot, Prot
# Following AN4248, rotating the magnetometer to the horizontal involves 
# pre-multiplying by inverse of pitch, then inverse of roll.
magVecH = t(Prot) %*% t(Rrot) %*% magVecCal # multiply by inverses of Prot, Rrot
# Eq 22 solves for yaw (heading). When the magVec data are in the horizontal
# plane, the heading angle relative to magnetic north is just the angle between
# the y and x vectors, which we solve for using the inverse tangent function.
h = atan2(-magVecH[2],magVecH[1])	# result in radians

cat('Calibrated roll: ',round(r*180/pi),'degrees\n')
cat('Calibrated pitch: ',round(p*180/pi), 'degrees\n')
cat('Input heading: ',round(h*180/pi),'degrees\n')


#####################################################################
# Read in some test data
#df = read.csv("D:/temp/forward_pitch_north.csv", header = FALSE)
#df = read.csv("D:/temp/LSM303D_rotation_test3.csv", header = FALSE)
#df = read.csv("D:/temp/downward_pitch_360.csv", header = FALSE)
#df = read.csv("D:/temp/upward_pitch_360.csv", header = FALSE)
#df = read.csv("D:/temp/roll_360.csv", header = FALSE)
df = read.csv("D:/temp/rotate_heading_360.csv", header = FALSE)
#df = read.csv("D:/temp/pitch90_south_test.csv", header = FALSE)
#df = read.csv("D:/temp/pitch60_wobble.csv", header = FALSE)
#df = read.csv("D:/temp/pitch60_wobble2.csv", header = FALSE)

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
myhpr = hpr(accVecCal,magVecCal, W)
myhpr = myhpr * 180/pi
plothpr(myhpr)
# Plot un-rotated data (though still calibrated) as tiny dots
myhpr2 = hpr(accVecCal,magVecCal)
myhpr2 = myhpr2 * 180/pi
points(myhpr2[,3], pch = '.', col = 'black')
points(myhpr2[,2], pch = '.', col = 'green')
points(myhpr2[,1], pch = '.', col = 'red')



