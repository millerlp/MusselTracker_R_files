# orientation_functions.R
# 
# Author: Luke Miller Aug 3, 2015
###############################################################################


# Load the Rmagneto.dll library (for Windows) to gain access to the 
# calibrate function. This file needs to be in the current working directory.
dyn.load('Rmagneto.dll')

###########################################
# Simple function to unload the dll, mainly useful during 
# troubleshooting of the C functions that underlay Rmagneto
detach.Rmagneto = function(){
	dyn.unload('Rmagneto.dll')
}

################################################################################
#	Function	processRawLSM303D
#	Use this function to convert raw LSM303D accleration and magnetometer data
# 	into real-world units (milli-G, where 1G = 9.8 m/s^-2, and nanoTesla) 

processRawLSM303D = function(accVec,magVec,accRange,magRange){
	# Start by making sure input accel and mag data are matrices
	if (class(accVec) != 'matrix'){
		if (class(accVec) == 'data.frame') {
			accVec = as.matrix(accVec)
			magVec = as.matrix(magVec)
		} else if (class(accVec == 'numeric' & length(accVec) == 3)){
			accVec = matrix(accVec, nrow = 1, byrow = TRUE)
			magVec = matrix(magVec, nrow = 1, byrow = TRUE)
		} else {
			warning("acceleration or magnetometer data can't be converted to nx3 matrix")
			break
		}	
	}
# Apply conversions to the raw outputs of the accelerometer and magnetometer
# to convert them to real units. Conversion factors are taken from the 
# LSM303D data sheet
	if (accRange == 2) {
		accVec = accVec * 0.061 # +/- 2 g
	} else if (accRange == 4) {
		accVec = accVec * 0.122 # +/- 4 g
	} else if (accRange == 6) {
		accVec = accVec * 0.183 # +/- 6 g
	} else if (accRange == 8) {
		accVec = accVec * 0.244 # +/- 8 g
	} else if (accRange == 16) {
		accVec = accVec * 0.732	# +/- 16 g
	}
	
	if (magRange == 2) {
		magVec = magVec * 0.08	# +/- 2 gauss
	} else if (magRange == 4) {
		magVec = magVec * 0.160	# +/- 4 gauss
	} else if (magRange == 8) {
		magVec = magVec * 0.320	# +/- 8 gauss
	} else if (magRange == 12) {
		magVec = magVec * 0.479	# +/- 12 gauss
	}
# Convert the magVec values from milligauss to nanoTesla
	magVec = magVec * 100
	
# Based on the North-East-Down reference frame outlined in AN4248, 
# rearrange the LSM303D accel and magnetometer axes to conform to the frame. 
	accVec[,1] = -1 * accVec[,1] # invert accelerometer x axis
	magVec[,2] = -1 * magVec[,2] # invert mag y-axis so east is positive
	magVec[,3] = -1 * magVec[,3] # invert mag z-axis so down is positive
	
	result = list(accVec = accVec, magVec = magVec)
	result	# return result as a list
}



###############################################################################
# Function calibrate
# Supply 3 vectors of X, Y, and Z data from a magnetometer or accelerometer,
# along with the local field total norm (for magnetometer) or the local 
# gravitational norm (which is always just 1000 milli-g). 
# 
# Output: a list containing the entries:
#	bias: 3 coefficients, to be subtracted from X, Y, and Z values respectively
#	softiron: a 3x3 matrix of coefficients used to scale the X, Y, Z axes 
# 	Calibrated: The original input vectors X, Y, Z, with the bias and softiron
#				calibrations applied. A 3 x n matrix. 
#
calibrate = function(X,Y,Z,localnorm){
	# X, Y, Z should be the raw data, re-scaled to magnetic units (nano-Tesla)
	# or acceleration units (millli-Galileo, milli-g). 
	# localnorm is the local Total norm of the magnetic field for your site
	# during the collection dates: http://www.ngdc.noaa.gov/geomag-web/#igrfwmm
	# or the norm of gravity, which we'll usually assume is simply 1000 milli-g
	# (i.e. 1 g). 
	
	# Generate the other input arguments and output vectors
	nrows = length(X) 	
	b = numeric(3)
	Xcorr = numeric(3)
	Ycorr = numeric(3)
	Zcorr = numeric(3)
	
	# The C function 'calibrate' is expecting the following arguments:
	# void calibrate(int *X, int *Y, int *Z, int *nrows, double *localnorm, 
	# 					double *b, double *Xcorr, double *Ycorr, double *Zcorr)
	
	
	# Call the C routine 'calibrate' from the Rmagneto.dll library
	# This needs to be have been loaded previously using a function like
	# dyn.load('Rmagneto.dll')
	res = .C('calibrate', as.double(X), as.double(Y), as.double(Z),
			as.integer(nrows), as.double(localnorm), as.double(b),
			as.double(Xcorr), as.double(Ycorr), as.double(Zcorr))
	# Keep the correction values
	results = list(bias = res[[6]], 
			softiron = as.matrix(cbind(res[[7]],res[[8]],res[[9]])))
	names(results$bias) = c('X','Y','Z')
	colnames(results$softiron) = c('X','Y','Z')
	
	# By convention we subtract the bias corrections from the original values
	Xoff = X - results$bias[1]
	Yoff = Y - results$bias[2]
	Zoff = Z - results$bias[3]
	# Then apply the scaling factors to the bias-corrected values
	Xcal = Xoff * results$softiron[1,1] + 
			Yoff * results$softiron[2,1] + 
			Zoff * results$softiron[3,1]
	Ycal = Xoff * results$softiron[1,2] +
			Yoff * results$softiron[2,2] + 
			Zoff * results$softiron[3,2]
	Zcal = Xoff * results$softiron[1,3] + 
			Yoff * results$softiron[2,3] + 
			Zoff * results$softiron[3,3]
	results$Calibrated = as.matrix(cbind(Xcal,Ycal,Zcal))
	colnames(results$Calibrated) = c('X','Y','Z')
	
	results	# return result at end of function. 
}

################################################################################
# Function 		RollMatrix
# Calculate pitch rotation matrix based off of AN4248 Eq 5
# Inputs:
#	roll	a n x 1 vector of roll angles, in radians.
# Outputs:
#	R		a 3 x 3 or 3 x 3 x n array used to rotate the roll angle to 
#			horizontal. 
RollMatrix = function(roll){
	if (length(roll) == 1){
		R = matrix(c( 1, 			0, 			 	0,
						0, 	cos(roll),		sin(roll),
						0, 	-sin(roll), 	cos(roll)),
				nrow = 3, byrow = TRUE)		
	} else if (length(roll) > 1){
		R = array(NA, dim = c(3,3,length(roll)))
		for (i in 1:length(roll)){
			R[,,i] = matrix(c( 1, 				0, 				 0,
							0, 		cos(roll[i]),	sin(roll[i]),
							0, 		-sin(roll[i]), 	cos(roll[i])),
					nrow = 3, byrow = TRUE)	
		}
	}
	R	# return roll rotation matrix
}

################################################################################
# Function 		PitchMatrix
# Calculate Pitch rotation matrix using AN4248 Eq 6
# Inputs:
# 	pitch	a n x 1 vector of pitch values, in radians
# Outputs:
#	P		a 3 x 3 matrix or 3 x 3 x n matrix to rotate input pitch values to
# 			the horizontal 
PitchMatrix = function(pitch){
	if (length(pitch) == 1){
		P = matrix(c( 	cos(pitch),   		0,  -sin(pitch),
						0,         1,     		  0,
						sin(pitch),   		0,  cos(pitch)),
				nrow = 3, byrow = TRUE)	
	} else if (length(pitch) > 1) {
		P = array(NA, dim = c(3,3,length(pitch)))
		# Need to create a 3rd dimension to house multiple
		# 3x3 rotation matrices
		for (i in 1:length(pitch)){
			P[,,i] = matrix(c( 	cos(pitch[i]),   		0,  -sin(pitch[i]),
							0,         1,     		  0,
							sin(pitch[i]),   		0,  cos(pitch[i])),
					nrow = 3, byrow = TRUE)	
		}
	}
	P	# return matrix or array
}
################################################################################
# Function		YawMatrix
# Calculate a rotation matrix to rotate a raw yaw value to magnetic north
# Based on AN4248 Eq 7
# Inputs:
#	yaw		a n x 1 vector of yaw values, in radians
# Outputs:
#	Y		a 3 x 3 or 3 x 3 x n array to be used as a rotation matrix to bring
# 			yaw values magnetic north. 
YawMatrix = function(yaw){
	if (length(yaw) == 1){
		Y = matrix(c(	cos(yaw),	sin(yaw),		0,
						-sin(yaw),	 cos(yaw),		0,
								0,			0,		1),
				nrow = 3, byrow = TRUE)
	} else if (length(yaw) > 1){
		Y = array(NA, dim = c(3,3,length(yaw)))
		Y = matrix(c(	cos(yaw[i]),	sin(yaw[i]),		0,
						-sin(yaw[i]),	 cos(yaw[i]),		0,
									0,			0,			1),
				nrow = 3, byrow = TRUE)
	}
	Y	# return y rotation matrix
}

################################################################################
# Function		calcWrotation
# A function to generate the sensor-to-frame rotation matrix W 
# Inputs:
# 	accVec	a n x 3 matrix of acceleration data from a time period when the
# 			animal is in a known position (horizontal, anterior facing to 
#			magnetic north). May be a single set of xyz values, or multiple
# 			values which will be averaged to generate average xyz values. Values
#			should already have offset and scale calibrations applied. Units
#			should be in milli-G.
#	magVec	a n x 3 matrix of magnetometer data from a time period when the
# 			animal is in a known position (horizontal, anterior facing to 
#			magnetic north). May be a single set of xyz values, or multiple
# 			values which will be averaged to generate average xyz values. Values
#			should already have offset and scale calibrations applied. Units
# 			should be in nanoTesla.
# Outputs:
# 	W		a 3 x 3 rotation matrix that can be used to rotate sensor data
#			from the sensor frame to the animal frame
calcWrotation = function(accVec, magVec){
	# Deal with multi-row matrices that will need to be averaged
	if (class(accVec) == 'matrix') {
		if (sum(dim(accVec) == c(3,1)) > 1){
			# If fed a column vector, transpose it to a row vector
			accVec = t(accVec)
			magVec = t(magVec)
		}
		accVec = colMeans(accVec)
		magVec = colMeans(magVec)
	} else if (class(accVec) == 'data.frame'){
		accVec = matrix(accVec, ncol = 3, nrow = nrow(accVec), byrow = TRUE)
		accVec = colMeans(accVec)
		magVec = matrix(magVec, ncol = 3, nrow = nrow(magVec), byrow = TRUE)
		magVec = colMeans(magVec)
	} else {
		warning("Couldn't coerce data into n x 3 matrix")
		break
	}
	# Normalize values before calculating pitch, roll, heading of new data
	accVecNorm = accVec / norm(accVec,"2")
	magVecNorm = magVec / norm(magVec,"2")
	
	# Now calculate the pitch and roll values of the sensor, again assuming that
	# the animal frame is horizontal and north-facing
	# Calculate roll angle as the inverse tangent of the animal/phone frame 
	# y axis and z axis (roll is registered as a change on y and z axes)
	# Eq 13
	r = atan2(accVecNorm[2],accVecNorm[3])	# result in radians
	# Calculate pitch angle as the inverse tangent of accelerometer x axis 
	# divided by the accel y + z axes and roll angle (in radians)
	# Eq 15 uses the atan function so that the result is restricted in the range
	# from -90 to +90. 
	p = atan(-accVecNorm[1]/ (accVecNorm[2]*sin(r) + accVecNorm[3]*cos(r)))

	# Restrict pitch to range -90 to 90 degrees
	p[p>(90*180/pi)] = (180*pi/180) - p[p>(90*180/pi)]
	p[p<(-90*180/pi)] = (-180*pi/180) - p[p<(-90*180/pi)]
	

	# With pitch and roll
	# Now that you know pitch and roll angles, the magnetometer readings can be
	# rotated to the horizontal plane, correcting for the animal/phone 
	# orientation. AN4248 Eq 17 allows you to use the roll and pitch rotation 
	# matrices to convert the magnetic vector in the animal/phone frame into a 
	# magnetic vector if the animal/phone were horizontal. In AN4248 this is the
	# change from the Bp vector to the Bf (flat) vector in Eq 19. 
	Rrot = RollMatrix(r)
	Prot = PitchMatrix(p)
	# Multiply by inverses of Prot, Rrot to bring magVecNorm to horizontal
	magVecH = t(Prot) %*% t(Rrot) %*% magVecNorm 
	
	# Eq 22 solves for yaw (heading). When the magVec data are in the horizontal
	# plane, the heading angle relative to magnetic north is just the angle 
	# between the y and x vectors, which we solve for using the inverse tangent 
	# function.
	hEst = atan2(-magVecH[2],magVecH[1])	# result in radians
	# Calculate the rotation matrix needed to bring this sensor heading around 
	# to north. 
	Yrot = YawMatrix(hEst)
	# Calculate the rotation matrix to go from sensor to animal frame.
	# Order matters, so the rotation will be done as the inverse of Yrot by
	# inverse of Prot by inverse of Rrot (Yaw then Pitch then Rotation), and 
	# will ultimately be used to pre-multiply the input vectors of the 
	# accelerometer and magnetometer to get all sensor readings into the animal 
	# frame. 
	W = t(Yrot) %*% t(Prot) %*% t(Rrot)	# should be right
	
	W	# return W, a 3x3 rotation matrix
}




################################################################################
################################################################################
# Function		hpr  	(heading, pitch, roll)
# Method based on application note AN4248 Rev 3 2012 
# from Freescale Semiconductor
# Inputs:
#	accVec	n x 3 matrix of acceleration data, in milli-G
#	magVec	n x 3 matrix of magnetometer data, in nanoTesla
#	W		3 x 3 rotation matrix, usually to rotate from sensor frame to 
#			animal frame
# Outputs:
#	results		a n x 3 matrix of heading, pitch, and roll angles, in radians, 
#				in the animal frame (rotated by W) for the supplied data
hpr = function(accVec,magVec, W = NULL) {
	
	# If a rotation matrix W was not supplied, use a simple identity matrix
	# which will not cause any rotation of the input vectors. 
	if (is.null(W)){
		W = matrix(c(1, 	0,  	0,
						0,		1,		0,
						0,		0,		1), nrow = 3, byrow = TRUE)
	}
	
	# Rotate the acceleration vector and magnetometer vector by the rotation
	# matrix W, derived from calibration procedures
	accVecA = matrix(NA, ncol = 3, nrow = nrow(accVec))
	magVecA = matrix(NA, ncol = 3, nrow = nrow(accVec))
	for (i in 1:nrow(accVec)){
		# Go through each row, multiply by W rotation matrix, and put
		# into a matrix. Original version
		accVecA[i,] = t(W %*% accVec[i,])
		magVecA[i,] = t(W %*% magVec[i,])
		# Alternate
#		accVecA[i,] = t(W) %*% accVec[i,]
#		magVecA[i,] = t(W) %*% magVec[i,]
	}
	# Output should be a n x 3 matrix  where the columns are x,y,z axes
	# and the values have been rotated to bring sensor-frame values into the
	# animal frame, as specified by the rotation matrix W
	
	accNorms = numeric(length=nrow(accVecA))
	magNorms = numeric(length = nrow(magVecA))
	for (i in 1:nrow(accVecA)){
		accNorms[i] = norm(accVecA[i,],"2")
		magNorms[i] = norm(magVecA[i,],"2")
	}
	
	accVecA = accVecA / accNorms
	magVecA = magVecA / magNorms

	# Calculate roll angle as the inverse tangent of the animal/phone frame 
	# y axis and z axis (roll is registered as a change on y and z axes)
	# Eq 13
	r = atan2(accVecA[,2],accVecA[,3])	# result in radians
	
	# Calculate pitch angle as the inverse tangent of accelerometer x axis 
	# divided by the accel y + z axes and roll angle (in radians)
	# Eq 15 uses the atan function so that the result is restricted in the range
	# from -90 to +90. 
	p = atan(-accVecA[,1]/ (accVecA[,2]*sin(r) + accVecA[,3]*cos(r)))

	# Restrict pitch to range -90 to 90 degrees
	p[p>(90*180/pi)] = (180*pi/180) - p[p>(90*180/pi)]
	p[p<(-90*180/pi)] = (-180*pi/180) - p[p<(-90*180/pi)]
	
	# Compute rotation matrices for Pitch and Roll using p and r values
	# Apply the RollMatrix and PitchMatrix functions to calculate a set of
	# Roll and Pitch rotation matrices for every entry in magVec
	Rrot = RollMatrix(r)
	Prot = PitchMatrix(p)
	
	# Pre-multiply magVecA by inverse of Rrot and Prot to get magnetometer data 
	# into the horizontal plane
	magVecH = matrix(NA, nrow = nrow(magVecA), ncol = ncol(magVecA),
			byrow = TRUE)
	if (class(Prot) == 'matrix') {
		magVecH = t(Prot) %*% t(Rrot) %*% t(magVecA)
		if (ncol(magVecH) == 1) {
			magVecH = t(magVecH) # convert back to horizontal vector
		}
	} else if (class(Prot) == 'array') {
		for (i in 1:nrow(magVecH)){
			magVecH[i,] = t(Prot[,,i]) %*% t(Rrot[,,i]) %*% magVecA[i,]	
		}	
	}
	# AN4248 Eq 22 solves for yaw (heading). When the magVec data are in the 
	# horizontal plane, the heading angle relative to magnetic north is just the
	# angle between the y and x vectors, which we solve for using the inverse 
	# tangent function. We use the atan2 version to allow values between
	# -pi and pi (-180 to +180)
	h = atan2(-magVecH[,2],magVecH[,1])	# result in radians
	
	results = cbind(h, p, r)
	colnames(results) = c('Heading','Pitch','Roll')
	results	# return results vector
}

plothpr = function(hpr){
	plot(hpr[,3], type = 'n', las = 1, ylim = c(-180,180), 
			ylab = 'Degrees', yaxt = 'n')
	rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4], col = 'grey90')
	axis(side = 2, at = seq(-180,180,by=30), las = 1)
	abline(h = seq(-180,180, by = 30), col = 'grey70')
	# Plot roll
	points(hpr[,3], col = 'black', pch = 20)
	# Plot pitch
	points(hpr[,2], col = 'green', pch = 20)
	# Plot heading
	points(hpr[,1], col = 'red', pch = 20)
	legend("topleft", legend = c("roll","pitch","heading"), 
			col = c('black','green','red'), lty = 1, lwd = 2, bty = 'o',
			bg = 'white')
}

