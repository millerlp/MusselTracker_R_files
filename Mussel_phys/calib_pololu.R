# calib_pololu.R
# 
# Author: Luke Miller Aug 3, 2015
###############################################################################
setwd('D:/R/Mussel_phys/')
require(rgl)
require(RColorBrewer)
require(plotrix)
source("D:/R/Mussel_phys/calib_functions.R")

# Calibration script for the Pololu board with LSM303D on it. 

#fname = "D:/temp/pololu_calib.csv"
fname = 'D:/temp/pololu_calibration_jig_20150805.csv' # Big calibration done on gimbal jig
dfo = read.csv(fname, header = FALSE)

# Extract raw acceleration and magnetometer vectors
accVec = dfo[,c(1,2,3)]
magVec = dfo[,c(4,5,6)]

accRange = 2	# +/- 2g
magRange = 4	# +/- 4gauss



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
# rearrange the accel and magnetometer axes to conform to the frame. 
accVec[,1] = -1 * accVec[,1] # invert accelerometer x axis
magVec[,2] = -1 * magVec[,2] # invert magnetometer y-axis so east is positive
magVec[,3] = -1 * magVec[,3] # invert magnetometer z-axis so down is positive

accNorm = 1000 # 1000 milli-g local acceleration (gravity) norm
magNorm = 47894.9	# nanoTesla local magnetic field total norm

##########################################
# Use calibrate function to produce scaling and offset factors 
accCal = calibrate(accVec[,1],accVec[,2],accVec[,3], accNorm)
magCal = calibrate(magVec[,1],magVec[,2],magVec[,3], magNorm)

# Extract the calibrated data vectors. Units should still be in milli-G and
# nanoTesla
dfcal = data.frame(accCal$Calibrated)
dfcal = cbind(dfcal,data.frame(magCal$Calibrated))
names(dfcal) = c('accX','accY','accZ','magX','magY','magZ')


# Also grab the original uncalibrated data
df = accVec
df = cbind(df, magVec)
names(df) = c('accX','accY','accZ','magX','magY','magZ')

df2 = df
for (i in 1:nrow(df)){
	# Subtract bias value for each axis
	df2[i,1] = df[i,1] - accCal$bias[1]
	df2[i,2] = df[i,2] - accCal$bias[2]
	df2[i,3] = df[i,3] - accCal$bias[3]
	softiron = accCal$softiron
	# Pre-multiply the bias-corrected values by the softiron matrix
	df2[i,1:3] =  softiron %*% t(as.matrix(df2[i,1:3], nrow = 1, ncol = 3))
	# Subtract bias value for each axis
	df2[i,4] = df[i,4] - magCal$bias[1]
	df2[i,5] = df[i,5] - magCal$bias[2]
	df2[i,6] = df[i,6] - magCal$bias[3]
	softiron = magCal$softiron
	# Pre-multiply the bias-corrected values by the softiron matrix
	df2[i,4:6] =  softiron %*% t(as.matrix(df2[i,4:6], nrow = 1, ncol = 3))
}



######################################################
# Use rgl package to make 3D plot

plotAcc = TRUE
plotBoth = TRUE

# Input data should be in milli-G and nanoTesla
plot3Dcalib = function(df, dfcal, plotAcc = TRUE, plotBoth = TRUE){
	open3d(windowRect = c(30,30,630,630))
	pal = brewer.pal(3, 'Pastel2' )
	 
	# If plotting accel data, assume X,Y,Z are in columns 1, 2, 3
	if (plotAcc) {
		colX = 1
		legendtitle = 'Acceleration'
	} else {
		# If plotting magnetometer data, assume X,Y,Z are in columns 4,5,6
		colX = 4
		legendtitle = 'Magnetic field'
	}
	# Determine range limits
	if (plotBoth) {
		xlims = range(c(df[,c(colX:(colX+2))], dfcal[,c(colX:(colX+2))]))
		zlims = ylims = xlims
	} else if (!plotBoth) {
		xlims = range(dfcal[,c(colX:(colX+2))])		
		zlims = ylims = xlims
	}
	
# When using Pastel2 palette, Green will be the unedited points, Orange will
# be the calibrated points.
	if (plotBoth) {
		# Plot both uncalibrated and calibrated data
		plot3d(df[,colX],df[,(colX+1)],df[,(colX+2)], 
				xlab = 'X', ylab='Y', zlab='Z', col = pal[1], 
				aspect = c(1,1,1), 
				xlim = xlims, ylim = ylims, zlim = zlims)
		plot3d(dfcal[,colX],dfcal[,(colX+1)],dfcal[,(colX+2)], 
				add = TRUE, col = pal[2], 
				aspect = c(1,1,1))
	} else if (!plotBoth) {
		# Only plot the calibrated data
		plot3d(dfcal[,colX],dfcal[,(colX+1)],dfcal[,(colX+2)], 
				xlab = 'X', ylab='Y', zlab='Z', col = pal[2], 
				aspect = c(1,1,1), 
				xlim = xlims, ylim = ylims, zlim = zlims)
	}
	grid3d(side = c("x","y","z")) # plot grid lines
	# get bounding box coordinates (xmin, xmax, ymin, ymax, zmin, zmax)
	lims = par3d('bbox')

	if (plotBoth) {
		# Plot X vs Y on lower Z plane
		plot3d(df[,colX],df[,(colX+1)], z = rep(lims[5],nrow(df)), add = TRUE,
				col = pal[1], aspect = c(1,1,1))
		plot3d(dfcal[,colX],dfcal[,(colX+1)], z = rep(lims[5],nrow(dfcal)), 
				add = TRUE,
				col = pal[2], aspect = c(1,1,1))
	# Plot X vs Z on lower Y plane
		plot3d(df[,colX],y = rep(lims[3],nrow(df)),df[,(colX+2)],
				add = TRUE,
				col = pal[1], aspect = c(1,1,1))
		plot3d(dfcal[,colX],y = rep(lims[3],nrow(dfcal)),dfcal[,(colX+2)],
				add = TRUE,
				col = pal[2], aspect = c(1,1,1))
	# Plot Y vs Z on lower X plane
		plot3d(x = rep(lims[1],nrow(dfcal)),df[,(colX+1)],df[,(colX+2)], 
				add = TRUE,
				col = pal[1], aspect = c(1,1,1))
		plot3d(x = rep(lims[1],nrow(dfcal)),dfcal[,(colX+1)],dfcal[,(colX+2)],
				add = TRUE,
				col = pal[2], aspect = c(1,1,1))
		legend3d('topright',c('Raw','Corrected'), col = c(pal[1],pal[2]), 
				pch = 20, title = legendtitle,
				cex = 2)
	} else if (!plotBoth) {
		# Plot X vs Y on lower Z plane
		plot3d(dfcal[,colX],dfcal[,(colX+1)], z = rep(lims[5],nrow(dfcal)), 
				add = TRUE,
				col = pal[2], aspect = c(1,1,1))
	# Plot X vs Z on lower Y plane
		plot3d(dfcal[,colX],y = rep(lims[3],nrow(dfcal)),dfcal[,(colX+2)], 
				add = TRUE,
				col = pal[2], aspect = c(1,1,1))
	# Plot Y vs Z on lower X plane
		plot3d(x = rep(lims[1],nrow(dfcal)),dfcal[,(colX+1)],dfcal[,(colX+2)],
				add = TRUE,
				col = pal[2], aspect = c(1,1,1))
		legend3d('topright',c('Corrected'), col = c(pal[2]), pch = 20,
				title = legendtitle,
				cex = 2)
	}
} # end of plot3Dcalib function

plot3Dcalib(df,dfcal,plotAcc = TRUE, plotBoth = TRUE)
plot3Dcalib(df,dfcal,plotAcc = FALSE, plotBoth = TRUE)

###############################################################################
# Plot the data in 2-D plots

par(mfrow=c(2,2))
# X vs Y
plot(df$X, df$Y, type = 'p', 
		xlim = c(-1,1),
		ylim = c(-1,1),
		asp = 1,
		las = 1)
grid()

#points(dfoff$X, dfoff$Y, col = 'blue') # Plot offset-corrected values
points(dfcal$X, dfcal$Y, col = 'red', pch = 20) # Plot bias-corrected values
draw.circle(0,0, radius = 0.569, border = 'blue', lwd = 2)
# Y vs Z
plot(df$Y, df$Z, type = 'p', 
		xlim = c(-1,1),
		ylim = c(-1,1),
		asp = 1,
		las = 1)
grid()
#points(dfoff$Y, dfoff$Z, col = 'blue') # Plot offset-corrected values
points(dfcal$Y, dfcal$Z, col = 'red', pch = 20) # Plot bias-corrected values
draw.circle(0,0, radius = 0.569, border = 'blue', lwd = 2)
# X vs Z
plot(df$X, df$Z, type = 'p', 
		xlim = c(-1,1),
		ylim = c(-1,1),
		asp = 1,
		las = 1)
grid()
#points(dfoff$X, dfoff$Z, col = 'blue') # Plot offset-corrected values
points(dfcal$X, dfcal$Z, col = 'red', pch = 20) # Plot bias-corrected values
draw.circle(0,0, radius = 0.569, border = 'blue', lwd = 2)
plot.new() # for the 4th empty plot