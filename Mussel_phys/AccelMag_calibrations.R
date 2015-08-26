# AccelMag_calibrations.R
# A script to load accelerometer and magnetometer calibration files, clean the
# data, and calculate offset + scaling factors and a rotation matrix that can 
# be used to properly orient the LSM303D sensor and mussel in space.   
# Author: Luke Miller Aug 24, 2015
###############################################################################

######################################################
# Specify which board and mussel you want to work with here:
mussel = "2" # Specify the mussel, either 1 or 2, as a string
sn = "SN01"	# Specify the board serial number
#####################################################

setwd("D:/R/Mussel_phys/")
#wdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/mag_calibs_20150712/"
#wdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/mag_calibs_20150806/"
#
#fnames = dir(wdir, "*.csv")
#
## Get a vector of mussel numbers from the file names
#musselNum = substr(fnames,4,4)
## Get a vector of serial numbers
#snNames = substr(fnames,23,26)
#
## Find the matching mussel and serial number
#fnum = which( (musselNum == mussel) & (snNames == sn) )
#
## If there is more than 1 matching file, ask the user which they prefer
#if (length(fnum) > 1) {
#	cat("Enter the number of a file from the list\n")
#	for (i in 1:length(fnum)) {
#		cat(i,': ')
#		cat(fnames[fnum[i]],'\n')
#	}
#	choice = scan(what = integer(), nmax = 1)
#	fnum = choice # overwrite fnum
#}
## Get the chosen file name
#datfile = fnames[fnum]

# Interactive file chooser
datfile = file.choose()


# Read the file's metadata in the first row
fileMeta = scan(datfile, what = character(), sep = ',', nlines = 1)
# The 5th entry is the accelerometer setting (usually +/- 4g), and the 7th 
# entry is the magnetometer setting, usually +/- 8 gauss.

# Load the file into memory, skipping the 1st row which has metadata
df = read.csv(datfile, skip = 1)

# Convert acceleration data to milli-g
if (fileMeta[5] == '4 g'){
	df[,2:4] = df[,2:4] * 0.122
}
# Convert magnetometer data to milli-gauss (and then to nanoTesla)
if (fileMeta[7] == '8') {
	df[,5:7] = df[,5:7] * 0.320 * 100
}
# For reference, the total magnetic field norm at HMS is around 47,894 nanoTesla
# magNorm = 47894.9	# nanoTesla local magnetic field total norm

###############################################################################
# Plot the data in 2-D plots
library(plotrix)
par(mfcol=c(3,2))
# X vs Y
plot(df$a.x, df$a.y, type = 'p', 
#		xlim = c(-1,1),
#		ylim = c(-1,1),
		asp = 1,
		las = 1)
points(df$a.x[1:1000], df$a.y[1:1000], col = 'blue')
grid()

# Y vs Z
plot(df$a.y, df$a.z, type = 'p', 
#		xlim = c(-1,1),
#		ylim = c(-1,1),
		asp = 1,
		las = 1)
points(df$a.y[1:1000], df$a.z[1:1000], col = 'blue')
grid()

# X vs Z
plot(df$a.x, df$a.y, type = 'p', 
#		xlim = c(-1,1),
#		ylim = c(-1,1),
		asp = 1,
		las = 1)
points(df$a.x[1:1000], df$a.z[1:1000], col = 'blue')
grid()

# Mag X vs Y
plot(df$m.x, df$m.y, type = 'p', asp = 1, las = 1)
points(df$m.x[1:1000], df$m.y[1:1000], col = 'blue')
grid()
# Mag Y vs Z
plot(df$m.y, df$m.z, type = 'p', asp = 1, las = 1)
points(df$m.y[1:1000], df$m.z[1:1000], col = 'blue')
grid()
# Mag X vs Z
plot(df$m.x, df$m.z, type = 'p', asp = 1, las = 1)
points(df$m.x[1:1000], df$m.z[1:1000], col = 'blue')
grid()

##################################################################
#  A function to identify and overplot chose points, used in functions below
identifyPch <- function(x, y = NULL, n = length(x), pch = 19, ...)
{
	xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
	sel <- rep(FALSE, length(x)); res <- integer(0)
	while(sum(sel) < n) {
		ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
		if(!length(ans)) break
		ans <- which(!sel)[ans]
		points(x[ans], y[ans], pch = pch, col = 'red')
		sel[ans] <- TRUE
		res <- c(res, ans)
	}
	res
}

# A function to plot accelerometer or magnetometer data and use the identifyPch
# function to choose points to be excluded from the dataset.
filterPoints = function(dataF, axisA = 'a.x', axisB = 'a.y'){
	par(mfrow=c(1,1))
	# Plot the input data
	plot(dataF[,axisA],dataF[,axisB], type = 'p', asp = 1, las = 1,
			xlab = axisA,
			ylab = axisB)
	points(dataF[1:1000,axisA],dataF[1:1000,axisB], col = 'blue')
	grid()
	# The user can now use the identifyPch function to click on
	# points that they want to remove
	removeIndex = identifyPch(dataF[,axisA], dataF[,axisB])
	# If at least one point is selected, remove it from the input data set
	if (length(removeIndex) > 0){
		# Recall that the first 10 seconds of calibration should be fixed with
		# the animal in a known position. These points should always be saved
		# Sampling rate was typically 100Hz, so we'll want to save 100Hz * 10s
		# samples at the start.
		removeIndex = removeIndex[removeIndex > 1000]
		# Remove the identified points from the data frame
		filtdataF = dataF[-removeIndex,]
	} else {
		# If user quits out early without choosing any points
		filtdataF = dataF
	}

	# Return the pared-down dataset
	filtdataF
}

# A function to allow the user to choose points to exclude, then loop back and
# view the filtered values, and offer the chance to filter further on the same
# axes or different axes. 
filterLoop = function(dat1){
	repeatFlag = TRUE
	df2 = dat1 # Make a copy of the input data
	# Ask the user if they want to filter the accelerometer data
	cat("Filter data? y or n\n")
	answer = scan(what = character(), n = 1)
	# If they answer no, kill the session immediately
	if (answer == 'n') { repeatFlag = FALSE}
	while (repeatFlag){

		# If they answer yes, start the filtering process
		if (answer == 'y'){
			cat('Choose one axis: ')
			cat(colnames(df2)[-1])
			cat('\n')
			answer1 = scan(what = character(), n = 1)
			cat('Choose 2nd axis: ')
			cat(colnames(df2)[-1])
			cat('\n')
			answer2 = scan(what = character(), n = 1 )
			# Use filterPoints on the x and y accel data
			# This overwrites the original copy in df2
			df2 = filterPoints(df2, axisA = answer1, axisB = answer2)
		}
		# Plot the result
		plot(df2[,answer1], df2[,answer2], type = 'p', asp = 1, las = 1, 
				main = 'Filtered data',
				xlab = answer1,
				ylab = answer2)
		points(df2[1:1000,answer1], df2[1:1000,answer2], col = 'blue')
		grid()
		############################
		cat("Repeat for another round? y or n\n")
		answer = scan(what = character(), n = 1)
		if (answer == 'y'){
			repeatFlag = TRUE
		} else {
			repeatFlag = FALSE # kill the loop
		}
	}  # end of while loop
	df2 # Return the df2 data frame of filtered values
}
###############################################################################
# Utilize the functions above to have the user filter the bad data points
filteredData = filterLoop(df) 
###############################################################################


# Give the user the option to save the output
cat("Save data frame? y on n\n")
answer = scan(what = character(), n = 1)
if (answer == 'y') {
	# Write the header metadata to the file first
	write(fileMeta, file = paste0("Filtered_",datfile), 
			ncolumns = length(fileMeta), sep = ',')
	# Now write the rest of the filtered data
	write.table(filteredData, file = paste0("Filtered_",datfile), 
			row.names = FALSE, append = TRUE, sep = ",")
}








