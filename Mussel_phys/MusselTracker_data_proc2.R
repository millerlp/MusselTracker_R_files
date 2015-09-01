# MusselTracker_data_proc.R
# Procedures to concatenate data files from a MusselTracker board into a single
# continuous time series, including filtering of data, removing bad data, 
# making small interpolations for small gaps (less than 30 seconds) etc. 
# Author: Luke Miller Jul 9, 2015
###############################################################################

# Enter the serial number of the board files you want to process below:
##########################################
board = 'SN03'
##########################################


require(RColorBrewer)
require(signal) # For low-pass filtering
require(zoo) # For merging time series

Sys.setenv(TZ = 'PST8PDT')	# Set the time zone for the mussel tracker data
Sys.getenv('TZ') # Report back to show that time zone setting worked


## Show fractional seconds when printing to terminal. This still often doesn't
## properly show the fractional seconds when displaying the formatted time
options(digit.secs = 6, digits = 12) 


# This assumes that there will be a folder of the same name (SN14 for example)
# in which the script will find all of the files. Specify the path to the 
# directory below.

# Specify the directory holding the field data files
fdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Field_data/"
# Specify the directory holding the thermocouple calibration data file
tcdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Thermocouple_calibrations/20150630_calibrations/"
# Specify the file with the timestamps that need to be excised
exfile = paste0(fdir,'Timestamps_excise_list.csv')


# Specify a start and end time cutoff for the dataset. The plates were 
# deployed in the field by 2015-07-15 7:30PDT, prior points are in the water table
cutoffstart = as.POSIXct('2015-07-15 07:30:00', tz = 'PST8PDT')
cutoffend = as.POSIXct('2015-08-06 09:20', tz = 'PST8PDT')

fname1 = paste0(fdir,board)
fnames = dir(fname1, pattern = '*.csv', full.names = TRUE)

# Open all the field data files and drop them in one large data frame
tot = length(fnames)
cat('Opening ',tot,'files\n')
for (i in 1:length(fnames)){
	# Read the first line of the file, which contains board serial number and
	# other info
	info = scan(fnames[i], nlines = 1, what = 'character', sep = ',')
	sn = info[1] # Extract the serial number
	# Only proceed if the serial number in the header file matches the 
	# value the user specified in 'board' variable. 
	if (sn == board) {
		dftemp = read.csv(fnames[i], skip = 1,
				colClasses = c('integer','POSIXct',rep('integer',7),
						'numeric',rep('integer',7),'numeric','integer'),
				na.strings = c('NA','nan','NaN'))
		cat('Processed file ', i,'\n'); flush.console()
		if (i == 1) {
			df = dftemp
		} else {
			df = rbind(df,dftemp)
		}
	} else {
		warning('Serial number mismatch: you asked for ', board,'\n')
		warning(basename(fnames[i]),' returned ', sn, '\n')
	}

}
cat('Finished importing files\n')

# Add the fractional seconds value onto the unix timestamps
df$POSIXt = df$POSIXt + (df$fractional.Second / 100)

df$DateTime = df$DateTime + (df$fractional.Second/100)
# Note that although the printed time stamps will be in the local time zone,
# the internal representation of the POSIXct value (a numeric value of seconds
# since the start of the epoch (origin) will be in the UTC time zone, which can
# be observed by doing as.numeric(df$DateTime[1]), and plugging the result into
# a website like http://www.unixtimestamp.com/index.php to show how that value
# is converted into the UTC time zone. Or simply compare the value in 
# df$POSIXt[1] with the result of as.numeric(df$DateTime[1]). 

# Reorder the data frame based on the POSIXt column, since we've pasted several
# input data files together
df = df[order(df$POSIXt),]


# Add a column to hold the board serial number, insert it near the front of
# the data frame
df = cbind(df[,1:2], data.frame(SerialNumber=rep(sn,nrow(df))), df[3:ncol(df)])

dfo = df # keep a copy of the unmolested input data frame in dfo
###############################################################################
# Filter suspect time points. Suspect time points can have ms values that aren't
# equal to 0, 25, 50, or 75. 
# Make a matrix to hold the test results
filt = matrix(0,nrow = nrow(df), ncol = 5)
filt[,1] = rep(0L,nrow(df)) == df[,'fractional.Second'] # compare every ms value to 0
filt[,2] = rep(25L,nrow(df)) == df[,'fractional.Second'] # compare every ms value to 25
filt[,3] = rep(50L,nrow(df)) == df[,'fractional.Second'] # compare every ms value to 50
filt[,4] = rep(75L,nrow(df)) == df[,'fractional.Second'] # compare every ms value to 75
# Sum every row of filt. If none of the tests above returned true, the result
# in the 5th column will be 0 (== FALSE). Otherwise it will be TRUE. The goal
# is to have one of the 4 values match so you get a TRUE in column 5.
filt[,5] = rowSums(filt[,1:4]) 
# Now find every row that has a suspect ms value. A 0 value in filt[,5] 
# indicates that none of the fractional seconds values matched. 
badrows = which(!(filt[,5]))
if (length(badrows)>0){
	# Remove any rows with suspect ms values
	df = df[-badrows,]	
	badrowtot = length(badrows)
} else {
	badrowtot = 0
}

# Look for spurious future dates. Take the last timestamp in the dataset, which
# we'll assume is real, and look for any timestamps in the
# body of the dataset larger than it. This works on the DateTime column, 
# but in most cases the POSIXt column values will be spurious small values, so
# this test wouldn't work on the POSIXt column
badrows = which(df$DateTime > df$DateTime[nrow(df)])
if (length(badrows)> 0) {
	df = df[-badrows,]
	badrowtot = badrowtot + length(badrows)
}
# Look for any spurious early dates. We haven't seen any of these yet in the
# data, but it's worth checking anyways. 
badrows = which(df$DateTime < df$DateTime[1])
if (length(badrows) > 0) {
	df = df[-badrows,]
	badrowtot = badrowtot + length(badrows)
}
# Some of the spurious rows have a sensible timestamp and fractional.Seconds
# value, so the only real indication that they're bad is that the timestamp
# exactly matches the timestamp in the previous second's worth of data. 
# Look for any duplicated POSIXt timestamps, the duplicated() function will
# return the row index of the 2nd entry, which in our datasets will be the 
# spurious row
badrows = which(duplicated(df$POSIXt))
if (length(badrows) > 0){
	df = df[-badrows,]
	badrowtot = badrowtot + length(badrows)
}

cat('Removed',badrowtot,'suspect rows from data.\n')
# That should take care of any spurious values
#******************************************************************************

###############################################################################
# The data frame df should already be ordered. We can create a larger data
# frame of all time points from start to finish, including gaps in the files,
# by using the first POSIXt value and the last, with a 0.25sec interval.
cat('Expanding dataset to entire time period\n')
st = df$POSIXt[1]
et = df$POSIXt[nrow(df)]
# Make the entire sequence of POSIXt time stamps into a data frame
dftimes = data.frame(POSIXt = seq(st,et, by = 0.25))
# Merge the two dataframes, and missing values will be filled with NA's
df2 = merge(df,dftimes, by = 'POSIXt', all = TRUE)
# Since some DateTime values will be missing in the expanded dataframe, fill
# them in by converting the POSIXt values into POSIXct values in DateTime column
# Start by shifting all time values forward 7 hours, since the conversion to 
# a POSIXct value will unfailingly assume the numeric value is in GMT, so 
# even if you specify a time zone, it will be shifting that your times to fit
# that time zone. 
df2$POSIXt = df2$POSIXt + (60 * 60 * 7) # shift time values forward 7 hrs
df2$DateTime = as.POSIXct(df2$POSIXt, origin = '1970-1-1', tz = "PST8PDT")
# Fill in the board serial number column to fill any missing rows
df2$SerialNumber = board

# Remove data from before and after the field deployment, if present.
df2 = df2[df2$DateTime >= cutoffstart,]
df2 = df2[df2$DateTime <= cutoffend,]

###############################################################################
# Recall that the temperature and hall effect data were only sampled 1 time
# per second, with the same value being written to all 4 sub-second samples.
# As a result, we ought to pull out the "true" 1-second values, resulting in a 
# much smaller data set, and also leave the accelerometer/magnetomter values
# out, since they were sampled at 4Hz. 

df3 = df2[seq(1,(nrow(df2)-3), by = 4),c('POSIXt','DateTime','SerialNumber',
				'Temp1','Hall1', 'Temp2','Hall2')]
df3$SerialNumber = as.factor(df3$SerialNumber)

# Remove the columns of df2 pertaining to Temp and Hall data, retaining the
# accel and magnetometer data
df2 = df2[c(1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18)]
df2$SerialNumber = as.factor(df2$SerialNumber)



################################################################################
# The time periods with questionable data are listed in the file
# Timestamps_excise_list.csv. Open it and remove any data from df2 and df3 that 
# are suspect. 

excisedf = read.csv(exfile)
excisedf$StartIgnore = as.POSIXct(excisedf$StartIgnore, 
		format = '%m/%d/%Y %H:%M')
excisedf$EndIgnore = as.POSIXct(excisedf$EndIgnore, 
		format = '%m/%d/%Y %H:%M')

##############################################################
# exciseTempHall: A function to excise questionable rows of data from the 1-Hz 
# temperature and Hall effect sensor data frame
# Inputs:
#	temphall: a data frame of temperature and hall effect data from the
#			Mussel Tracker raw data files
#	excisedf: a data frame of start and end times to ignore in the data set
# Output:
# 	temphall: a modified version of the input data frame with the questionable
#			rows replaced with NAs. 

exciseTempHall = function(temphall, excisedf) {
	# Find the rows in excisedf that have the same serial number as the current
	# board, and have time points to remove for each channel and sensor combo
	# Start with Ch1 Temperature
	matches = which(excisedf$SerialNumber == levels(temphall$SerialNumber)[1] &
					excisedf$Mussel == 'Ch1' & excisedf$Sensor == 'Temp')
	if (length(matches) > 0){
		for (i in 1:length(matches)){
			startrow = which.min(abs(temphall$DateTime - 
									excisedf$StartIgnore[matches[i]]))
			endrow = which.min(abs(temphall$DateTime - 
									excisedf$EndIgnore[matches[i]]))
			# blank out the suspect rows of data
			temphall$Temp1[startrow:endrow] = NA
		}
	}
	# Do the same for Channel 2 temperature
	matches = which(excisedf$SerialNumber == levels(temphall$SerialNumber)[1] &
					excisedf$Mussel == 'Ch2' & excisedf$Sensor == 'Temp')
	if (length(matches) > 0){
		for (i in 1:length(matches)){
			startrow = which.min(abs(temphall$DateTime - 
									excisedf$StartIgnore[matches[i]]))
			endrow = which.min(abs(temphall$DateTime - 
									excisedf$EndIgnore[matches[i]]))
			# blank out the suspect rows of data
			temphall$Temp2[startrow:endrow] = NA
		}
	}
	# Next remove bad data from Hall sensor Ch1
	matches = which(excisedf$SerialNumber ==levels(temphall$SerialNumber)[1] &
					excisedf$Mussel == 'Ch1' & excisedf$Sensor == 'Hall')
	if (length(matches) > 0){
		for (i in 1:length(matches)){
			startrow = which.min(abs(temphall$DateTime - 
									excisedf$StartIgnore[matches[i]]))
			endrow = which.min(abs(temphall$DateTime - 
									excisedf$EndIgnore[matches[i]]))
			# blank out the suspect rows of data
			temphall$Hall1[startrow:endrow] = NA
		}
	}
	# And remove bad data from Hall sensor Ch2
	matches = which(excisedf$SerialNumber == levels(temphall$SerialNumber)[1] &
					excisedf$Mussel == 'Ch2' & excisedf$Sensor == 'Hall')
	if (length(matches) > 0){
		for (i in 1:length(matches)){
			startrow = which.min(abs(temphall$DateTime - 
									excisedf$StartIgnore[matches[i]]))
			endrow = which.min(abs(temphall$DateTime - 
									excisedf$EndIgnore[matches[i]]))
			# blank out the suspect rows of data
			temphall$Hall2[startrow:endrow] = NA
		}
	}
	temphall # return data frame as output
} # End of exciseTempHall function

###############################################################################
# Do the same removal procedure for the suspect magnetometer and accelerometer
# data
# exciseAccelMag: a function to excise questionable rows of data from 4Hz 
# accelerometer and magnetometer data frame
# Inputs:
# 	accelmag: a data frame of accel and magnetometer data from the 
# 			Mussel Tracker raw data files.
#	excisedf: a data frame of start and end times for data that should be 
# 			ignored
# Output:
# 	accelmag: a modified version of the input data frame with the spurious data
# 				rows replaced with NAs.
exciseAccelMag = function(accelmag, excisedf) {
	# Start with Ch1 accel/mag
	matches = which(excisedf$SerialNumber == levels(accelmag$SerialNumber)[1] &
					excisedf$Mussel == 'Ch1' & excisedf$Sensor == 'accelmag')
	if (length(matches) > 0){
		for (i in 1:length(matches)){
			startrow = which.min(abs(accelmag$DateTime - 
									excisedf$StartIgnore[matches[i]]))
			endrow = which.min(abs(accelmag$DateTime - 
									excisedf$EndIgnore[matches[i]]))
			# blank out the suspect rows of data
			accelmag$a1.x[startrow:endrow] = NA
			accelmag$a1.y[startrow:endrow] = NA
			accelmag$a1.z[startrow:endrow] = NA
			accelmag$m1.x[startrow:endrow] = NA
			accelmag$m1.y[startrow:endrow] = NA
			accelmag$m1.z[startrow:endrow] = NA
		}
	}
	# Then do Ch2 accel/mag 
	matches = which(excisedf$SerialNumber == levels(accelmag$SerialNumber)[1] &
					excisedf$Mussel == 'Ch2' & excisedf$Sensor == 'accelmag')
	if (length(matches) > 0){
		for (i in 1:length(matches)){
			startrow = which.min(abs(accelmag$DateTime - 
									excisedf$StartIgnore[matches[i]]))
			endrow = which.min(abs(accelmag$DateTime - 
									excisedf$EndIgnore[matches[i]]))
			# blank out the suspect rows of data
			accelmag$a2.x[startrow:endrow] = NA
			accelmag$a2.y[startrow:endrow] = NA
			accelmag$a2.z[startrow:endrow] = NA
			accelmag$m2.x[startrow:endrow] = NA
			accelmag$m2.y[startrow:endrow] = NA
			accelmag$m2.z[startrow:endrow] = NA
		}
	}
	accelmag	# Return the modified accelmag data frame
}
###############################################################################
# Apply the excise functions to replace spurious data with NA's. 
cat('Removing suspect timepoints, replacing with NAs\n')
df2 = exciseAccelMag(df2, excisedf)
df3 = exciseTempHall(df3, excisedf)

#######*********************************************************************
# There will be small gaps in the dataset of a few seconds, along with some
# longer gaps. Use the functions in the zoo package to fill in the brief 
# missing data with the prior values, defining a maxgap beyond which values
# won't be filled in. 
cat('Filling small gaps in the dataset\n')
# Convert the data frame to a zoo object, leaving behind the POSIXt, DateTime
# and SerialNumber columns so that only numeric data are included
z = zoo(df3[, -c(1,2,3)], order.by = df3$POSIXt)
# For brief missing data, carry forward the last non-NA value, up to a maximum
# of maxgap rows. For these 1-second data, 30 = thirty second gaps max.
z2 = na.locf(z, maxgap = 30)
# Copy the filled data back into df3
df3[,4:7] = z2
# At this point there should only be significant gaps of missing data during
# times when the datalogger was not recording for some reason, due to battery
# failure, interruptions for downloading, etc. 

###############################################################################
cat('Applying temperature calibrations\n')
#*********************************
# Find any spurious 0C temperature values and convert them to NAs
df3$Temp1[df3$Temp1 < 0.25] = NA
df3$Temp2[df3$Temp2 < 0.25] = NA
df3$Temp1[df3$Temp1 > 50] = NA
df3$Temp2[df3$Temp2 > 50] = NA

# Apply the calibration data to the temperature records produced by
# Thermocouple_calibration.R
tcfile = 'Calibration_coefficients_20150630.csv'

tc = read.csv(paste0(tcdir,tcfile))

# Find the matching intercept and slope coefficients in the calibration file
# so that we can convert each thermocouple temperature into a true 
# corrected temperature. 

# Compare the value of 'board' defined above with the available values in the
# Serial column of the tc data frame. 
# Extract the intercept and slope values for the thermocouple calibration. 
tc1int = tc[as.character(tc$Serial) == board & tc$Channel == 1, 'Intercept']
tc1sl =  tc[as.character(tc$Serial) == board & tc$Channel == 1, 'Slope']
tc2int = tc[as.character(tc$Serial) == board & tc$Channel == 2, 'Intercept']
tc2sl =  tc[as.character(tc$Serial) == board & tc$Channel == 2, 'Slope']

applyTCcalib = function(temps, calibInt, calibSlope){
	# Apply the calibration slope and intercept to all input temperatures
	temps2 = (temps * calibSlope) + calibInt
	# The "calibrated" temperatures will now have greater precision than the
	# original sensor was capable of. We need to round every value to the 
	# nearest 0.25°C. This can be done by multiplying the value by 4, rounding
	# to the nearest whole number, and then dividing that whole number by 4.
	# For example: 25.24 * 4 = 100.95, rounds to 101, 101/4 = 25.25
	temps2 = round(temps2 * 4) / 4
	temps2 	# return the rounded, calibrated values
}

df3$Temp1calib = applyTCcalib(df3$Temp1, tc1int, tc1sl)
df3$Temp2calib = applyTCcalib(df3$Temp2, tc2int, tc2sl)



##*****************************************************************************
# Use the hall effect data to calculate a % of max opening based on the
# maximum and minimum values of the hall effect data for a mussel

# A function to determine where there are gaps of missing data, so that we can 
# apply a filter to the good chunks of data. Most filters fail if there are NAs
# present, and our remaining gaps in the data are big enough that they probably
# shouldn't be smoothed-over anyhow. 

gapBounds = function (values){
	# This function returns a data frame of 2 columns, start + end, that give
	# the row indices of the start and end of each run of good data in the 
	# input vector. It should handle data streams that start with NAs and
	# should handle streams that end with NAs. 
	
	# Use run length encoding function to find NA gaps
	gaps = rle(is.na(values))
	# If the first entry in gaps is TRUE, this indicates that the data started
	# with NA values.
	if (gaps$values[1] == TRUE) {
		startNA = TRUE
	} else {	
		startNA = FALSE	 # Data started with real values
	}
	
	# The number of gaps with value = TRUE is the number of good runs of data
	# and gaps with value == FALSE are the starts of runs of NAs. This will
	# get the number of FALSE (negated to TRUE) values in gaps. A dataset
	# ending in good data will typically have a last entry in gaps of FALSE 
	numgaps = sum(!gaps$values)
	# Create the output data frame
	results = data.frame(Start = integer(numgaps), End = integer(numgaps))

	# The first entry should always be 1 (1st row) for data that start with 
	# real values. If there are no gaps, the
	# contents of gaps$lengths will just be the length of the values vector
	results$Start[1] = 1
	results$End[1] = gaps$lengths[1]
	# However, if the dataset starts with NAs, the first entry in gaps will
	# be the index of the first good data, while the 2nd entry will be the 
	# start of the next stretch of NAs
	if (startNA) {
		results$Start[1] = gaps$lengths[1]+1
		results$End[1] = sum(gaps$lengths[1:2])
	}
	
	# If there is more than 1 entry in gaps$lengths, process the other gaps
	j = 2; # counter
	if (numgaps > 1){
		if (numgaps %% 2 == 0){
			if (!startNA){		
				# Even number of gaps, dataset ends on NAs
				for (i in seq(2, length(gaps$lengths)-1, by = 2)){
					nextstart = sum(gaps$lengths[1:i]) + 1
					nextend = sum(gaps$lengths[1:(i+1)])
					results$Start[j] = nextstart
					results$End[j] = nextend
					j = j + 1
				}
			} else if (startNA){
				# Even number of gaps, dataset started on NAs. The 1st entry
				# is taken care of above, so we start with the 3rd entry which
				# should be the next stretch of good data (i.e. it should be
				# a TRUE entry in gaps$values
				for (i in seq(3, length(gaps$lengths)-1, by = 2)){
					nextstart = sum(gaps$lengths[1:i]) + 1
					nextend = sum(gaps$lengths[1:(i+1)])
					results$Start[j] = nextstart
					results$End[j] = nextend
					j = j + 1
				}
			}
		} else if (numgaps %% 2 != 0){
			# Odd number of gaps
			if (!startNA){
				for (i in seq(2,length(gaps$lengths), by = 2)){
					# The first value in gaps$lengths should represent the end 
					# of a
					# a good run of data, and the 2nd value should represent the
					# length of a run of NAs. As such, the start of the next run
					# of good data is the sum of all the previous gaps$lengths 
					# values so far
					nextstart = sum(gaps$lengths[1:i]) + 1
					nextend = sum(gaps$lengths[1:(i+1)])
					results$Start[j] = nextstart
					results$End[j] = nextend
					j = j + 1
				}
			} else if (startNA){
				# Handle the case where the data stream started with NA values
				for (i in seq(3,length(gaps$lengths), by = 2)){
					# When starting with NAs, the first value in gaps$lengths 
					# should represent the start of a good run of data (TRUE). 
					# As such, the next run of NAs starts with the 3rd entry
					nextstart = sum(gaps$lengths[1:i]) + 1
					nextend = sum(gaps$lengths[1:(i+1)])
					results$Start[j] = nextstart
					results$End[j] = nextend
					j = j + 1
				}
			}
		}	# end of else if statement
	} # end of if (numgaps > 1)
	
	results	# return the results dataframe, 2 columns Start and End
}

###########################################
### hallFilter function
hallFilter = function(hallData){
	# Input a vector of hall effect data, which can contain NAs. This function
	# will first figure out where the good runs of data are, then apply a 
	# Butterworth order 1 low-pass filter to those runs, and output a 
	# vector of filtered data with the NA gaps still intact.
	
	# Define a butterworth low-pass filter, 1st order, to filter at 1/10 the 
	# sampling rate (which was 1Hz)
	myfilter = butter(1,0.1,type='low', plane = 'z')
	# Find the start and end of any gaps in the data set using the gapBounds
	# function defined earlier
	mygaps = gapBounds(hallData)
	# Apply the filter to the good runs of data
	for (i in 1:nrow(mygaps)){
		# Extract the run of good data
		dats = hallData[mygaps$Start[i]:mygaps$End[i]]
		# Find the offset of the data from zero
		offset = dats[1]
		# Subtract the offset off of all of the data in this run
		# (The butterworth filter returns large transient values at the start
		# of the run if the value is much different from zero)
		dats = dats - offset
		# Call the filter routine to apply the filter
		yfiltered = signal:::filter(myfilter,dats)
		# Add the offset back on so the data are back on their original scale
		yfiltered = yfiltered + offset
		# Write the filtered data back into the data vector
		hallData[mygaps$Start[i]:mygaps$End[i]] = yfiltered
	}
	# Round the filtered data back to the nearest whole number, since these
	# represent ADC count data.
	hallData = round(hallData)
	hallData	# return vector of filtered data
}


################################################
### percentileRange function
# Define a function to calculate the average reading for the upper and lower
# x percentile of the hall effect data.  
percentileRange = function (Hallvec, percentileLim = 0.01){
	# Remove any NA's
	temp = Hallvec[!is.na(Hallvec)]
	# Reorder the values smallest to largest
	temp = temp[order(temp)]
	# Get the index of the entry closest to percentileLim
	indx = round(percentileLim * length(temp))
	# Calculate the mean value for the lower 1% of closed valve values
	closedVal = mean(temp[1:indx])
	
	# Now do the same for the other end of the range of hall effect values
	# These would normally represent "fully open" readings near 512 if the 
	# magnet and sensor are situated so that a nearby magnet drives the signal
	# below 512. 
	indx = round((1 - percentileLim) * length(temp))
	# Calculate the mean value for the upper 1% of open valve values
	openVal = mean(temp[indx:length(temp)])
	result = c(closedVal, openVal)
	result # return the two values, always smallest then largest
}

############################################
### calcPercentGape function
# Use this function to calculate the percentage gape (0-100%)
# This calculates the baseline fully-closed and fully-open values based
# on the upper and lower 1% of hall effect sensor data for the given range
# of rows in the input dataset. 
# Supply a vector of Hall effect sensor data. This function will find any 
# gaps in the vector (NA's) and calculate a new baseline for each section of
# good data. 
calcPercentGape = function (hallData){
	# Get the row indices for the major gaps that exist now
	mygaps = gapBounds(hallData)
	# Create an empty output vector of the same length as the input data. 
	outputData = vector(mode = 'numeric', length = length(hallData))
	outputData[] = NA # Convert all the values to NA to begin with
	
	# Now for each section of the input data, calculate the percent gape. This
	# involves using the entries in mygaps to subset the input data and 
	# calculate individual percentile values for each contiguous section of
	# data. If there is a significant gap in the data (usually > 30 sec), the
	# percentages will be re-calculated for the each good section of data. This
	# should accomodate any sensor reposition issues. 
	for (i in 1:nrow(mygaps)){
		st = mygaps$Start[i]
		end = mygaps$End[i]
		# First calculate the upper and lower 1% values of the raw Hall readings
		myrange = percentileRange(hallData[st:end], percentileLim = 0.01)
		# I think filtering by the 1% values is no longer necessary with the 
		# new butteworth filter applied earlier
#		myrange = range(hallData[st:end])
		if (myrange[2] < 512){
			# If the hall effect and magnet are situated so that the value goes 
			# down as the mussel closes, then calculate % opening using the 
			# maximum value. If the maximum value is less than 512, we'll assume
			# that closing drives the value down.
			# Calculate the width of the range of values. If the readings are 
			# less than 512, the smaller value in meanVals represents the 
			# fully-closed reading while the larger value is the fully open 
			# reading. 
			rangeval = myrange[2] - myrange[1]
			# Calculate a percentage opening based on the difference between 
			# each reading and the minimum value (fully closed), divided by the 
			# range
			outputData[st:end] = ((hallData[st:end] - myrange[1]) / 
						rangeval) * 100
		} else if (myrange[2] > 512){
			# Calculate the width of the range of values. If the readings are 
			# larger than 512, the smaller value in myrange represents the 
			# fully-open reading while the larger value is the fully-closed 
			# reading. 
			rangeval = myrange[2] - myrange[1]
			# Calculate a percentage opening based on the difference between 
			# each reading and the maximum value (fully closed), divided by the 
			# range
			outputData[st:end] = ((myrange[2] - hallData[st:end]) / 
						rangeval) * 100
		}
	}	
	outputData	# return the output data
}
 
cat('Filtering Hall effect data\n')
df3$Hall1filt = hallFilter(df3$Hall1)
df3$Hall2filt = hallFilter(df3$Hall2)
cat('Calculating percent gape\n')
df3$Hall1.percent = calcPercentGape(df3$Hall1filt)
df3$Hall2.percent = calcPercentGape(df3$Hall2filt)

# Update column names in df3 to be more descriptive
names(df3) = c('POSIXt','DateTimePDT','SerialNumber','Temp1raw','Hall1raw',
		'Temp2raw','Hall2raw','Temp1calib','Temp2calib','Hall1filt',
		'Hall2filt','Hall1.percent','Hall2.percent')
# Also make DateTime column name more descriptive in df2
names(df2)[2] = 'DateTimePDT'

cat('\aFinished processing input files\n')


##***************************************************************************
# Plot the data streams for mussel 1
plot2 = TRUE # set to false to only show one mussel's data

#########################################################################
# Function to plot one or two channels of temperature data from the board
plotTemps = function(df3, start = 1, end = nrow(df3), plot2 = TRUE){
	if (class(start)[1] == 'character'){
		# Handle character timestamps
		start = as.POSIXct(start, origin = '1970-1-1')
		end = as.POSIXct(end, origin = '1970-1-1')
	}
	if (class(start)[1] == 'POSIXct' & class(end)[1] == 'POSIXct'){
		# Find nearest row for each time stamp
		startrow = which.min(abs(start - df3$DateTimePDT))
		endrow = which.min(abs(end - df3$DateTimePDT))
	}
	if (class(start)[1] == 'numeric'){
		startrow = start
		endrow = end
	}
	# Now assign the result of those if statements to two variables sr and er
	# start row, end row
	sr = startrow
	er = endrow

	cols = brewer.pal(3,'Set1')
# Plot temperature first
	if (plot2){
		ylims = range(c(df3$Temp1calib[sr:er],df3$Temp2calib[sr:er]), 
				na.rm = TRUE)
		ylims[2] = ylims[2] + 2 # bump the upper temp limit up slightly
		plot(df3$DateTimePDT[sr:er], df3$Temp1calib[sr:er], type = 'n', 
				ylab = expression(Temperature*','~degree*C),
				xlab = 'Time',
				col = cols[1],
				ylim = ylims,
				las = 1)
		grid()
		lines(df3$DateTimePDT[sr:er], df3$Temp1calib[sr:er], col = cols[1])
		lines(df3$DateTimePDT[sr:er], df3$Temp2calib[sr:er], col = cols[2])	
		legend('topleft', legend = c('Mussel 1','Mussel 2'), lty = 1,
				col = cols[1:2], bty = 'n')
	} else {
		# Else just plot the 1st mussel data
		plot(df3$DateTimePDT[sr:er], df3$Temp1calib[sr:er], type = 'l', 
				ylab = expression(Temperature*','~degree*C),
				xlab = 'Time',
				col = cols[1],
				ylim = ylims,
				las = 1)
	}
}	# end of plotTemps function

##################################################################
# Function to plot one channel of Hall effect data at a time
plotHall = function(df3, Ch = 1, startT = '2015-07-15 07:30', 
		endT = '2015-08-06 9:20'){
	cols = brewer.pal(6,'Set1')
	# Establish start and end times for the plot
	startT = as.POSIXct(startT)
	endT = as.POSIXct(endT)
	# Find nearest row indices in the data set
	startTind = which.min(abs(df3$DateTimePDT - startT))
	endTind = which.min(abs(df3$DateTimePDT - endT))
	ylims = range(-5,110)
	if (Ch == 1) {
		plot(df3$DateTimePDT[startTind:endTind], 
				df3$Hall1.percent[startTind:endTind], type = 'n',
				ylab = 'Gape width, %',
				xlab = 'Time',
				col = cols[1],
				ylim = ylims,
				las = 1,
				main = 'Mussel 2')
		grid()
		lines(df3$DateTimePDT[startTind:endTind],
				df3$Hall1.percent[startTind:endTind], col = cols[1])
	} else if (Ch == 2) {
		plot(df3$DateTimePDT[startTind:endTind], 
				y = df3$Hall2.percent[startTind:endTind], type = 'n',
				ylab = 'Gape width, %',
				xlab = 'Time',
				col = cols[2],
				ylim = ylims,
				las = 1,
				main = 'Mussel 2')
		grid()
		lines(df3$DateTimePDT[startTind:endTind],
				df3$Hall2.percent[startTind:endTind], col = cols[2])
	}
}
################################################################################
## Plot the accelerometer + magnetometer data

plotAccelMag = function(df, Ch = 1, startT = '2015-07-15 07:30', 
		endT = '2015-08-06 9:20'){
	cols = brewer.pal(6,'Set1')
	# Establish start and end times for the plot
	startT = as.POSIXct(startT)
	endT = as.POSIXct(endT)
	# Find nearest row indices in the data set
	startTind = which.min(abs(df$DateTimePDT - startT))
	endTind = which.min(abs(df$DateTimePDT - endT))
	yrange = c(-32767, 32767) # The accel/mag sensor range should be
								# the limits of a 16-bit signed integer
	if (Ch == 1) {
		plot(df$DateTimePDT[startTind:endTind], df$a1.x[startTind:endTind],
				type = 'n',
				ylim = yrange,
				ylab = 'Accelerometer', 
				xlab = 'Time',
				main = 'Mussel 1')
		rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4],col = 'grey90')
		grid()
		lines(df$DateTimePDT[startTind:endTind], df$a1.x[startTind:endTind], 
				col = cols[1])
		lines(df$DateTimePDT[startTind:endTind], df$a1.y[startTind:endTind], 
				col = cols[2])
		lines(df$DateTimePDT[startTind:endTind], df$a1.z[startTind:endTind], 
				col = cols[3])
		lines(df$DateTimePDT[startTind:endTind], df$m1.x[startTind:endTind], 
				col = cols[4])
		lines(df$DateTimePDT[startTind:endTind], df$m1.y[startTind:endTind], 
				col = cols[5])
		lines(df$DateTimePDT[startTind:endTind], df$m1.z[startTind:endTind], 
				col = cols[6])
		legend('top', legend = c('accel.x','accel.y','accel.z',
						'mag.x','mag.y','mag.z'),
				col = cols, lty = 1, lwd = 2, bty = 'n', horiz = TRUE)
	}
	else if (Ch == 2) {
		plot(df$DateTimePDT[startTind:endTind], df$a2.x[startTind:endTind],
				type = 'n',
				ylim = yrange,
				ylab = 'Accelerometer', 
				xlab = 'Time',
				main = 'Mussel 2')
		rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4],col = 'grey90')
		grid()
		lines(df$DateTimePDT[startTind:endTind], df$a2.x[startTind:endTind], 
				col = cols[1])
		lines(df$DateTimePDT[startTind:endTind], df$a2.y[startTind:endTind], 
				col = cols[2])
		lines(df$DateTimePDT[startTind:endTind], df$a2.z[startTind:endTind], 
				col = cols[3])
		lines(df$DateTimePDT[startTind:endTind], df$m2.x[startTind:endTind], 
				col = cols[4])
		lines(df$DateTimePDT[startTind:endTind], df$m2.y[startTind:endTind], 
				col = cols[5])
		lines(df$DateTimePDT[startTind:endTind], df$m2.z[startTind:endTind], 
				col = cols[6])
		legend('top', legend = c('accel.x','accel.y','accel.z',
						'mag.x','mag.y','mag.z'),
				col = cols, lty = 1, lwd = 2, bty = 'n', horiz = TRUE)
	}
}


################################################################################
saveMe = FALSE
# Only execute this when you're really ready to save the cleaned data
if (saveMe == TRUE){

	# Write out the 1-Hz temperature and hall effect data to its own file.
	write.csv(df3, file = paste0(sn,'_all_TempHall_filtered.csv'), 
			row.names = FALSE)
	# Write out the 4-Hz accel/magnetometer data to its own file. 
	write.csv(df2, file = paste0(sn,'_all_AccelMag_raw.csv'),
			row.names = FALSE)
	cat('\aSaved files as\n')
	cat(paste0(sn,'_all_TempHall_filtered.csv'),'\n')
	cat(paste0(sn,'_all_AccelMag_raw.csv'),'\n')
}
################################################################################







