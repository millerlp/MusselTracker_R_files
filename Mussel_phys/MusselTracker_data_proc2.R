# MusselTracker_data_proc.R
# Procedures to concatenate data files from a MusselTracker board into a single
# continuous time series, including filtering of data, removing bad data, 
# making small interpolations for small gaps (less than 30 seconds) etc. 
# Author: Luke Miller Jul 9, 2015
###############################################################################
require(zoo)
require(RColorBrewer)
require(signal)
require(zoo)

Sys.setenv(TZ = 'PST8PDT')	# Set the time zone for the mussel tracker data
Sys.getenv('TZ')

## Show fractional seconds when printing to terminal. This still often doesn't
## properly show the fractional seconds when displaying the formatted time
options(digit.secs = 6, digits = 12) 

board = 'SN14'

# Specify a start and end time cutoff for the dataset. The plates were 
# deployed in the field by 2015-07-15 7:30PDT, prior points are in the water table
cutoffstart = as.POSIXct('2015-07-15 0:00:00', tz = 'PST8PDT')
cutoffend = as.POSIXct('2015-08-06 09:20', tz = 'PST8PDT')

# Specify the directory holding the field data files
fdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Field_data/"
# Specify the directory holding the thermocouple calibration data file
tcdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Thermocouple_calibrations/20150630_calibrations/"



fname1 = paste0(fdir,board)
fnames = dir(fname1, pattern = '*.csv', full.names = TRUE)

#fname = '20150710_0000_00_SN02.csv'

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
# The data frame df should already by ordered. We can create a larger data
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

# Remove data from before the field deployment if present.
df2 = df2[df2$DateTime >= cutoffstart,]

# Overwrite df with the new df2 version
df = df2

###############################################################################
# Recall that the temperature and hall effect data were only sampled 1 time
# per second, with the same value being written to all 4 sub-second samples.
# As a result, we ought to pull out the "true" 1-second values, resulting in a 
# much smaller data set, and also leave the accelerometer/magnetomter values
# out, since they were sampled at 4Hz. 

df3 = df[seq(1,(nrow(df)-3), by = 4),c('POSIXt','DateTime','SerialNumber',
				'Temp1','Hall1', 'Temp2','Hall2')]




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

df3$Temp1c = applyTCcalib(df3$Temp1, tc1int, tc1sl)
df3$Temp2c = applyTCcalib(df3$Temp2, tc2int, tc2sl)

#######*********************************************************************
# There will be small gaps in the dataset of a few seconds, along with some
# longer gaps. Use the functions in the zoo package to fill in the brief 
# missing data with the prior values, defining a maxgap above which values
# won't be filled in. 
cat('Filling small gaps in the dataset\n')
# Convert the data frame to a zoo object, leaving behind the POSIXt, DateTime
# and SerialNumber columns so that only numeric data are included
z = zoo(df3[, -c(1,2,3)], order.by = df3$POSIXt)
# For brief missing data, carry forward the last non-NA value, up to a maximum
# of maxgap rows. For these 1-second data, 30 = thirty second gaps max.
z2 = na.locf(z, maxgap = 10)
# Copy the filled data back into df3
df3[,4:9] = z2
# At this point there should only be significant gaps of missing data during
# times when the datalogger was not recording for some reason, due to battery
# failure, interruptions for downloading, etc. 


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
	# input vector. It does not currently handle data streams that start off
	# with NA values, but should handle streams that end with NAs. 
	
	# Use run length encoding function to find NA gaps
	gaps = rle(is.na(values))
	# The number of gaps with value = TRUE is the number of good runs of data
	numgaps = sum(!gaps$values)
	# Find the indices of the entries that are FALSE (negated to TRUE) that
	# represent the start and end points of good data
	truegaps = which(!gaps$values)
	# Create the output data frame
	results = data.frame(Start = integer(numgaps), End = integer(numgaps))
	# The first entry should always be 1 (1st row). If there are no gaps, the
	# contents of gaps$lengths will just be the length of the values vector
	results$Start[1] = 1
	results$End[1] = gaps$lengths[1]
	# If there is more than 1 entry in gaps$lengths, process the other gaps
	j = 2; # counter
	if (numgaps> 1){
		if (numgaps %% 2 == 0){
			# Even number of gaps, dataset ends on NAs
			for (i in seq(2, length(gaps$lengths)-1, by = 2)){
				nextstart = sum(gaps$lengths[1:i]) + 1
				nextend = sum(gaps$lengths[1:(i+1)])
				results$Start[j] = nextstart
				results$End[j] = nextend
				j = j + 1
			}
			
		} else if (numgaps %% 2 != 0){
			# Odd number of gaps
			for (i in seq(2,length(gaps$lengths), by = 2)){
				# The first value in gaps$lengths should represent the end of a
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
		yfiltered = yfiltered +offset
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
# x percentile of the hall effect data. This will hopefully remove the influence
# of any single spurious high or low readings. 
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
# Supply a vector of Hall effect sensor data
calcPercentGape = function (hallData){
	# First calculate the upper and lower 1% values of the raw Hall readings
	myrange = percentileRange(hallData, percentileLim = 0.01)
	if (myrange[2] < 512){
		# If the hall effect and magnet are situated so that the value goes down 
		# as the mussel closes, then calculate % opening using the maximum 
		# value. If the maximum value is less than 512, we'll assume that 
		# closing drives the value down.
		# Calculate the width of the range of values. If the readings are less
		# than 512, the smaller value in meanVals represents the fully-closed 
		# reading while the larger value is the fully open reading. 
		rangeval = myrange[2] - myrange[1]
		# Calculate a percentage opening based on the difference between each 
		# reading and the minimum value (fully closed), divided by the range
		Hall.percent = ((hallData - myrange[1]) / 
					rangeval) * 100

	} else if (myrange[2] > 512){
		# Calculate the width of the range of values. If the readings are larger
		# than 512, the smaller value in myrange represents the fully-open 
		# reading while the larger value is the fully-closed reading. 
		rangeval = myrange[2] - myrange[1]
		# Calculate a percentage opening based on the difference between each 
		# reading and the maximum value (fully closed), divided by the range
		Hall.percent = ((myrange[2] - hallData) / 
					rangeval) * 100
	}
	Hall.percent	# return the set of values
}



# TODO: make a version of the above that can split datasets and deal with 
# times where a magnetic sensor had to be repositioned, creating a new 
# baseline. For board SN12 Hall1 this is necessary. 

cat('Filtering Hall effect data\n')
df3$Hall1filt = hallFilter(df3$Hall1)
df3$Hall2filt = hallFilter(df3$Hall2)

df3$Hall1.percent = calcPercentGape(df3$Hall1filt)
df3$Hall2.percent = calcPercentGape(df3$Hall2filt)

# Note that some of the sensors need to have portions of their data excised 
# due to spurious values resulting from sensor dislodgement or other failures
if (board == 'SN12'){
	# Enter a set of starting and ending time values that you want to cut
	# out of the hall effect data for a particular sensor
	timevals = c( # 	 start times	 	end times
					'2015-07-17 00:00', '2015-07-18 12:00',
					'2015-07-20 10:00', '2015-07-20 14:00',
					'2015-07-20 18:00',	'2015-07-20 21:00',
					'2015-07-22 07:00', '2015-07-22 12:00'
					)
	timevals = as.POSIXct(timevals) # convert to timestamps
	# Find the closest matching row in df3 by searching for the smallest time
	# difference between each timestamp in timevals and each DateTime in df3
	for (i in 1:length(timevals)){
		timevals2[i] = which.min(abs(timevals[i] - df3$DateTime))	
	}
	# Create a data frame by taking the start values and end values of 
	# row indices and splitting them into two columns
	excisedf = data.frame(Start = timevals2[seq(1,length(timevals2),by=2)],
			End = timevals2[seq(2,length(timevals2),by=2)])
	
	# Go through each set of row indices in excisedf, overwrite spurious data
	# as NA's.
	for (i in 1:nrow(excisedf)){
		# Convert all suspect data to NAs
		st = excisedf$Start[i]
		end = excisedf$End[i]
		df3$Hall1[st:end] = NA
	}
	# Apply the butterworth filter
	df3$Hall1filt = hallFilter(df3$Hall1)
	
	# Get the row indices for the major gaps that exist now
	mygaps = gapBounds(df3$Hall1filt)
	# For each entry in mygaps, use the calcPercentGape function to calculate
	# the percent gape opening. This will calculate a new baseline closed 
	# value after every major gap, which often represents a repositioned 
	# magnet or sensor. 
	df3$Hall1.percent = NA	# initialize output column
	for (i in 1:nrow(mygaps)){
		st = mygaps$Start[i]
		end = mygaps$End[i]
		df3$Hall1.percent[st:end] = calcPercentGape(df3$Hall1filt[st:end])
	}	
}

if (board == 'SN11'){
	# SN11 has gape sensor issues on both mussels, probably for separate reasons
	# Mussel 1 also has bad thermocouple data after the 29th. 
	
}

#plot(Hall1.percent~DateTime, data = df3[,], type = 'l')

##***************************************************************************
# Plot the data streams for mussel 1
plot2 = TRUE # set to false to only show one mussel's data

plotTemps = function(df3, start = 1, end = nrow(df3), plot2 = TRUE){
	if (class(start)[1] == 'character'){
		# Handle character timestamps
		start = as.POSIXct(start, origin = '1970-1-1')
		end = as.POSIXct(end, origin = '1970-1-1')
	}
	if (class(start)[1] == 'POSIXct' & class(end)[1] == 'POSIXct'){
		# Find nearest row for each time stamp
		startrow = which.min(abs(start - df3$DateTime))
		endrow = which.min(abs(end - df3$DateTime))
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
		ylims = range(c(df3$Temp1c[sr:er],df3$Temp2c[sr:er]), na.rm = TRUE)
		ylims[2] = ylims[2] + 2 # bump the upper temp limit up slightly
		plot(df3$DateTime[sr:er], df3$Temp1c[sr:er], type = 'n', 
				ylab = expression(Temperature*','~degree*C),
				xlab = 'Time',
				col = cols[1],
				ylim = ylims,
				las = 1)
		grid()
		lines(df3$DateTime[sr:er], df3$Temp1c[sr:er], col = cols[1])
		lines(df3$DateTime[sr:er], df3$Temp2c[sr:er], col = cols[2])	
		legend('topleft', legend = c('Mussel 1','Mussel 2'), lty = 1,
				col = cols[1:2], bty = 'n')
	} else {
		# Else just plot the 1st mussel data
		plot(df3$DateTime[sr:er], df3$Temp1c[sr:er], type = 'l', 
				ylab = expression(Temperature*','~degree*C),
				xlab = 'Time',
				col = cols[1],
				ylim = ylims,
				las = 1)
	}
}	# end of plotTemps function

par(mfrow = c(2,1), mar = c(4,5,1,1))
cols = brewer.pal(3,'Set1')
# Plot temperature first
if (plot2){
	ylims = range(c(df3$Temp1c,df3$Temp2c), na.rm = TRUE)
	ylims[2] = ylims[2] + 2 # bump the upper temp limit up slightly
	plot(df3$DateTime, df3$Temp1c, type = 'n', 
			ylab = expression(Temperature*','~degree*C),
			xlab = 'Time',
			col = cols[1],
			ylim = ylims,
			las = 1)
	grid()
	lines(df3$DateTime, df3$Temp1c, col = cols[1])
	lines(df3$DateTime, df3$Temp2c, col = cols[2])	
	legend('topleft', legend = c('Mussel 1','Mussel 2'), lty = 1,
			col = cols[1:2], bty = 'n')
} else {
	 # Else just plot the 1st mussel data
	plot(df3$DateTime, df3$Temp1c, type = 'l', 
			ylab = expression(Temperature*','~degree*C),
			xlab = 'Time',
			col = cols[1],
			ylim = ylims,
			las = 1)
}

# Plot the hall effect sensor data
if (plot2) {
	ylims = range(-5,110)
	plot(df3$DateTime, df3$Hall1.percent, type = 'n',
			ylab = 'Gape width, %',
			xlab = 'Time',
			col = cols[1],
			ylim = ylims,
			las = 1)
	grid()
	lines(df3$DateTime, df3$Hall1.percent, col = cols[1])
	lines(df3$DateTime, df3$Hall2.percent, col = cols[2])
#	legend('topleft', legend = c('Mussel 1','Mussel 2'), col = cols[1:2],
#			lty = 1)
} else {
	# Else just plot the 1st mussel data
	ylims = range(-5,110)
	plot(df3$DateTime, df3$Hall1.percent, type = 'n',
			ylab = 'Gape width, %',
			xlab = 'Time',
			col = cols[1],
			ylim = ylims,
			las = 1)
	grid()
	lines(df3$DateTime, df3$Hall1.percent, col = cols[1])
}

## Plot the accelerometer + magnetometer data
#yrange = range(c(df$a1.x,df$a1.y,df$a1.z,df$m1.x,df$m1.y,df$m1.z),na.rm=TRUE)
#plot(df$DateTime, df$a1.x, type = 'l',
#		ylim = yrange,
#		ylab = 'Accelerometer',
#		xlab = 'Time'
#		)
#lines(df$DateTime, df$a1.y, col = 2)
#lines(df$DateTime, df$a1.z, col = 3)
#legend('top', legend = c('accel.x','accel.y','accel.z'),
#		col = 1:3, lty = 1, bty = 'n', horiz = TRUE)
#plot(df$DateTime, df$m1.x, type = 'l',
#		col = 4,
#		ylim = yrange,
#		ylab = 'Magnetometer',
#		xlab = 'Time')
#lines(df$DateTime, df$m1.y, col = 5)
#lines(df$DateTime, df$m1.z, col = 6)
#legend('top', legend = c('mag.x','mag.y','mag.z'),
#		col = 4:6, lty = 1, bty = 'n', horiz = TRUE)
#
###****************************
## Plot the data streams for mussel 2
##par(mfrow = c(4,1))
## Plot temperature first
#plot(df$DateTime, df$Temp2, type = 'l', 
#		ylab = 'Temperature, C',
#		xlab = 'Time',
#		col = cols[2],
#		las = 1)
#plot(df$DateTime, df$Hall2.percent, type = 'l',
#		ylab = 'Gape width, %',
#		xlab = 'Time',
#		col = cols[2],
#		las = 1)
## Plot the accelerometer + magnetometer data
mystart = as.POSIXct('2015-07-15 00:00', tz = 'PST8PDT')
myend = as.POSIXct('2015-07-15 9:30', tz = 'PST8PDT')
mystart = which.min(abs(df$DateTime - mystart))
myend = which.min(abs(df$DateTime - myend))

yrange = range(c(df$m1.x,df$m1.y,df$m1.z,df$m1.x,df$m1.y,df$m1.z), na.rm=TRUE)
plot(df$DateTime[mystart:myend], df$m1.x[mystart:myend], type = 'l',
		ylim = yrange,
		ylab = 'Accelerometer',
		xlab = 'Time'
)
lines(df$DateTime[mystart:myend], df$m1.y[mystart:myend], col = 2)
lines(df$DateTime[mystart:myend], df$m1.z[mystart:myend], col = 3)
#legend('top', legend = c('accel.x','accel.y','accel.z'),
#		col = 1:3, lty = 1, bty = 'n', horiz = TRUE)
#plot(df$DateTime, df$m2.x, type = 'l',
#		col = 4,
#		ylim = yrange,
#		ylab = 'Magnetometer',
#		xlab = 'Time')
#lines(df$DateTime, df$m2.y, col = 5)
#lines(df$DateTime, df$m2.z, col = 6)
#legend('top', legend = c('mag.x','mag.y','mag.z'),
#		col = 4:6, lty = 1, bty = 'n', horiz = TRUE)

######################################################################
require(ggplot2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	library(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

#
#mussT = ggplot(df, aes(x = DateTime)) + 
#		geom_line(aes(y = Temp1, color='Mussel 1')) +
#		geom_line(aes(y = Temp2, color = 'Mussel 2'))+
#		labs(x = 'Date', y = 'Temperature, C') + 
#		guides(color = guide_legend(title = 'Mussel')) + 
#		ggtitle("Mussel Temperature")
#
#mussH = ggplot(df, aes(x = DateTime)) + 
#		geom_line(aes(y = Hall1.percent, color = "Mussel 1")) +
#		geom_line(aes(y = Hall2.percent, color = "Mussel 2")) +
#		labs(x = 'Date', y = 'Percent opening') + 
#		guides(color = guide_legend(title = 'Mussel')) +
#		ggtitle("Mussel gaping")
#
#multiplot(mussT, mussH, cols = 1)		
		
#m2T = ggplot(df)

