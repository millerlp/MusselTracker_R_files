# Thermocouple_calibration.R
# 
# Author: Luke Miller Jul 24, 2015
###############################################################################

tspan = 100 # number of samples to extract for each temperature setting

tcdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Thermocouple_calibrations/"
tcdir = paste0(tcdir,'20150630_calibrations/')

fnames = dir(tcdir, pattern = '*.csv', full.names = TRUE)

# Open the file containing time points and reference thermometer temperatures
times = read.table(paste0(tcdir,'calib_time_temp2.txt'))
times$DateTime = as.POSIXct(paste(times$V1,times$V2),origin = '1970-1-1',
		format = "%m/%d/%Y %H:%M")
times = times[,c(4,3)]
names(times)[2] = 'TempC'

calibs = data.frame(SerialNumber = character(0), 
		Channel = numeric(0), Intercept = numeric(0), Slope = numeric(0))

for (i in 1:length(fnames)){
	# Open each file in turn, reading the 1st line and parsing the serial number
	info = scan(fnames[i], nlines = 1, what = 'character', sep = ',')
	sn = as.character(info[2]) # grab just the serial number
	# Now re-open and read the file, skipping the 1st header line
	df = read.csv(fnames[i], skip = 1)
	df$DateTime = as.POSIXct(df$DateTime, origin  = '1970-1-1')
	
	# Now go through the two data frames, times and df, find the matching
	# time stamps, extract a set of temperature readings from df, and regress
	# those against the known temperatures in the 'times' data frame.
	df2 = data.frame(RefTemp = numeric(0),TCTemp1 = numeric(0),
			TCTemp2=numeric(0))
	
	for (j in 1:nrow(times)){
		# Find closest time point in df that matches times$DateTime[j]
		indx = which.min(abs(df$DateTime - times$DateTime[j]))
		reftemp = times[j, 'TempC'] # Extract the reference temp
		meastemps1 = df[indx:(indx+(tspan-1)),'Temp1']
		meastemps2 = df[indx:(indx+(tspan-1)),'Temp2']
		tempdf = data.frame(RefTemp = rep(reftemp, tspan), 
				TCTemp1 = meastemps1, TCTemp2 = meastemps2)
		df2 = rbind(df2,tempdf)
	}
	# At this point df2 should contain a series of RefTemp values and the
	# associated thermocouple temperatures at the matching time points. 
	# Now we can regress those to get an equation for the line that relates
	# measured temperature to actual (reference) temperature.
	tc1 = lm(RefTemp ~ TCTemp1, data = df2)
	tc2 = lm(RefTemp ~ TCTemp2, data = df2)
	
	tempdf2 = data.frame(Serial = rep(sn,2), Channel = c(1,2),
			Intercept = c(tc1$coefficients[1], tc2$coefficients[1]),
			Slope = c(tc1$coefficients[2], tc2$coefficients[2]))
	calibs = rbind(calibs,tempdf2)
}

write.csv(calibs,file = paste0(tcdir,'Calibration_coefficients_20150630.csv'),
		row.names = FALSE)