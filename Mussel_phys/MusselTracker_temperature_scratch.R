# MusselTracker_temperature_scratch.R
# 
# Author: Luke Miller Sep 2, 2015
###############################################################################
# Experimenting with the temperature data files from Summer 2015
library(zoo)
library(RColorBrewer)
# Set the current working environment to use Pacific Standard/Daylight time zone
# since the MusselTracker data timestamps are all local daylight savings time
Sys.setenv(TZ='PST8PDT')

##############################################################################
idir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/ibuttontemps/"

ibuttdf = read.csv(paste0(idir,'ibutton_mussels_20150725-20150806.csv'),
		colClasses = c('POSIXct',rep('numeric',5)))

##############################################################################
# Load tide height data
tidedir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Field_data/"
tides = read.csv(paste0(tidedir,'Tides_20150714-20150807.csv'),
		colClasses = c('POSIXct','numeric','character','numeric'))
waves = read.csv(paste0(tidedir,'CDIP_158_buoy_data.csv'), 
		colClasses = c('POSIXct',rep('numeric',5)))

###############################################################################
fdir = 'D:/Miller_projects/Mussel_phys/MusselTracker_data/Field_data_concatenated/'

# Low Shore boards
lowshore = c('SN01','SN02','SN03','SN13','SN14','SN15')
# Step through each board serial number and load the data into a data frame
# with the board's serial number as the name
# This will take several minutes

#for (i in 1:length(lowshore)){
#	assign(lowshore[i], 
#			as.data.frame(read.csv(paste0(fdir,lowshore[i],
#									'_all_TempHall_filtered.csv'), 
#			colClasses = c('numeric','POSIXct','factor','numeric','integer',
#					'numeric','integer','numeric','numeric','integer','integer',
#					'numeric','numeric'))))
#}
#
#lowdf = data.frame(DateTimePDT = SN01$DateTimePDT,
#		SN01Temp1calib = SN01$Temp1calib,
#		SN01Temp2calib = SN01$Temp2calib,
#		SN02Temp1calib = SN02$Temp1calib,
#		SN02Temp2calib = SN02$Temp2calib,
#		SN03Temp1calib = SN03$Temp1calib,
#		SN03Temp2calib = SN03$Temp2calib,
#		SN13Temp1calib = SN13$Temp1calib,
#		SN13Temp2calib = SN13$Temp2calib,
#		SN14Temp1calib = SN14$Temp1calib,
#		SN14Temp2calib = SN14$Temp2calib,
#		SN15Temp1calib = SN15$Temp1calib,
#		SN15Temp2calib = SN15$Temp2calib)

# High Shore boards
highshore = c('SN04','SN05','SN06','SN10','SN11','SN12')
# Step through each board serial number and load the data into a data frame
# with the board's serial number as the name
# This will take several minutes

#for (i in 1:length(highshore)){
#	assign(highshore[i], 
#			as.data.frame(read.csv(paste0(fdir,highshore[i],
#									'_all_TempHall_filtered.csv'), 
#				colClasses = c('numeric','POSIXct','factor','numeric','integer',
#					'numeric','integer','numeric','numeric','integer','integer',
#							'numeric','numeric'))))
#}
#
#highdf = data.frame(DateTimePDT = SN04$DateTimePDT,
#		SN04Temp1calib = SN04$Temp1calib,
#		SN04Temp2calib = SN04$Temp2calib,
#		SN05Temp1calib = SN05$Temp1calib,
#		SN05Temp2calib = SN05$Temp2calib,
#		SN06Temp1calib = SN06$Temp1calib,
#		SN06Temp2calib = SN06$Temp2calib,
#		SN10Temp1calib = SN10$Temp1calib,
#		SN10Temp2calib = SN10$Temp2calib,
#		SN11Temp1calib = SN11$Temp1calib,
#		SN11Temp2calib = SN11$Temp2calib,
#		SN12Temp1calib = SN12$Temp1calib,
#		SN12Temp2calib = SN12$Temp2calib)

########################################################################
# Tidepool boards
tidepool = c('SN07','SN08','SN09')

# Step through each board serial number and load the data into a data frame
# with the board's serial number as the name
# This will take several minutes

#for (i in 1:length(tidepool)){
#	assign(tidepool[i], 
#			as.data.frame(read.csv(paste0(fdir,tidepool[i],
#									'_all_TempHall_filtered.csv'), 
#				colClasses = c('numeric','POSIXct','factor','numeric','integer',
#					'numeric','integer','numeric','numeric','integer','integer',
#					'numeric','numeric'))))
#}
#
#pooldf = data.frame(DateTimePDT = SN07$DateTimePDT,
#		SN07Temp1calib = SN07$Temp1calib,
#		SN07Temp2calib = SN07$Temp2calib,
#		SN08Temp1calib = SN08$Temp1calib,
#		SN08Temp2calib = SN08$Temp2calib,
#		SN09Temp1calib = SN09$Temp1calib,
#		SN09Temp2calib = SN09$Temp2calib)

##############################################################################
# Optionally, just load the combined data frames for each plate (high,low,pool)
# from the csv files I made
bigdir = "D:/Miller_projects/Mussel_phys/MusselTracker_data/Working_files/"
lowdf = read.csv(paste0(bigdir,'Low_shore_temps.csv'),
		colClasses = c('POSIXct',rep('numeric',12)))
highdf = read.csv(paste0(bigdir,"High_shore_temps.csv"),
		colClasses = c('POSIXct',rep('numeric',12)))
pooldf = read.csv(paste0(bigdir,'Tidepool_temps.csv'),
		colClasses = c('POSIXct',rep('numeric',6)))

##############################################################################


##########################################
# Convert lowshore values to a zoo object
lowdfz = zoo(lowdf[,-1], order.by = lowdf[,'DateTimePDT'])

# Daily aggregation. Using as.Date strips off the time values from z, leaving 
# only date values that aggregate will use. The data for each plate remain
# separate in their own columns
lowz.day = aggregate(lowdfz,as.Date(time(lowdfz)), FUN = mean, na.rm = TRUE)
lowz.daily.max = aggregate(lowdfz,as.Date(time(lowdfz)), max, na.rm = TRUE)
lowz.daily.min = aggregate(lowdfz,as.Date(time(lowdfz)), min, na.rm = TRUE)

# Replace any -Inf values that have crept in with NAs
for (i in 1:ncol(lowz.daily.max)){
	lowz.daily.max[!is.finite(lowz.daily.max[,i]),i] = NA	
}
for (i in 1:ncol(lowz.day)){
	lowz.day[!is.finite(lowz.day[,i]),i] = NA	
}
for (i in 1:ncol(lowz.daily.min)){
	lowz.daily.min[!is.finite(lowz.daily.min[,i]),i] = NA	
}

# For each day, calculate the range of maximum temperature achieved by the
# mussels on that plate. 
lowz.daily.max.range = apply(lowz.daily.max,MARGIN = 1, FUN = range, na.rm=TRUE)
# Put the daily max, min, and range in a data frame. The values stored here
# are the lowest maximum, lowest maximum, and difference between those values
# within a given day for the mussels on the low shore plate. 
lowdftemps = data.frame(Date = index(lowz.daily.max),
		MaxTemp = lowz.daily.max.range[2,],
		MinTemp = lowz.daily.max.range[1,],
		TempRange = lowz.daily.max.range[2,] - lowz.daily.max.range[1,])



# Calculate daily temperature range by subtracting max and min values
lowz.daily.range = lowz.daily.max - lowz.daily.min
# Calculate standard deviation of temperatures within a day, for each separate 
# plate
#lowz.daily.sd = aggregate(lowdfz,as.Date(time(lowdfz)), FUN = sd, na.rm = TRUE)
#lowz.daily.var = aggregate(lowdfz,as.Date(time(lowdfz)), FUN = var, na.rm = TRUE)

# Calculate overall mean for each plate during the whole experiment
low.av.daily = colMeans(lowz.day, na.rm = TRUE)
low.av.daily.max = colMeans(lowz.daily.max, na.rm = TRUE) 
low.av.daily.min = colMeans(lowz.daily.min, na.rm = TRUE)
low.av.daily.range = colMeans(lowz.daily.range, na.rm = TRUE)
#av.daily.sd = colMeans(lowz.daily.sd, na.rm = TRUE)
#av.daily.var = colMeans(lowz.daily.var, na.rm = TRUE)


# Calculate total degree hours
low.tot.degree = colSums(lowdfz, na.rm = TRUE) # sum each column's temperatures
low.tot.degree = (low.tot.degree  / 60) / 60 # divide by 60 sec and divide by 60 min.

# Calculate the number of missing data points in the data set. This has entries
# for every column in z.
low.nas = which(is.na(lowdfz))
# Calculate the total number of temperature records in the data set, excluding
# missing data (NA's). 
low.totaltemps = (nrow(lowdfz) * ncol(lowdfz)) - length(low.nas)

#######################################################################
# Convert highshore values to a zoo object
highdfz = zoo(highdf[,-1], order.by = highdf[,'DateTimePDT'])

# Daily aggregation. Using as.Date strips off the time values from z, leaving 
# only date values that aggregate will use. The data for each plate remain
# separate in their own columns
highz.day = aggregate(highdfz,as.Date(time(highdfz)), FUN = mean, na.rm = TRUE)
highz.daily.max = aggregate(highdfz,as.Date(time(highdfz)), max, na.rm = TRUE)
highz.daily.min = aggregate(highdfz,as.Date(time(highdfz)), min, na.rm = TRUE)

# Replace any -Inf values that have crept in with NAs
for (i in 1:ncol(highz.daily.max)){
	highz.daily.max[!is.finite(highz.daily.max[,i]),i] = NA	
}
for (i in 1:ncol(highz.day)){
	highz.day[!is.finite(highz.day[,i]),i] = NA	
}
for (i in 1:ncol(highz.daily.min)){
	highz.daily.min[!is.finite(highz.daily.min[,i]),i] = NA	
}

# For each day, calculate the range of maximum temperature achieved by the
# mussels on that plate. 
highz.daily.max.range = apply(highz.daily.max,MARGIN = 1, FUN = range, na.rm=TRUE)
# Put the daily max, min, and range in a data frame. The values stored here
# are the highest maximum, lowest maximum, and difference between those values
# within a given day for the mussels on the high shore plate. 
highdftemps = data.frame(Date = index(highz.daily.max),
		MaxTemp = highz.daily.max.range[2,],
		MinTemp = highz.daily.max.range[1,],
		TempRange = highz.daily.max.range[2,] - highz.daily.max.range[1,])

# Calculate daily temperature range by subtracting max and min values
# This will primarily give the range of temperatures between ocean temperature
# and the daily max temperature (but cold nights may be the minimum on some 
# days). This still gives values for individual mussels
highz.daily.range = highz.daily.max - highz.daily.min
# Calculate standard deviation of temperatures within a day, for each separate 
# plate
#highz.daily.sd = aggregate(highdfz,as.Date(time(highdfz)), FUN = sd, na.rm = TRUE)
#highz.daily.var = aggregate(highdfz,as.Date(time(highdfz)), FUN = var, na.rm = TRUE)

# Calculate overall mean for each mussel during the whole experiment
high.av.daily = colMeans(highz.day, na.rm = TRUE)
high.av.daily.max = colMeans(highz.daily.max, na.rm = TRUE) 
high.av.daily.min = colMeans(highz.daily.min, na.rm = TRUE)
high.av.daily.range = colMeans(highz.daily.range, na.rm = TRUE)
#av.daily.sd = colMeans(highz.daily.sd, na.rm = TRUE)
#av.daily.var = colMeans(highz.daily.var, na.rm = TRUE)


# Calculate total degree hours
high.tot.degree = colSums(highdfz, na.rm = TRUE) # sum each column's temperatures
high.tot.degree = (high.tot.degree  / 60) / 60 # divide by 60 sec and divide by 60 min.

# Calculate the number of missing data points in the data set. This has entries
# for every column in z.
high.nas = which(is.na(highdfz))
# Calculate the total number of temperature records in the data set, excluding
# missing data (NA's). 
high.totaltemps = (nrow(highdfz) * ncol(highdfz)) - length(high.nas)

#######################################################################
# Convert poolshore values to a zoo object
pooldfz = zoo(pooldf[,-1], order.by = pooldf[,'DateTimePDT'])

# Daily aggregation. Using as.Date strips off the time values from z, leaving 
# only date values that aggregate will use. The data for each plate remain
# separate in their own columns
poolz.day = aggregate(pooldfz,as.Date(time(pooldfz)), FUN = mean, na.rm = TRUE)
poolz.daily.max = aggregate(pooldfz,as.Date(time(pooldfz)), max, na.rm = TRUE)
poolz.daily.min = aggregate(pooldfz,as.Date(time(pooldfz)), min, na.rm = TRUE)

# Replace any -Inf values that have crept in with NAs
for (i in 1:ncol(poolz.daily.max)){
	poolz.daily.max[!is.finite(poolz.daily.max[,i]),i] = NA	
}
for (i in 1:ncol(poolz.day)){
	poolz.day[!is.finite(poolz.day[,i]),i] = NA	
}
for (i in 1:ncol(poolz.daily.min)){
	poolz.daily.min[!is.finite(poolz.daily.min[,i]),i] = NA	
}

# For each day, calculate the range of maximum temperature achieved by the
# mussels on that plate. 
poolz.daily.max.range = apply(poolz.daily.max,MARGIN = 1, FUN = range, na.rm=TRUE)
# Put the daily max, min, and range in a data frame. The values stored here
# are the poolest maximum, lowest maximum, and difference between those values
# within a given day for the mussels on the pool shore plate. 
pooldftemps = data.frame(Date = index(poolz.daily.max),
		MaxTemp = poolz.daily.max.range[2,],
		MinTemp = poolz.daily.max.range[1,],
		TempRange = poolz.daily.max.range[2,] - poolz.daily.max.range[1,])

# Calculate daily temperature range by subtracting max and min values
# This will primarily give the range of temperatures between ocean temperature
# and the daily max temperature (but cold nights may be the minimum on some 
# days). This still gives values for individual mussels
poolz.daily.range = poolz.daily.max - poolz.daily.min

# Calculate overall mean for each plate during the whole experiment
pool.av.daily = colMeans(poolz.day, na.rm = TRUE)
pool.av.daily.max = colMeans(poolz.daily.max, na.rm = TRUE) 
pool.av.daily.min = colMeans(poolz.daily.min, na.rm = TRUE)
pool.av.daily.range = colMeans(poolz.daily.range, na.rm = TRUE)
#av.daily.sd = colMeans(poolz.daily.sd, na.rm = TRUE)
#av.daily.var = colMeans(poolz.daily.var, na.rm = TRUE)


# Calculate total degree hours
pool.tot.degree = colSums(pooldfz, na.rm = TRUE) # sum each column's temperatures
pool.tot.degree = (pool.tot.degree  / 60) / 60 # divide by 60 sec and divide by 60 min.

# Calculate the number of missing data points in the data set. This has entries
# for every column in z.
pool.nas = which(is.na(pooldfz))
# Calculate the total number of temperature records in the data set, excluding
# missing data (NA's). 
pool.totaltemps = (nrow(pooldfz) * ncol(pooldfz)) - length(pool.nas)
##########################################################################
# Do the same for ibutton mussel mimics
ibuttdfz = zoo(ibuttdf[,-1], order.by = ibuttdf[,'DateTimePDT'])
# Daily aggregation. Using as.Date strips off the time values from z, leaving 
# only date values that aggregate will use. The data for each plate remain
# separate in their own columns
ibuttz.day = aggregate(ibuttdfz,as.Date(time(ibuttdfz)), FUN = mean, na.rm = TRUE)
ibuttz.daily.max = aggregate(ibuttdfz,as.Date(time(ibuttdfz)), max, na.rm = TRUE)
ibuttz.daily.min = aggregate(ibuttdfz,as.Date(time(ibuttdfz)), min, na.rm = TRUE)

# Replace any -Inf values that have crept in with NAs
for (i in 1:ncol(ibuttz.daily.max)){
	ibuttz.daily.max[!is.finite(ibuttz.daily.max[,i]),i] = NA	
}
for (i in 1:ncol(ibuttz.day)){
	ibuttz.day[!is.finite(ibuttz.day[,i]),i] = NA	
}
for (i in 1:ncol(ibuttz.daily.min)){
	ibuttz.daily.min[!is.finite(ibuttz.daily.min[,i]),i] = NA	
}

# For each day, calculate the range of maximum temperature achieved by the
# mussels on that plate. 
ibuttz.daily.max.range = apply(ibuttz.daily.max,MARGIN = 1, FUN = range, na.rm=TRUE)
# Put the daily max, min, and range in a data frame. The values stored here
# are the ibuttest maximum, lowest maximum, and difference between those values
# within a given day for the mussels on the ibutt shore plate. 
ibuttdftemps = data.frame(Date = index(ibuttz.daily.max),
		MaxTemp = ibuttz.daily.max.range[2,],
		MinTemp = ibuttz.daily.max.range[1,],
		TempRange = ibuttz.daily.max.range[2,] - ibuttz.daily.max.range[1,])
# Calculate daily temperature range by subtracting max and min values
# This will primarily give the range of temperatures between ocean temperature
# and the daily max temperature (but cold nights may be the minimum on some 
# days). This still gives values for individual mussels
ibuttz.daily.range = ibuttz.daily.max - ibuttz.daily.min

# Calculate overall mean for each plate during the whole experiment
ibutt.av.daily = colMeans(ibuttz.day, na.rm = TRUE)
ibutt.av.daily.max = colMeans(ibuttz.daily.max, na.rm = TRUE) 
ibutt.av.daily.min = colMeans(ibuttz.daily.min, na.rm = TRUE)
ibutt.av.daily.range = colMeans(ibuttz.daily.range, na.rm = TRUE)
#av.daily.sd = colMeans(ibuttz.daily.sd, na.rm = TRUE)
#av.daily.var = colMeans(ibuttz.daily.var, na.rm = TRUE)


# Calculate total degree hours
ibutt.tot.degree = colSums(ibuttdfz, na.rm = TRUE) # sum each column's temperatures
ibutt.tot.degree = (ibutt.tot.degree  / 60) / 60 # divide by 60 sec and divide by 60 min.

# Calculate the number of missing data points in the data set. This has entries
# for every column in z.
ibutt.nas = which(is.na(ibuttdfz))
# Calculate the total number of temperature records in the data set, excluding
# missing data (NA's). 
ibutt.totaltemps = (nrow(ibuttdfz) * ncol(ibuttdfz)) - length(ibutt.nas)

####################################################
# Plot high shore temps
cols = brewer.pal(12,'Set3')
sn = which.min(abs(highdf[,1] - as.POSIXct('2015-07-31 00:00')))
en = which.min(abs(highdf[,1] - as.POSIXct('2015-08-06 09:00')))
plotcols = c(2,12) # which data columns to plot
ylims = range(highdf[sn:en,2:ncol(highdf)],na.rm = TRUE)
plot(highdf[sn:en,1],
		y = highdf[sn:en,2],
		type = 'n',
		ylim = c(0,40),
		las = 1,
		ylab = expression(Temperature*','~degree*C),
		xlab = 'Date',
		main = 'High shore')
rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4], col = 'grey30')
grid(col = 'white')
for (i in plotcols[1]:plotcols[2]){
	lines(highdf[sn:en,1], y = highdf[sn:en,i],
			col = cols[i-1])
}
# Plot tide height at the bottom
lines(tides$DateTimePDT,tides$TideHT.m, col = 'blue')
lines(waves$DateTimePDT,waves$Hs.m, col = 'green')
legend('topleft',legend = colnames(lowdf[plotcols[1]:plotcols[2]]), 
		col = cols[(plotcols[1]:plotcols[2])-1], lty = 1,
		lwd = 2)

#####################################################################
# Plot low shore temps
sn = which.min(abs(lowdf[,1] - as.POSIXct('2015-08-01 12:00')))
en = which.min(abs(lowdf[,1] - as.POSIXct('2015-08-01 20:00')))
plotcols = c(2,12) # which data columns to plot
ylims = range(lowdf[,2:ncol(lowdf)],na.rm = TRUE)
plot(lowdf[sn:en,1],
		y = lowdf[sn:en,2],
		type = 'n',
		ylim = c(-1,45),
		las = 1,
		ylab = expression(Temperature*','~degree*C),
		xlab = 'Date',
		main = 'Low shore')
rect(par()$usr[1],par()$usr[3],par()$usr[2],par()$usr[4], col = 'grey30')
grid(col = 'white')
for (i in plotcols[1]:plotcols[2]){
	lines(lowdf[sn:en,1], y = lowdf[sn:en,i],
			col = cols[i-1])
}
# Plot tide height at the bottom
lines(tides$DateTimePDT,tides$TideHT.m, col = 'blue')
lines(waves$DateTimePDT,waves$Hs.m, col = 'green')
legend('topleft',legend = colnames(lowdf[plotcols[1]:plotcols[2]]), 
		col = cols[(plotcols[1]:plotcols[2])-1], lty = 1,
		lwd = 2)