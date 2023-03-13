#### Code to clean Germ Pheno data #####

library(tidyverse)
library(readxl)
library(lubridate)

# Raw Data File
germdat = read.csv('./Germination.Timing/Raw.Data/Germ Pheno Master Datasheet.csv')

# change column name to match previously used code
germdatx=germdat
germdatx = rename(germdat, germdate = "Germ.date")

##### Fix errors in germination dates #####

germdatx2 = germdatx %>%
  #fix errors in dates.  There are less clunky ways to do this (with fewer lines of code), but the date format makes this longer and clunkier
  mutate(germdate = replace(germdate,germdate =="11/820", "2020-11-08")) %>%
  mutate(germdate = replace(germdate,germdate =="1930-11-26", "2020-11-26")) %>% 
  mutate(germdate = replace(germdate,germdate =="2010-12-15", "2020-12-15")) %>%
  mutate(germdate = replace(germdate,germdate =="2020-01-04", "2021-01-04")) %>% 
  mutate(germdate = replace(germdate,germdate =="2020-01-25", "2021-01-25")) %>%
  mutate(germdate = replace(germdate,germdate =="2021-10-20", "2020-10-20")) %>%
  mutate(germdate = replace(germdate,germdate =="2022-10-20", "2020-10-20")) %>%
  mutate(germdate = replace(germdate,germdate =="2021-12-06", "2020-12-06")) %>%
  mutate(germdate = replace(germdate,germdate =="2021-12-09", "2020-12-09")) %>%
  mutate(germdate = replace(germdate,germdate =="2021-12-28", "2020-12-28"))

# change germdate for row 945 only. Can't use mutate like above because multiple have 2020-10-01 germdate
germdatx2[945,6]="2020-10-11"

# these germ dates were originally missing, but individuals had fitness
# Still 16 individuals with fitness data but no germ data
germdatx2[2202,6]="2020-11-15"
germdatx2[3442,6]="2020-12-24"
germdatx2[3532,6]="2020-12-24"
germdatx2[2045,6]="2020-11-07"
germdatx2[2101,6]="2020-11-08"
germdatx2[2345,6]="2020-11-08"
germdatx2[2627,6]="2020-11-22"
germdatx2[2957,6]="2020-11-22"
germdatx2[2632,6]="2020-11-23"

# format germ date so that we can work with it.
germdatx2 = mutate(germdatx2, germdate = ymd(germdate))

##### Adding planting date to file #####

# load lookup table with planting dates for each cohort
# This data is in the raw data file
planted = read_xlsx("./Germination.Timing/Raw.Data/Fall 2020 germ experiment schedules.xlsx")
planted = mutate(planted, plantdate = ymd(plantdate))

# merge germ data file and planted date file
# both datasets have "Cohort" so it should recognize this
germdaty = left_join(germdatx2, planted, by = c("Cohort"))


# format dates
germdaty2 = germdaty %>%
  mutate(germjul = as.Date(germdate)) %>%  #, origin = "2020-09-01" didn't seem to work
  mutate(germjul = yday(germjul))  %>% #this gives julian date (yday function)
  mutate(plantjul = as.Date(plantdate)) %>%
  mutate(plantjul = yday(plantjul))

# calculate time to germ and time since sept 1
germdaty3 = germdaty2 %>%
  mutate(days2germ = difftime(germdate, plantdate)) %>% #difftime function calculates difference of days between 2 times
  mutate(plantSep = difftime(plantdate, ymd("2020-09-01"))) %>% #this calculates the difference in time from Sept 1 to the planting date
  mutate(germSep = difftime(germdate, ymd("2020-09-01"))) %>%  #this calculates the difference in time from Sept 1 to the germ date
  mutate(days2germ2 = as.numeric(germSep - plantSep)) %>% #checking that the numeric version works well.  It does!
  mutate(planted = 1) %>% #variable that allows us to count each individual in the data
  mutate(germinated = if_else(is.na(germdate) == TRUE, 0,1))

#### Final germination phenology data frame #### 
# can be found in Formatted.Data

write.csv(germdaty3, file="germ.pheno.formatted.csv")

#### Formatting germination phenology Round 2 data ####

# Raw Data File
germdat = read.csv('./Germination.Timing/Raw.Data/Germ_Pheno_Round_2_Germination_Final Copy.csv')

# change column name to match previously used code
germdatx=germdat
germdatx = rename(germdat, germdate = "Germination.Date")

# format germ date so that we can work with it.
germdatx$germdate = as.Date(germdatx$germdate)

##### Fix errors in germination dates #####

germdatx = germdatx %>%
  #fix errors in dates.  There are less clunky ways to do this (with fewer lines of code), but the date format makes this longer and clunkier
  mutate(germdate = replace(germdate,germdate =="2021-05-10", "2021-10-05"))
  
##### Adding planting date to file #####

# load lookup table with planting dates for each cohort
# This data is in the raw data file
planted = read_xlsx("./Germination.Timing/Raw.Data/Fall 2020 germ experiment schedules.xlsx")
planted = mutate(planted, plantdate = ymd(plantdate))

# merge germ data file and planted date file
# change "Old.Cohort in germdatx to "Cohort" so can be merged with plantdate

colnames(germdatx)[2]="Cohort"

# both datasets now have "Cohort" so it should recognize this
germdaty = left_join(germdatx, planted, by = c("Cohort"))

# format dates
germdaty2 = germdaty %>%
  mutate(germjul = as.Date(germdate)) %>%  #, origin = "2020-09-01" didn't seem to work
  mutate(germjul = yday(germjul))  %>% # this gives julian date (yday function)
  mutate(plantjul = as.Date(plantdate)) %>%
  mutate(plantjul = yday(plantjul))

# calculate time to germ and time since sept 1
germdaty3 = germdaty2 %>%
  mutate(days2germ = difftime(germdate, plantdate)) %>% #difftime function calculates difference of days between 2 times
  mutate(plantSep = difftime(plantdate, ymd("2020-09-01"))) %>% #this calculates the difference in time from Sept 1 to the planting date
  mutate(germSep = difftime(germdate, ymd("2020-09-01"))) %>%  #this calculates the difference in time from Sept 1 to the germ date
  mutate(days2germ2 = as.numeric(germSep - plantSep)) %>% #checking that the numeric version works well.  It does!
  mutate(planted = 1) %>% #variable that allows us to count each individual in the data
  mutate(germinated = if_else(is.na(germdate) == TRUE, 0,1))

# calculate time to germ from when watering started again in Round 2: 09/15/2021
# add water date column to df

germdaty3$water.date=as.Date("2021-09-15")

germdaty4 = germdaty3 %>%
  mutate(RegermSep = difftime(germdate, ymd("2021-09-15"))) %>%  #this calculates the difference in time from rewater date to the germ date
  mutate(Redays2germ = as.numeric(RegermSep)) # makes the number of days numeric

#### Remove cones that actually germinated in Round 1 ####

# See: round.2.removals in ./Raw.Data for complete list of removals

germdaty5=germdaty4[-c(1178,1508,2012,663,1442,2372,541,2017,626,132,1576,17,1294,568,1528,1709),]

#### Final Round 2 germination phenology data frame #### 
# can be found in Formatted.Data

write.csv(germdaty5, file="germ.pheno.round.2.formatted.csv")






  