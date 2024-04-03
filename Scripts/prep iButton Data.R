# Code to prep ibutton data

#rm(list = ls())
library("tidyverse")
library("lubridate")
library("purrr")

##### reading in ibutton data for each block and adding cohort and block columns ####
#Cohort 1
#Block 1
C1B1 = read.csv('./Germination.Timing/Raw.Data/Cohort 1 Block 1 iButton Data.csv', skip = 14)
head(C1B1)
tail(C1B1)
#C1B1= C1B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C1B1$cohort <- 1
C1B1$block <- 1
C1B1$Date.Time <- mdy_hm(C1B1$Date.Time) #change this because original format didn't have seconds and the rest of the files do,
#even though format says hm, check the values, it gives seconds (as 00)

#Cohort 1
#Block 1 Part 2
C1B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 1 Block 1 (2) iButton Data.csv", skip = 14)
#C1B3= C1B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C1B1.2$cohort <- 1
C1B1.2$block <- 1
C1B1.2$Date.Time <- mdy_hms(C1B1.2$Date.Time)


#Cohort 1
#Block 2
C1B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 1 Block 2 iButton Data.csv", skip = 14)
head(C1B2)
#C1B2= C1B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C1B2$cohort <- 1
C1B2$block <- 2
C1B2$Date.Time <- mdy_hms(C1B2$Date.Time)

#Cohort 1
#Block 2 Part 2
C1B2.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 1 Block 2 (2) iButton Data.csv", skip = 14)
#C1B3= C1B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C1B2.2$cohort <- 1
C1B2.2$block <- 2
C1B2.2$Date.Time <- mdy_hms(C1B2.2$Date.Time)

#Cohort 1
#Block 3
C1B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 1 Block 3 iButton Data.csv", skip = 14)
#C1B3= C1B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C1B3$cohort <- 1
C1B3$block <- 3
C1B3$Date.Time <- mdy_hms(C1B3$Date.Time)

#Cohort 1
#Block 3 Part 2
C1B3.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 1 Block 3 (2) iButton Data.csv", skip = 14)
#C1B3= C1B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C1B3.2$cohort <- 1
C1B3.2$block <- 3
C1B3.2$Date.Time <- mdy_hms(C1B3.2$Date.Time)

#Cohort 2 # No dataframe in Raw.Data?
#Block 1
#C2B1 = read.csv("./Germination.Timing/Raw.Data/Cohort 2 Block 1 iButton Data.csv", skip = 14)
#C2B1= C2B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
#C2B1$cohort <- "2"
#C2B1$block <- "1"

#Cohort 2
#Block 1 Part 2
C2B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 2 Block 1 (2) iButton Data.csv", skip = 14)
#C2B1.2= C2B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C2B1.2$cohort <- 2
C2B1.2$block <- 1
C2B1.2$Date.Time <- mdy_hms(C2B1.2$Date.Time)

#Cohort 2
#Block 2
C2B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 2 Block 2 iButton Data.csv", skip = 14)
#C2B2= C2B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C2B2$cohort <- 2
C2B2$block <- 2
C2B2$Date.Time <- mdy_hms(C2B2$Date.Time)

#Cohort 2
#Block 2 Part 2
C2B2.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 2 Block 2 (2) iButton Data.csv", skip = 14)
#C2B2= C2B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C2B2.2$cohort <- 2
C2B2.2$block <- 2
C2B2.2$Date.Time <- mdy_hms(C2B2.2$Date.Time)

#Cohort 2
#Block 3
C2B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 2 Block 3 iButton Data.csv", skip = 14)
#C2B3= C2B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C2B3$cohort <- 2
C2B3$block <- 3
C2B3$Date.Time <- mdy_hm(C2B3$Date.Time)

#Cohort 3
#Block 1
C3B1 = read.csv("./Germination.Timing/Raw.Data/Cohort 3 Block 1 iButton Data.csv", skip = 14)
#C3B1= C3B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C3B1$cohort <- 3
C3B1$block <- 1
C3B1$Date.Time <- mdy_hms(C3B1$Date.Time)

#Cohort 3
#Block 1 Part 2
C3B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 3 Block 1 (2) iButton Data.csv", skip = 14)
#C3B1= C3B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C3B1.2$cohort <- 3
C3B1.2$block <- 1
C3B1.2$Date.Time <- mdy_hms(C3B1.2$Date.Time)

#Cohort 3
#Block 2
C3B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 3 Block 2 iButton Data.csv", skip = 14)
#C3B2= C3B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C3B2$cohort <- 3
C3B2$block <- 2
C3B2$Date.Time <- mdy_hms(C3B2$Date.Time)

#Cohort 3
#Block 2 Part 2
C3B2.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 3 Block 2 (2) iButton Data.csv", skip = 14)
#C3B2= C3B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C3B2.2$cohort <- 3
C3B2.2$block <- 2
C3B2.2$Date.Time <- mdy_hms(C3B2.2$Date.Time)

#Cohort 3
#Block 3
C3B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 3 Block 3 iButton Data.csv", skip = 14)
#C3B3= C3B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C3B3$cohort <- 3
C3B3$block <- 3
C3B3$Date.Time <- mdy_hms(C3B3$Date.Time)

#Cohort 3
#Block 3 Part 2
C3B3.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 3 Block 3 (2) iButton Data.csv", skip = 14)
#C3B3= C3B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C3B3.2$cohort <- 3
C3B3.2$block <- 3
C3B3.2$Date.Time <- mdy_hms(C3B3.2$Date.Time)

#Cohort 4
#block 1
C4B1 = read.csv("./Germination.Timing/Raw.Data/Cohort 4 Block 1 iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C4B1$cohort <- 4
C4B1$block <- 1
C4B1$Date.Time <- mdy_hms(C4B1$Date.Time)

#Cohort 4
#block 1 Part 2
C4B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 4 Block 1 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C4B1.2$cohort <- 4
C4B1.2$block <- 1
C4B1.2$Date.Time <- mdy_hms(C4B1.2$Date.Time)

#Cohort 4
#block 2
C4B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 4 Block 2 iButton Data.csv", skip = 14)
#C4B2= C4B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C4B2$cohort <- 4
C4B2$block <- 2
C4B2$Date.Time <- mdy_hms(C4B2$Date.Time)

#Cohort 4
#block 2 Part 2
C4B2.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 4 Block 2 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C4B2.2$cohort <- 4
C4B2.2$block <- 2
C4B2.2$Date.Time <- mdy_hms(C4B2.2$Date.Time)

#Cohort 4
#block 3
C4B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 4 Block 3 iButton Data.csv", skip = 14)
#C4B3= C4B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C4B3$cohort <- 4
C4B3$block <- 3
C4B3$Date.Time <- mdy_hms(C4B3$Date.Time)

#Cohort 4
#block 3 Part 2
C4B3.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 4 Block 3 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C4B3.2$cohort <- 4
C4B3.2$block <- 3
C4B3.2$Date.Time <- mdy_hms(C4B3.2$Date.Time)

#Cohort 5
#block 1
C5B1 = read.csv("./Germination.Timing/Raw.Data/Cohort 5 Block 1 iButton Data.csv", skip = 14)
#C5B1= C5B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C5B1$cohort <- 5
C5B1$block <- 1
C5B1$Date.Time <- mdy_hms(C5B1$Date.Time)

#Cohort 5
#block 1 Part 2
C5B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 5 Block 1 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C5B1.2$cohort <- 5
C5B1.2$block <- 1
C5B1.2$Date.Time <- mdy_hms(C5B1.2$Date.Time)

#Cohort 5
#block 2
C5B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 5 Block 2 iButton Data.csv", skip = 14)
#C5B2= C5B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C5B2$cohort <- 5
C5B2$block <- 2
C5B2$Date.Time <- mdy_hms(C5B2$Date.Time, truncated = 3) 
# Added truncated = 3, so it would parse correctly. Without this is changes 1/1/2021 to 1/1/2020 for some reason

#Cohort 5
#block 3
C5B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 5 Block 3 iButton Data.csv", skip = 14)
#C5B3= C5B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C5B3$cohort <- 5
C5B3$block <- 3
C5B3$Date.Time <- mdy_hms(C5B3$Date.Time)

#Cohort 5
#block 3 Part 2
C5B3.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 5 Block 3 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C5B3.2$cohort <- 5
C5B3.2$block <- 3
C5B3.2$Date.Time <- mdy_hms(C5B3.2$Date.Time)

#Cohort 6
#block 1
C6B1 = read.csv("./Germination.Timing/Raw.Data/Cohort 6 Block 1 iButton Data.csv", skip = 14)
#C6B1= C6B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C6B1$cohort <- 6
C6B1$block <- 1
C6B1$Date.Time <- mdy_hms(C6B1$Date.Time)

#Cohort 6
#block 1 Part 2
C6B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 6 Block 1 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C6B1.2$cohort <- 6
C6B1.2$block <- 1
C6B1.2$Date.Time <- mdy_hms(C6B1.2$Date.Time)

#Cohort 6
#block 2
C6B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 6 Block 2 iButton Data.csv", skip = 14)
#C6B2= C6B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C6B2$cohort <- 6
C6B2$block <- 2
C6B2$Date.Time <- mdy_hms(C6B2$Date.Time)

#Cohort 6
#block 2 Part 2
C6B2.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 6 Block 2 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C6B2.2$cohort <- 6
C6B2.2$block <- 2
C6B2.2$Date.Time <- mdy_hms(C6B2.2$Date.Time)

#Cohort 6
#block 3
C6B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 6 Block 3 iButton Data.csv", skip = 14)
#C6B3= C6B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C6B3$cohort <- 6
C6B3$block <- 3
C6B3$Date.Time <- mdy_hms(C6B3$Date.Time)

#Cohort 6
#block 3 Part 2
C6B3.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 6 Block 3 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C6B3.2$cohort <- 6
C6B3.2$block <- 3
C6B3.2$Date.Time <- mdy_hms(C6B3.2$Date.Time)

#Cohort 7
#block 1
C7B1 = read.csv("./Germination.Timing/Raw.Data/Cohort 7 Block 1 iButton Data.csv", skip = 14)
#C7B1= C7B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C7B1$cohort <- 7
C7B1$block <- 1
C7B1$Date.Time <- mdy_hms(C7B1$Date.Time)

#Cohort 7
#block 1 Part 2
C7B1.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 7 Block 1 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C7B1.2$cohort <- 7
C7B1.2$block <- 1
C7B1.2$Date.Time <- mdy_hms(C7B1.2$Date.Time)

#Cohort 7
#block 2
C7B2 = read.csv("./Germination.Timing/Raw.Data/Cohort 7 Block 2 iButton Data.csv", skip = 14)
#C7B2= C7B2 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C7B2$cohort <- 7
C7B2$block <- 2
C7B2$Date.Time <- mdy_hms(C7B2$Date.Time)

#Cohort 7
#block 2 Part 2
C7B2.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 7 Block 2 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C7B2.2$cohort <- 7
C7B2.2$block <- 2
C7B2.2$Date.Time <- mdy_hms(C7B2.2$Date.Time)


#Cohort 7
#block 3
C7B3 = read.csv("./Germination.Timing/Raw.Data/Cohort 7 Block 3 iButton Data.csv", skip = 14)
#C7B3= C7B3 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C7B3$cohort <- 7
C7B3$block <- 3
C7B3$Date.Time <- mdy_hms(C7B3$Date.Time)

#Cohort 7
#block 3 Part 2
C7B3.2 = read.csv("./Germination.Timing/Raw.Data/Cohort 7 Block 3 (2) iButton Data.csv", skip = 14)
#C4B1= C4B1 %>% separate("Date.Time", into = c("Date", "Time"), sep = " ", remove = TRUE)
C7B3.2$cohort <- 7
C7B3.2$block <- 3
C7B3.2$Date.Time <- mdy_hms(C7B3.2$Date.Time)

##### joining ibutton block datasets ####

ibuttonall=list(C1B1, C1B2, C1B3, C2B2, C2B3, C3B1, C3B2, C3B3, 
                C4B1, C4B2, C4B3, C5B1, C5B2, C5B3, C6B1, C6B2, C6B3, C7B1,
                C7B2, C7B3, C1B1.2, C1B2.2, C1B3.2, C2B1.2,C2B2.2, C3B1.2,
                C3B2.2, C3B3.2, C4B1.2, C4B2.2, C4B3.2, C5B1.2, C5B3.2, C6B1.2,
                C6B2.2, C6B3.2, C7B1.2, C7B2.2, C7B3.2) %>% 
  reduce(full_join)
dim(ibuttonall)
summary(ibuttonall)
#Jenny note 7/21/22 - changed code to keep date.time format, which will be useful for plotting, then will separate and format below
# formatting date to read in as dates instead of characters
#ibuttonall$Date <- mdy(ibuttonall$Date)
#ibuttonall$Date  <- as.Date(ibuttonall$Date)

#ibuttonall = rename(ibuttonall, temp = 'Value')

#Jenny adding dplyr/tidyverse way to calculate/reformat variables

ibuttonall = ibuttonall %>%
             rename( temp = 'Value')   %>% 
             mutate(hour = hour(Date.Time)) %>% #this gives just the hour
             mutate(minute = minute(Date.Time)) %>% #%>% 
            # mutate(Date = ymd(Date.Time)) %>%#%>% #example of how to do the above Date formatting in tidy form with pipes
             mutate(Date = as.Date(Date.Time))
summary(ibuttonall)

head(ibuttonall)
str(ibuttonall)
dim(ibuttonall)
summary(ibuttonall)

###### get max, min, mean daily temp for each cohort #####
dailytemp_coh = ibuttonall %>% 
  group_by(Date, cohort) %>% 
  summarize(max_dailytemp = max(temp), min_dailytemp = min(temp), mean_dailytemp = mean(temp))
  
  #samplesize = sum(planted, na.rm=T), %>% #sample size is the total # planted for each cohort and pop
  #mutate(se_mean_dailytemp = sd_dailytemp/sqrt(samplesize), %>%
  #mutate(cohortfact = as.factor(Cohort))

dailytemp_coh

###### write csv files ######
#writing full dataset to a CSV
write.csv(ibuttonall, file = "Germination.Timing/Formatted.Data/ibutton_allcohorts.csv", row.names = FALSE)

#writing daily temps to a CSV
write.csv(dailytemp_coh,  file = "Germination.Timing/Formatted.Data/ibutton_dailytemps.csv", row.names = FALSE)

### formatting ibutton data for round 2 of Germination Phenology ####
# reading in ibutton data for each block and adding cohort and block columns 
#Round 2
#Block 4
R2B4 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 04 iButton Data.csv', skip = 14)
head(R2B4)
tail(R2B4)
R2B4$block <- 4
R2B4$Date.Time <- mdy_hm(R2B4$Date.Time) #change this because original format didn't have seconds and the rest of the files do,
#even though format says hm, check the values, it gives seconds (as 00)

#Round 2
#Block 4 Part 2
R2B4.2 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 04 iButton Data (2).csv', skip = 14)
head(R2B4.2)
tail(R2B4.2)
R2B4.2$block <- 4
R2B4.2$Date.Time <- mdy_hms(R2B4.2$Date.Time, truncated = 3)

#Round 2
#Block 7
R2B7 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 07 iButton Data.csv', skip = 14)
head(R2B7)
tail(R2B7)
R2B7$block <- 7
R2B7$Date.Time <- mdy_hm(R2B7$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 7 Part 2
R2B7.2 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 07 iButton Data (2).csv', skip = 14)
head(R2B7.2)
tail(R2B7.2)
R2B7.2$block <- 7
R2B7.2$Date.Time <- mdy_hm(R2B7.2$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 8
R2B8 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 08 iButton Data.csv', skip = 14)
head(R2B8)
tail(R2B8)
R2B8$block <- 8
R2B8$Date.Time <- mdy_hm(R2B8$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 8 Part 2
R2B8.2 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 08 iButton Data (2).csv', skip = 14)
head(R2B8.2)
tail(R2B8.2)
R2B8.2$block <- 8
R2B8.2$Date.Time <- mdy_hm(R2B8.2$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 11
R2B11 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 11 iButton Data.csv', skip = 14)
head(R2B11)
tail(R2B11)
R2B11$block <- 11
R2B11$Date.Time <- mdy_hm(R2B11$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 11 Part 2
R2B11.2 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 11 iButton Data (2).csv', skip = 14)
head(R2B11.2)
tail(R2B11.2)
R2B11.2$block <- 11
R2B11.2$Date.Time <- mdy_hm(R2B11.2$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 13
R2B13 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 13 iButton Data.csv', skip = 14)
head(R2B13)
tail(R2B13)
R2B13$block <- 13
R2B13$Date.Time <- mdy_hm(R2B13$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

#Round 2
#Block 13 Part 2
R2B13.2 = read.csv('./Germination.Timing/Raw.Data/Germ Pheno 2.0 iButton Data/Round 2 Block 13 iButton Data (2).csv', skip = 14)
head(R2B13.2)
tail(R2B13.2)
R2B13.2$block <- 13
R2B13.2$Date.Time <- mdy_hm(R2B13.2$Date.Time) #change this because original format didn't have seconds and the rest of the files do,

##### joining ibutton block datasets ####

ibuttonall=list(R2B4,R2B4.2,R2B7,R2B7.2,R2B8,R2B8.2,R2B11,R2B11.2,R2B13,R2B13.2) %>% 
  reduce(full_join)
dim(ibuttonall)
summary(ibuttonall)

ibuttonall = ibuttonall %>%
  rename( temp = 'Value')   %>% 
  mutate(hour = hour(Date.Time)) %>% #this gives just the hour
  mutate(minute = minute(Date.Time)) %>% #%>% 
  # mutate(Date = ymd(Date.Time)) %>%#%>% #example of how to do the above Date formatting in tidy form with pipes
  mutate(Date = as.Date(Date.Time))

head(ibuttonall)
str(ibuttonall)
dim(ibuttonall)
summary(ibuttonall)


###### get max, min, mean daily temp #####
dailytemp_round_2 = ibuttonall %>% 
  group_by(Date) %>% 
  summarize(max_dailytemp = max(temp), min_dailytemp = min(temp), mean_dailytemp = mean(temp))

dailytemp_round_2

###### write csv files ######
#writing full dataset to a CSV
write.csv(ibuttonall, file = "Germination.Timing/Formatted.Data/ibutton_round_2.csv", row.names = FALSE)

#writing daily temps to a CSV
write.csv(dailytemp_round_2,  file = "Germination.Timing/Formatted.Data/ibutton_dailytemps_round_2.csv", row.names = FALSE)


