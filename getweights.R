
pollingdata = read.csv("PollingDataSetStateLevelOnly9-29-20.csv")
pollingdata$Undecided = as.numeric(pollingdata$Undecided)
library(rockchalk)
pollingdata$Undecided[pollingdata$Year == 2012 & is.na(pollingdata$Undecided)] =  100-pollingdata$Obama2012[pollingdata$Year == 2012 & is.na(pollingdata$Undecided)]-pollingdata$Romney[pollingdata$Year == 2012 & is.na(pollingdata$Undecided)]
pollingdata$Undecided[pollingdata$Year == 2008 & is.na(pollingdata$Undecided)] =  100-pollingdata$Obama2008[pollingdata$Year == 2008 & is.na(pollingdata$Undecided)]-pollingdata$McCain[pollingdata$Year == 2008 & is.na(pollingdata$Undecided)]
pollingdata$Undecided[pollingdata$Year == 2016 & is.na(pollingdata$Undecided)] =  100-pollingdata$Clinton[pollingdata$Year == 2016 & is.na(pollingdata$Undecided)]-pollingdata$Trump[pollingdata$Year == 2016 & is.na(pollingdata$Undecided)]
pollingdata$competitivefrompoll = ifelse(abs(pollingdata$PollMargin) < 5, "competitive", "noncompetitive")
pollingdata$sample_subpopulation = combineLevels(pollingdata$sample_subpopulation, c("LV", "Likely Voters"), "LV")
pollingdata$sample_subpopulation = combineLevels(pollingdata$sample_subpopulation, c("RV", "Registered Voters"), "RV")
pollingdata$sample_subpopulation = combineLevels(pollingdata$sample_subpopulation, c("A", "Adults"), "A")
pollingdata = pollingdata[pollingdata$sample_subpopulation != "A",]
pollingdata$sample_subpopulation = droplevels(pollingdata$sample_subpopulation)
pollingdata$pollminusPV = pollingdata$Dem2P - pollingdata$PVDemVote
polls2008 = subset(pollingdata, pollingdata$Year == 2008)
polls2012 = subset(pollingdata, pollingdata$Year == 2012)
polls2016 = subset(pollingdata, pollingdata$Year == 2016)
polldatatemp = subset(polldatatemp, polldata$EndDaysUntil < 100)
polls2008temp = subset(polldatatemp, Year == "2008")
polls2008temp = polls2008temp[complete.cases(polls2008temp[, c("Obama2008", "McCain", "observations","end_date")]),]
polls2012temp = subset(polldatatemp, Year == "2012")
polls2012temp = polls2012temp[complete.cases(polls2012temp[, c("Obama2012", "Romney", "observations", "end_date")]),]
polls2016temp = subset(polldatatemp, Year == "2016")
polls2016temp = polls2016temp[complete.cases(polls2016temp[, c("Trump", "Clinton", "observations", "end_date")]),]
pres_results = propnormdfreplace(pres_results, c(4,5))
pres_results$margin = pres_results$dem-pres_results$rep
label7 = c("Strong Red", "Red", "Lean Red", "Competitive",
           "Lean Blue", "Blue", "Strong Blue")
elect1976 = subset(pres_results , year == 1976)
elect1980 = subset(pres_results , year == 1980)
elect1984 = subset(pres_results , year == 1984)
elect1988 = subset(pres_results , year == 1988) 
elect1992  = subset(pres_results , year == 1992)
elect1996  = subset(pres_results , year == 1996)
elect2000  = subset(pres_results , year == 2000)
elect2004  = subset(pres_results , year == 2004)
elect2008  = subset(pres_results , year == 2008)
elect2012  = subset(pres_results , year == 2012)
elect2016  = subset(pres_results , year == 2016)
electdata = data.frame("state" = elect2008$state, "1976" = elect1976$margin, "1980" = elect1980$margin, "1984" = elect1984$margin, "1988" = elect1988$margin, "1992" =elect1992$margin, "1996" =elect1996$margin, "2000" = elect2000$margin, "2004" = elect2004$margin, "2008" =  elect2008$margin, "2012" = elect2012$margin, "2016" = elect2016$margin)
polls2008 = subset(pollingdata, pollingdata$Year == 2008)
polls2012 = subset(pollingdata, pollingdata$Year == 2012)
polls2016 = subset(pollingdata, pollingdata$Year == 2016)
polls_2012_2016 = rbind(polls2012, polls2016)
modelfor2020 =stan_glm(as.formula((VoteErrorAbs+0.001) ~  Dem2P + mode + MidDaysFactor + PVDemVote + competitivefrompoll + sample_subpopulation), family = Gamma(link = inverse), data = polls_2012_2016)