#clean RScript
library(bayesurvey)
library(politicaldata)
library(MCMCpack)
elect1992  = subset(pres_results , year == 1992)
elect1992$margin = elect1992$dem - elect1992$rep
elect1996  = subset(pres_results , year == 1996)
elect1996$margin = elect1996$dem - elect1996$rep
elect2000  = subset(pres_results , year == 2000)
elect2000$margin = elect2000$dem - elect2000$rep
elect2004 = subset(pres_results , year == 2004)
elect2004$margin = elect2004$dem - elect2004$rep
elect2008  = subset(pres_results , year == 2008)
elect2008$margin = elect2008$dem - elect2008$rep
elect2012 = subset(pres_results , year == 2012)
elect2012$margin = elect2012$dem - elect2012$rep
elect2016 = subset(pres_results , year == 2016)
elect2016$margin = elect2016$dem - elect2016$rep
electdata = data.frame("state" = elect2008$state, "1992" =elect1992$margin, "1996" =elect1996$margin, "2000" = elect2000$margin, "2004" = elect2004$margin, "2008" =  elect2008$margin, "2012" = elect2012$margin, "2016" = elect2016$margin)
source("dataprocessing.R")
#drop tracking polls
dataeconnew = dataecon[dataecon$pollster != "Morning Consult"
                       & dataecon$pollster != "USC Dornsife", ]
dataeconnew = dataeconnew[dataeconnew$MidDaysUntil < 100, ]
dataeconnew = dataeconnew[dataeconnew$state != "--",]
dataeconnew$id = 1:nrow(dataeconnew)
`%notin%` <- Negate(`%in%`)
dataeconnew$state = as.character(dataeconnew$state)
library(bayesurvey)
library(MCMCpack)
electdata = electdata[,c(1, 7,8)]
electdata$state = as.character(electdata$state)
fit = noniterativegaussianmodelprop(dataeconnew, stateloc = 1, proploc = c(10,11), 
                                    candidateloc = 10, nloc = 7, dateloc = 5, idloc = 21,  
                                    npolls = 10,invgamma = T, v0=1, a0=0.0001, b0=0.0001,  
                                    election_data = electdata, 
                                    cutoffs = c(-.25, -.125, -0.075, 0.075, .125, .25))
fit
getpriorassign(electdata)
npollsperstate = rep(NA, nrow(electdata))
for(i in 1:nrow(electdata)){
  npollsperstate[i] = nrow(dataeconnew[fit$State[i] == dataeconnew$state,])
  
}
fit$npollsperstate = npollsperstate
fit$groupassignment = getpriorassign(electdata, cutoffs = c(-.25, -.125,  -0.075, 0.075, .125, .25))[,2]
fit = fit[, - 8]
filename = paste0("fit", Sys.Date(), ".csv")
write.csv(fit, file = filename )
library(googledrive)
drive_upload(filename, path = as_dribble("/ElectionModel"))
