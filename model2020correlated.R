library(truncdist)
rinvgammatrun2 = function(n, a, b, trun){
  for(i in 1:n){
    count = 1
    samples = 1/rtrunc(n, "gamma", a = 1/trun, b = Inf, shape = a, scale = 1/b)
  }
  
  return(samples)
}


propnorm = function(data){
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("data is not matrix or data frame")
  }
  normdata = data / rowSums(data)
  
  return(normdata)
}


propnormdf = function(df, columns){
  data = df[columns]
  normdata = propnorm(data)
  colnames(normdata) = paste(colnames(data), "normalized")
  return(cbind(df,normdata))
}


propnormdfreplace = function(df, columns){
  normdata = propnorm(df[, columns])
  df[, columns] = normdata
  return(df)
}  

modelstatesmvn = function(poll_data, stateloc, proploc, candidateloc, idloc, varloc = NULL, nloc, npolls = NULL, dateloc = NULL, statepriormeans, statepriorcov, tauk = rep(1, nrow(poll_data)), a0i = 1.01, b0i = 1.01, trunc = .125, nMCMC = 11000, burnin = 1000){
  #normalizes the polling data
  newdf = propnormdfreplace(poll_data, proploc)
  #gets the list of state names
  statenum = length(statepriormeans[,1])
  statenames = statepriormeans[,1]
  
  #implementing last n polls
  #Step 1: Loop over states
  if(!is.null(npolls)){
    lastpolls = numeric(0)
    for(i in 1:length(statenames)){
      #Step 2: For each state get index of last n polls and add to array of last polls
      poll_temp = subset(newdf, newdf[, stateloc] == statenames[i])
      poll_temp = poll_temp[order(poll_temp[, dateloc], decreasing = F), ]
      if(nrow(poll_temp) > npolls){
        lastpolls = c(lastpolls, poll_temp[1:10, idloc])
      }
      else{
        lastpolls = c(lastpolls, poll_temp[, idloc])
      }
      
      #Step 3: redefine poll data
      #go through original 
      
    }
    tauk = tauk[newdf[, idloc] %in% lastpolls]
    newdf = newdf[newdf[, idloc] %in% lastpolls, ]
    
  }
  finaldf = newdf
  
  #initialize posterior mean and variance
  statemeansMCMC = matrix(NA, nrow = statenum, ncol = nMCMC) 
  statevarMCMC = matrix(NA, nrow = statenum, ncol = nMCMC)
  
  winprob = rep(NA, statenum)
  npollstate = rep(NA, statenum)
  statecovinverse = solve(statepriorcov[,-1])
  #statepriorcov = statepriorcov[colnames(statepriorcov) %in% statenames, (colnames(statepriorcov) %in% statenames)[-1]]
  priorstatemeans = statepriormeans[,2]
  #statecovinverse = statecovinverse[rownames(statecovinverse) %in% statenames, rownames(statecovinverse) %in% statenames]
  muitemp = priorstatemeans
  sigmaitemp = as.numeric(diag(as.matrix(statepriorcov)[, -1]))
  #invert sigmaC
  
  for(i in 1:statenum){
    #get number of polls in state i
    npollstate[i] = sum(finaldf[, stateloc] == statenames[i])
    tauk[finaldf[,stateloc] == statenames[i]] = npollstate[i]*tauk[finaldf[,stateloc] == statenames[i]]/sum(tauk[finaldf[,stateloc] == statenames[i]])
    #tauk = nrow(poll_temp)*poll_temp[,weightloc]/sum(poll_temp[, weightloc])
    
  }
  stateswithnopolls = statenames[which(npollstate == 0, arr.ind = T)]
  stateswithpolls = statenames[which(npollstate > 0, arr.ind = T)]
  for(k in 1:nMCMC){
    for(i in 1:statenum){
      if(npollstate[i] == 0){
        sumfornormal = sum(statecovinverse[i,-i]*(muitemp[-i]-priorstatemeans[-i]))
        normmeani = (-sumfornormal+statecovinverse[i,i]*priorstatemeans[i])*(statecovinverse[i,i])^(-1) 
        normvari = (statecovinverse[i,i])^(-1)
        muitemp[i] = rnorm(1, mean = normmeani, sd = sqrt(normvari))
        
        statemeansMCMC[i,k] = muitemp[i]
      }
      else{
        
        #normmeani = ((1/sigmajtemp[groupnameloc] + (1/sigmaitemp[i]*sum(tauk[poll_data[,stateloc] == statenames[i]])))^-1)*((1/sigmaitemp[i])*(sum(poll_data[poll_data[,stateloc] == statenames[i],candidateloc]/tauk[poll_data[,stateloc] == statenames[i]]))+mujtemp[groupnameloc]/sigmajtemp[groupnameloc])
        iga = a0i + npollstate[i]/2
        igb = b0i + 0.5*sum((finaldf[finaldf[,stateloc] == statenames[i],candidateloc]- muitemp[i])^2/tauk[finaldf[,stateloc] == statenames[i]])
        #print(iga)
        #print(igb)
        #print(i)
        sigmaitemp[i] = rinvgammatrun2(1, iga, igb, trunc)
        
        statevarMCMC[i,k] = sigmaitemp[i]
        sumfornormal = sum(statecovinverse[i,-i]*(muitemp[-i]-priorstatemeans[-i]))
        normmeani= ((statecovinverse[i,i]+sum(1/(sigmaitemp[i]*tauk[finaldf[,stateloc] == statenames[i]])))^-1)*((sum(finaldf[finaldf[,stateloc] == statenames[i],candidateloc]/(sigmaitemp[i]*tauk[finaldf[,stateloc] == statenames[i]])))-sumfornormal+statecovinverse[i,i]*priorstatemeans[i])
        normvari = (statecovinverse[i,i]+sum(1/(sigmaitemp[i]*tauk[finaldf[,stateloc] == statenames[i]])))^-1
        muitemp[i] = rnorm(1, mean = normmeani, sd = sqrt(normvari))
        statemeansMCMC[i,k] = muitemp[i]
      }
    }
    if(k %% 100 == 0){
      print(k)
    }
  }
  poststatemeans = rowMeans(statemeansMCMC[,(burnin+1):nMCMC])
  poststatevar = rowMeans((statemeansMCMC[,(burnin+1):nMCMC]-poststatemeans)^2)
  poststatesd = sqrt(poststatevar)
  statewinmat = statemeansMCMC[,(burnin+1):nMCMC] > .50
  evvotemat = statewinmat*electoralvotes$Electoral.Votes
  winprob = mean(colSums(((statemeansMCMC[,(burnin+1):nMCMC] > .50)*electoralvotes$Electoral.Votes))>270)
  statewinprobability = rowMeans(statemeansMCMC[,(burnin+1):nMCMC] > .50)
  cilow = apply(statemeansMCMC[, (burnin+1):nMCMC], 1, quantile, prob = 0.025)
  cihigh = apply(statemeansMCMC[, (burnin+1):nMCMC], 1, quantile, prob = 0.975)
  
  df = data.frame(state = statenames, means = poststatemeans, poststatevar = poststatevar, poststatesd = poststatesd, cilow = cilow, cihigh=cihigh, statewinprobability = statewinprobability)
  
  return(list(df = df, statenames = statenames, poststatemeans = poststatemeans, poststatevar = poststatevar, winprob = winprob, statemeansMCMC = statemeansMCMC, statevarMCMC = statevarMCMC, statewinmat = statewinmat, evvotemat = evvotemat))
  
}
source("dataprocessing.R")
#drop tracking polls
dataeconnew = dataecon[dataecon$pollster != "Morning Consult"
                       & dataecon$pollster != "USC Dornsife" & dataecon$pollster != "Trafalgar Group", ]
dataeconnew = dataeconnew[dataeconnew$MidDaysUntil < 100, ]
dataeconnew = dataeconnew[dataeconnew$state != "--",]
dataeconnew$id = 1:nrow(dataeconnew)
`%notin%` <- Negate(`%in%`)
dataeconnew$state = as.character(dataeconnew$state)

dataeconnew$competitivefrompoll = ifelse(abs(dataeconnew$biden_margin) < 5, "competitive", "noncompetitive")
dataeconnew$MidDaysFactor = cut(as.numeric(dataeconnew$MidDaysUntil), breaks = c(seq(0, 84, by = 7), 150, 300, 1000), labels = c("Last week",  paste(2:12, rep("Weeks to Go", 11)), "Summer", "Primary", "Pre-Primary"))
stateswithpolls2020  = unique(dataeconnew$state)
library(dplyr)
dataeconnew$mode = recode_factor(dataeconnew$mode, 'IVR' = "Automated Phone", 'SurveyUSA' = "Internet", 'Online' = "Internet", 'Text/Online' = "Internet", 'Live Phone/Text' = "Live Phone/Online", 'Live Pone' = "Live Phone")
dataeconnew$sample_subpopulation = recode_factor(dataeconnew$population, 'a' = "A", 'lv' = "LV", 'lv '= "LV", 'rv ' = "RV", 'rv' = "RV")
dataeconnew= dataeconnew[dataeconnew$sample_subpopulation != "A", ]
library(politicaldata)
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
elect2016 = subset(pres_results , year == 2016)
elect2016$margin = elect2016$dem - elect2016$rep
elect2016$Dem2P = elect2016$dem/(elect2016$dem + elect2016$rep)
elect2016$Rep2P = elect2016$rep/(elect2016$dem + elect2016$rep)
dataeconnew$PVRepVote = rep(NA, nrow(dataeconnew))
dataeconnew$PVDemVote = rep(NA, nrow(dataeconnew))
dataeconnew$PVMarginVote = rep(NA, nrow(dataeconnew))
library(MCMCpack)
library(truncdist)
electdata = data.frame("state" = elect2008$state, "1976" = elect1976$margin, "1980" = elect1980$margin, "1984" = elect1984$margin, "1988" = elect1988$margin, "1992" =elect1992$margin, "1996" =elect1996$margin, "2000" = elect2000$margin, "2004" = elect2004$margin, "2008" =  elect2008$margin, "2012" = elect2012$margin, "2016" = elect2016$margin)
for(i in 1:length(stateswithpolls2020)){
  electtemp = elect2016[which(elect2016$state == stateswithpolls2020[i]),]
  dataeconnew$PVRepVote[which(dataeconnew$state == stateswithpolls2020[i])] = electtemp$Rep2P
  dataeconnew$PVDemVote[which(dataeconnew$state == stateswithpolls2020[i])] = electtemp$Dem2P
  dataeconnew$PVMargin[which(dataeconnew$state == stateswithpolls2020[i])] = electtemp$margin
}
dataeconnew$Dem2P = dataeconnew$biden/(dataeconnew$biden+dataeconnew$trump)
dataeconnew$Rep2P = dataeconnew$trump/(dataeconnew$biden+dataeconnew$trump)
dataeconnew = dataeconnew[complete.cases(dataeconnew[, c("biden", "trump", "mode")]), ]
dataeconnew$weight = predict(modelfor2020, newdata= dataeconnew)
dataeconnew$mode = recode_factor(dataeconnew$mode,'Onilne' = "Internet")
source("getweights.R")
weights2020 = predict(modelfor2020, newdata= dataeconnew)

electdata = electdata[,c(1, 7,8)]
electdata$state = as.character(electdata$state)
priorstatemeans2016 = elect2016[, c("state", "dem")]
library(rockchalk)

fitoverall = modelstatesmvn(dataeconnew, stateloc = 1,  proploc = c(29:30), candidateloc = 29, nloc = 6,  statepriormeans = priorstatemeans2016, statepriorcov = priorstatecov, tauk = weights2020, trunc  = 0.0625, nMCMC = 202500, burnin = 2500)

fit = fitoverall$df
getpriorassign(electdata)
npollsperstate = rep(NA, nrow(electdata))
for(i in 1:nrow(electdata)){
  npollsperstate[i] = nrow(dataeconnew[fit$state[i] == dataeconnew$state,])
  
}
fit$npollsperstate = npollsperstate
#fit$groupassignment = getpriorassign(electdata, cutoffs = c(-.25, -.125,  -0.075, 0.075, .125, .25))[,2]
fit$margin = fit$means - (1-fit$means)
filename = "lastmodelrun2020.csv"
library(googledrive)
drive_upload(filename, path = as_dribble("/ElectionModel"))
