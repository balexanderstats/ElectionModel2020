#download and process data
#get covariance matrix
#this code is taken from the economist potus model github
#license info
# 
# MIT License
# 
# Copyright (c) 2020 The Economist Newspaper, Andrew Gelman and Merlin Heidemanns (Columbia University)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

library(RCurl)
dataecon = read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ56fySJKLL18Lipu1_i3ID9JE06voJEz2EXm6JW4Vh11zmndyTwejMavuNntzIWLY0RyhA1UsVEen0/pub?gid=0&single=true&output=csv")


dataecon$end.date = as.Date(dataecon$end.date, tryFormats= ("%m/%d/%y"))
dataecon$start.date = as.Date(dataecon$start.date, tryFormats= ("%m/%d/%y"))
dataecon$EndDaysUntil = difftime(as.Date("2020-11-3"),dataecon$end.date, units = "days")
dataecon$StartDaysUntil = difftime(as.Date("2020-11-3"),dataecon$start.date, units = "days") 
dataecon$MidDaysUntil = (dataecon$EndDaysUntil+dataecon$StartDaysUntil)/2
  