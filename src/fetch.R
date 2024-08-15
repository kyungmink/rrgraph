#! /usr/bin/Rscript
source("src/functions.R")
fetch <- poolSick(T20Files, T40Files, 2002:2013)
save(fetch, file = "met/fetch.RData")
source("src/functions.R")
