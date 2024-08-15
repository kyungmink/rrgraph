#!/usr/bin/Rscript
load('rrgraph/data/incidences.rda')
incidences <- fetch
save(incidences, file = 'rrgraph/data/incidences2.rda', compress = 'xz')
