FROM r-base
WORKDIR /home/rrgraph
COPY install-igraph.R .
RUN R CMD BATCH --no-restore --no-save --no-timing install-igraph.R
COPY makehelp.R  .
CMD R CMD BATCH --no-restore --no-save --no-timing makehelp.R
