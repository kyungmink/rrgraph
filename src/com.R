#The community information class, that works like the Ca class 
#because similar information has to be drawn.
#Dependency
#functions.R
Co <- function(ig = NULL) {
    if(is(ig, "igraph")) {
        #The community information table should contain
        #id (integer)
        #col (color, hex)
        clu <- cluster_walktrap(ig)
        num <- sizes(clu)
        num <- sort(num, decreasing = TRUE)
        l <- length(num)
        id <- 1:l
        g <- length(which(num > sum(num) %/% 10))
        col <- substr(c(rainbow(g), rainbow(l - g, 0.4, 0.8)), 1, 7)
        r <- data.frame(id, col, num, stringsAsFactors = FALSE)
        m <- modularity(ig, membership(clu))
        return(structure(r, class = c("Co", "data.frame"), modularity = m, communities = clu))
    } else {
        return(read.csv("com.csv", stringsAsFactors = FALSE))
    }
}
write.Co <- function(...) UseMethod("write.Ca", file = "com.csv")
co <- Co(igs)
if(exists("n")) n <- updat(n, co)
