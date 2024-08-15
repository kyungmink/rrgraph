#The Ca class represents a KCD category information table
#The ca object is the default KCD category information table
#It is stored as met/cat.csv
#This script writes met/cat.csv and met/nodes.csv.
#Dependencies:
#functions.R
Ca <- function(new = FALSE) {

#The category color reference table should contain
#name (roman numerals)
#kcd (range of kcd values)
#desc (decription)
#col (color, hex)

    #size, transitivity.

    if (new) {

#Use kcd.csv, the KCD codes sheet saved as utf8 csv.

    #Use kcd.csv, the KCD codes sheet saved as utf8 csv.
    kcd <- read.csv("raw/kcd.csv", stringsAsFactors = FALSE)
    s <- subset(kcd, 분류.기준 == "대")[c("질병분류.코드", "한글명칭", "영문명칭")]
    name <- as.character(as.roman(1:dim(s)[[1]]))
    e <- regexpr("\\(", s$영문명칭)
    desc <- substr(s$영문명칭, 1, e - 1)
    desc[20] <- substr(desc[20], 4, nchar(desc[20]))
    col = substr(rainbow(length(name)), 1, 7)

        ca <- list(name = name, kcd = s$질병분류.코드, desc = desc, col = col)
        ca <- data.frame(ca, stringsAsFactors = FALSE)

    } else {
        ca <- read.csv("met/cat.csv", stringsAsFactors = FALSE)
    }
    return(structure(ca, class = c("Ca", "data.frame")))
}
write.Ca <- function(ca, file = "met/cat.csv", ...) {
    #Store the table as cat.csv
    write.csv(ca, file, row.names = FALSE, ...)
}
updat.Ca <- function(ca, x) {
    if(is(x, "N")) {
        s <- sizes(x$cat)
        ca$num <- ca$nui <- rep(0, dim(ca)[[1]])
        ca$num[as.integer(names(s))] <- s
        s <- sizes(subset(x, !is.na(com))$cat)
        ca$nui[as.integer(names(s))] <- s
        return(ca)
    } else if (is(x, "igraph")) {
        m <- modularity(x, as.integer(n$cat[match(V(x)$name, n$name)]), weights = E(x)$weight)
        attr(ca, "modularity") <- m
        return(ca)
    } else stop()
}
ca <- Ca(FALSE)
if(exists("n") && is(n, "N")) {
    n <- updat(n, ca)
    ca <- updat(ca, n)
    write(n)
}
write(ca)
