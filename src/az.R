kPaste <- function(x, y) {
    paste(x, y, sep = "")
}
kCodes <- as.vector(outer(outer(toupper(letters), 0:9, kPaste), 0:9, kPaste))
