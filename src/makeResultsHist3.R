datCom <- read.csv("pValue.csv")
library(igraph)
graSta <- graph_from_data_frame(subset(datCom, RR > 4 & pFdr < 0.001))
E(graSta)$weight <- E(graSta)$RR
clu <- cluster_walktrap(graSta)
datNod <- read.delim("nodes.cys.tsv", quote = "")
datNod$class <- substr(datNod$kcd, 1, 1)
datNod$community <- membership(clu)[match(datNod$kcd, names(membership(clu)))]
write.table(datNod, file = "nodes.cys.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(datCom, RR > 4 & pFdr < 0.001), file = "edges.cys.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
is.hierarchical(clu)
plot_dendrogram(clu)
datNod$prevalence <- datNod$number / 1014730 * 100000
datNod$sexRatio <- ifelse(datNod$sex1 > datNod$sex2, paste(datNod$sex1 / datNod$sex2, ": 1"), paste("1 :", datNod$sex2 / datNod$sex1))
datNod$meanAge <- datNod$ageMean * 5 - 3
length(which(is.na(datNod$ageMean)))
datNod[which(is.na(datNod$ageMean)), c("kcd", "eng")]
datNod$prevalence <- round(datNod$prevalence)
datNod$sexRatio <- ifelse(datNod$sex1 > datNod$sex2, paste(round(datNod$sex1 / datNod$sex2, 2), ": 1"), paste("1 :", round(datNod$sex2 / datNod$sex1, 2)))
datNod$meanAge <- round(datNod$meanAge, 1)
com2 <- subset(datNod, community == 2)[,c("kcd", "eng", "prevalence", "sexRatio", "meanAge")]
head(com2)
com2 <- com2[sort(com2$prevalence, decreasing = TRUE, index.return = TRUE)$ix,]
com5 <- subset(datNod, community == 5)[,c("kcd", "eng", "prevalence", "sexRatio", "meanAge")]
com4 <- subset(datNod, community == 4)[,c("kcd", "eng", "prevalence", "sexRatio", "meanAge")]
com1 <- subset(datNod, community == 1)[,c("kcd", "eng", "prevalence", "sexRatio", "meanAge")]
write.csv(com2, file = "community2.csv", row.names = FALSE)
com5 <- com5[sort(com5$prevalence, decreasing = TRUE, index.return = TRUE)$ix,]
com4 <- com4[sort(com4$prevalence, decreasing = TRUE, index.return = TRUE)$ix,]
com1 <- com1[sort(com1$prevalence, decreasing = TRUE, index.return = TRUE)$ix,]
write.csv(com5, file = "community5.csv", row.names = FALSE)
write.csv(com4, file = "community4.csv", row.names = FALSE)
write.csv(com1, file = "community1.csv", row.names = FALSE)
datNod$betweenness <- betweenness(graSta)[match(datNod$kcd, names(betweenness(graSta)))]
datNod$class <- factor(datNod$class)
datNodFac <- scan("kcd-labels.txt", what = "", sep = "\n")
levels(datNod$class) <- datNodFac
plot(datNod$class, log(datNod$betweenness), xlab = "log(betweenness)", horizontal = TRUE)
plot(datNod$class, log(datNod$betweenness), xlab = "log(betweenness)", horizontal = TRUE, las = 1)
plot(datNod$class, log(datNod$betweenness), xlab = "log(betweenness)", horizontal = TRUE, las = 1, mar = mar()$mar +c(0, 15, 0, 0))
plot(datNod$class, log(datNod$betweenness), xlab = "log(betweenness)", horizontal = TRUE, las = 1, mar = par()$mar +c(0, 15, 0, 0))
tiff("betw-vs-class.tiff")
plot(datNod$class, log(datNod$betweenness), xlab = "log(betweenness)", horizontal = TRUE, las = 1, mar = par()$mar +c(0, 15, 0, 0))
dev.off()
tiff("betw-vs-class.tiff")
par(mar = par()$mar + c(0, 15, 0, 0))
plot(datNod$class, log(datNod$betweenness), xlab = "log(betweenness)", horizontal = TRUE, las = 1)
dev.off()
tiff("betw2-vs-class.tiff")
par(mar = par()$mar + c(0, 15, 0, 0))
plot(datNod$class, datNod$betweenness, xlab = "betweenness", horizontal = TRUE, las = 1)
dev.off()
E(graCom)$weight <- E(graCom)$RR
graCom <- graCom - E(graCom)[is.na(E(graCom)$weight)]
datNod$outs <- strength(graCom, mode = "out")[as.character(datNod$kcd)]
datNod$ins <- strength(graCom, mode = "in")[as.character(datNod$kcd)]
datNod[head(sort(datNod$outs + 0.65 * datNod$ins, decreasing = TRUE, index.return = TRUE)$ix),]
ruq <- datNod[head(sort(datNod$outs + 0.65 * datNod$ins, decreasing = TRUE, index.return = TRUE)$ix),]
ruq <- ruq[,c("kcd", "eng", "community", "prevalence", "sexRatio", "meanAge")]
write.csv(ruq, file = "ruq.csv", row.names = FALSE)
datNod[head(sort(datNod$outs + 0.65 * datNod$ins, index.return = TRUE)$ix),]
llq <- datNod[head(sort(datNod$outs + 0.65 * datNod$ins, index.return = TRUE)$ix),]
llq <- llq[,c("kcd", "eng", "community", "prevalence", "sexRatio", "meanAge")]
write.csv(llq, file = "llq.csv", row.names = FALSE)
datNod[head(sort(datNod$ins - datNod$outs * 0.65, decreasing = TRUE, index.return = TRUE)$ix),]
luq <- datNod[head(sort(datNod$ins - datNod$outs * 0.65, decreasing = TRUE, index.return = TRUE)$ix),]
luq <- luq[,c("kcd", "eng", "community", "prevalence", "sexRatio", "meanAge")]
write.csv(luq, file = "luq.csv", row.names = FALSE)
datNod[head(sort(datNod$ins - datNod$outs * 0.65, index.return = TRUE)$ix),]
rlq <- datNod[head(sort(datNod$ins - datNod$outs * 0.65, index.return = TRUE)$ix),]
rlq <- rlq[,c("kcd", "eng", "community", "prevalence", "sexRatio", "meanAge")]
write.csv(rlq, file = "rlq.csv", row.names = FALSE)
# Creating barplot of communities class composition
x <- table(datNod$class, datNod$community)[, c(2, 5, 4, 1)]
x <- matrix(x, nrow = 21)
barplot(x / rep(colSums(x), rep(21, 4)), col = rainbow(21))
