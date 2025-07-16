
setwd("~/IDA Lab/Dissertation/Paper4-Cooperations/Codes_R/github")
source("mpcs.R")

set.seed(235)
csv_names <- c("schizophrenia.csv","schizophreniform_disorder.csv")

names <- stringr::str_remove(csv_names,".csv")
disorder_list <- vector(mode="list",length=length(csv_names))

for(i in 1:length(csv_names)){
  disorder_list[[i]] <- read.csv(paste0("./disorders/",csv_names[i]), header=T)
}

df_all_disorder <- plyr::join_all(disorder_list,type="full",match="first")
df_all_disorder[is.na(df_all_disorder)] <- 0
d1_len <- dim(disorder_list[[1]])[1]
d2_len <- dim(disorder_list[[2]])[1]
t1 <- as.matrix(df_all_disorder[1:d1_len,-1])
t2 <- as.matrix(df_all_disorder[(d1_len+1):(d1_len+d2_len),-1])

dis <- MPCS(t1, t2, "mean")
t1t2 <- rbind(t1, t2)
dists <- distCosSim(t1t2) #cosine/angular distance
# if points overlap (dist=0) sammon mapping is not possible
dists[which(dists==0)] <- dists[which(dists==0)] + 0.0001

# Sammon Mapping
test.sam <- MASS::sammon(dists, y=jitter(cmdscale(dists, 2)) )
par(mar=c(8.1, 4.1, 4.1, 2.1))
plot(test.sam$points, type="n", main="Sammon Mapping", xlab="", ylab="", bty="n", asp=1, xlim=c(-0.5,0.5),ylim=c(-0.8,0.8))#, xlim=c(-0.6,0.6), ylim=c(-0.8,0.8))
abline(h=seq(-10,10,by=0.1), v=seq(-10,10,by=0.1), col="gray", lty=2)
abline(h=0, v=0, lwd=1.5, col="black")
points(matrix(test.sam$points[1:nrow(t1),],ncol=2), pch=20, col="red", cex=3)
points(matrix(test.sam$points[(nrow(t1)+1):nrow(test.sam$points),],ncol=2), pch=20, col="blue", cex=3)
mtext(paste0("MPCS: ",round(dis,2)),side=3, cex=1.2)
#legend(0.1,0.9, names, col=c("red","blue"), pch=20, cex=2, bty="n")
legend("bottom", inset=-0.20, names, col=c("red","blue"), horiz=T, pch=20, pt.cex=3, cex=1.5, xpd=T, text.width=0.4, bty='n')

# save sammon points as data
#df_clusters <- data.frame(x=test.sam$points[,1], y=test.sam$points[,2], disorder=c(rep("schizophrenia",nrow(t1)), rep("schizophreniform",nrow(t2))))
#write.csv(df_clusters, "sammon_disorder_clusters.csv")

########## other dimension reduction algorithms ##########

# Classic MDS
#test.cmds <- cmdscale(dists, eig=T, k=2)
#plot(test.cmds$points, type="n", main="Classic MDS", xlab="", ylab="", bty="n", asp=1, xlim=c(-0.25,0.65), ylim=c(-0.25,0.25)) #xlim=c(-0.4,0.4), )
#abline(h=seq(-10,10,by=0.1), v=seq(-10,10,by=0.1), col="gray", lty=2)
#abline(h=0, v=0, lwd=1.5, col="black")
#points(matrix(test.cmds$points[1:nrow(t1),],ncol=2), pch=20, col="red", cex=2)
#points(matrix(test.cmds$points[(nrow(t1)+1):nrow(test.cmds$points),],ncol=2), pch=20, col="blue", cex=2)
#mtext(paste0("Cosine Similarity: ",round(dis,2)),side=3, cex=1.2)
#legend(0.05,-0.05, names, col=c("red","blue"), pch=20, cex=2, bty="n")
#legend("bottom", inset=-0.2, names, col=c("red","blue"), horiz=T, pch=20, pt.cex=3, cex=1.5, xpd=T, text.width=0.2, bty='n')

# Nonmetric MDS
#test.imds <- MASS::isoMDS(dists, k=2)
#plot(test.imds$points, type="n", main="Nonmetric MDS", xlab="", ylab="", bty="n")#, xlim=c(-0.4,0.4), ylim=c(-0.3,0.3))
#abline(h=seq(-10,10,by=0.1), v=seq(-10,10,by=0.1), col="gray", lty=2)
#abline(h=0, v=0, lwd=1.5, col="black")
#points(matrix(test.imds$points[1:nrow(t1),],ncol=2), pch=20, col="red", cex=2)
#points(matrix(test.imds$points[(nrow(t1)+1):nrow(test.imds$points),],ncol=2), pch=20, col="blue", cex=2)
#mtext(paste0("Cosine Similarity: ",round(dis,2)),side=3)
#legend(-0.1,0, names, col=c("red","blue"), pch=20, cex=2, bty="n")
