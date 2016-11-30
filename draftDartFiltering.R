#dart filtering. Perhaps in munge stage?

index.repro <- gl@other$metrics[,"RepAvg"] > 0.98
index.callrate <-  gl@other$metrics[,"CallRate"] > 0.90
index.coverage <- (gl@other$metrics[,"AvgCountRef"]+gl@other$metrics[,"AvgCountSnp"] ) > 20 
index.highhet <- gl@other$metrics[,"FreqHets"] <0.75
index.comb <- filter.dart(index.repro, index.callrate, index.coverage)

fffff<-gl.fixed.diff(gi, t=0.05)

