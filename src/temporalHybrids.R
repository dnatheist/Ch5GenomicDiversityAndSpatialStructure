
hyb<-data.frame(subset(larv$Species, grepl("^F", larv$Species), Day.of.Year))

hyb<-larv$Species
isHyb<-data.frame(subset(larv, grepl("^F", hyb), c(Day.of.Year, estimatedAge)))
notHyb<-data.frame(subset(larv, !grepl ("^F", hyb), c(Day.of.Year, estimatedAge)))

t.test(isHyb$Day.of.Year,notHyb$Day.of.Year)
t.test(isHyb$estimatedAge,notHyb$estimatedAge)

