## Here we import a propiertary dataset (processed data files not
## yet publicly available)

x <- rbind(read.csv("../../../iml2020/T/npf_train_spring2020.csv"),
	   read.csv("../../../iml2020/T/npf_test_spring2020.csv"))[,-c(1,2,4),]
x[,-1] <- scale(x[,-1])

y <- x[,sapply(colnames(x),
	       function(s,n=nchar(s)) 
		       if(n<4) TRUE else substr(s,n-3,n)!=".std" )]

saveRDS(x,"webdemo/data/hyytiala_with_std.rds")
saveRDS(y,"webdemo/data/hyytiala.rds")


