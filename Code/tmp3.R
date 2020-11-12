#setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)

Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry", "base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
	"TS","Base Tail LD","Tip Tail LD","Nose LD")
Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
"nose_lateral_displacement_phase")

names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Genotype'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))
levels(data_per_stride$Strain)[1] <- "C57BL/6NJ"
levels(data_per_stride$Strain)[119] <- "Mrps22<tm1.1(KOMP)Vlcg> -/+"

#Remove Strains
toMatch <- c("Esrrb", "<em2J>/J COIN", "IMPC")
matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
Strains <- setdiff(unique(data_per_stride$Strain), matches)
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% 
vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
BodyLength <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'BodyLength'][1])
label <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'NetworkFilename'][1])
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
data_per_animal <- cbind(label,Strain, BodyLength, data_per_animal)
data_per_animal$label <- gsub("%2F","/",data_per_animal$label)
Strains <- gsub("*./.*","",unique(Strain))
#Strains[35] <- 'C57BL/6NJ'
Strains <- setdiff(Strains,Strains[table(data_per_animal$Strain) == 1])

allnose.list <- list()
allbt.list <- list()
alltt.list <- list()

for (s in 1:5){

	cat("Analyzing Strain", paste0(Strains[s]), "\n")
	nose <- read.csv(paste0('../Data/Phase_dfs/nose/',Strains[s],'.csv'), header = TRUE)
	nose$stride_num <- as.factor(nose$stride_num)

	#Attaching MouseID to NetworkID
	MouseID <- sapply(seq(dim(nose)[1]), function(x) 
	data_per_animal[data_per_animal$label == nose$label[x], 'MouseID'][1])
	nose <- cbind(nose,MouseID)
	#nose <- nose[nose$MouseID %in% CtrlIDs,]
	nose$MouseID <- droplevels(nose$MouseID)
	
	#Calculating the Fréchet mean for all strides made by every MouseID
	df.list <- list()
	invisible(lapply(seq(length(unique(nose$MouseID))), function(x) {
		df <- nose[nose$MouseID == unique(nose$MouseID)[x],]
		df1 <- df[,names(df) %in% c('Percent.Stride','stride_num','Displacement')]
		tmp <- list()
		tmp[[1]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,-1]
		tmp[[2]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,1]
		names(tmp) <- c('f','time')
		tmp$f <- as.matrix(tmp$f)
		invisible(capture.output(out <- time_warping(tmp$f,tmp$time, method = 'median')))
		df.list[[x]] <<- out$fmean
	}))
	df.tmp <- do.call(cbind,df.list)

	#Calculating the Fréchet mean for the strain line
	invisible(capture.output(out <- time_warping(df.tmp, unique(df$Percent.Stride))))
	allnose.list[[s]] <- out$fmean

}

allnose.df <- do.call(cbind,allnose.list)
allnose.melt <- reshape::melt(allnose.df)
allnose.melt[,1] <- rep(unique(nose$Percent.Stride),5)
allnose.melt[,2] <- rep(Strains[1:5], each = 60)
colnames(allnose.melt) <- c('Percent.Stride','Strain','Displacement')
ggplot(allnose.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line()

s <- 31 #'C57B'
s <- 84 #'Pcdh9'

CtrlIDs <- c('J69863','J69864','J69866','J69867','J69868','J69870','J69871','J69872') 
#'J85758','J85759','J85761','J85762','J85763','J85765','J85766')
labels <- c("LL7-GRB3/2016-10-03/GRB3-1_12-02-54.avi",           
"LL7-GRB3/2016-10-03/GRB3-2_11-25-32.avi",           
"LL7-GRB3/2016-10-03/GRB3-2_12-02-54.avi",           
"LL7-GRB3/2016-10-03/GRB3-3_11-25-32.avi",           
"LL7-GRB3/2016-12-26/GRB3-2_12-26-44.avi",           
"LL7-GRB3/2016-12-26/GRB3-3_09-01-13.avi",           
"LL7-GRB3/2016-12-26/GRB3-4_09-01-13.avi",           
"LL7-GRB3/2016-12-26/GRB3-4_12-26-44.avi",          
"LL8-GRB3/Cropped/2016-10-03/GRB3-10_10-46-04.avi",
"LL8-GRB3/Cropped/2016-10-03/GRB3-10_11-26-45.avi",
"LL8-GRB3/Cropped/2016-10-03/GRB3-12_10-46-04.avi",
"LL8-GRB3/Cropped/2016-10-03/GRB3-9_11-26-45.avi", 
"LL8-GRB3/Cropped/2016-12-26/GRB3-10_09-03-28.avi",
"LL8-GRB3/Cropped/2016-12-26/GRB3-11_09-03-28.avi",
"LL8-GRB3/Cropped/2016-12-26/GRB3-9_12-28-59.avi" )
nose <- nose[nose$label %in% labels,]
nose$label <- droplevels(nose$label)
DF <- df.tmp
DF <- cbind(df.tmp, DF)
DF <- data.frame(DF) 
colnames(DF)[1:3] <- c('J67783','J81952','J81953')
colnames(DF)[4:11] <- c('J69863','J69864','J69866','J69867','J69868','J69870','J69871','J69872')

df <- DF
df.melt <- reshape::melt(df)
df.melt <- cbind(rep(unique(nose$Percent.Stride),11),df.melt) 
colnames(df.melt) <- c('Percent.Stride','MouseID','Displacement')
df.melt$Strain <- ifelse(df.melt$MouseID %in% c('J67783','J81952','J81953'), 'Pcdh9-/+', 'C57BL/6NJ')
ggplot(df.melt, aes(x = Percent.Stride, y = Displacement, color = Strain, fill = MouseID)) + geom_line()

p1 <- ggplot(df.melt, aes(x = Percent.Stride, y = Displacement, color = MouseID)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Data') + 
scale_color_manual(values=getPalette(length(unique(MouseID))))

MouseIDs <- c('J67783','J81952','J81953','J69863','J69864','J69866','J69867','J69868',
	'J69870','J69871','J69872')
df.list <- list()
df.list[["f"]] <- as.matrix(df)
df.list[['time']] <- unique(nose$Percent.Stride) 
out <- fdasrvf::time_warping(df.list$f,df.list$time, method = 'median',MaxItr = 100)
dfa <- as.matrix(out$fn)
colnames(dfa) <- c('J67783','J81952','J81953','J69863','J69864','J69866','J69867','J69868',
	'J69870','J69871','J69872')
dfp <- as.matrix(out$gam)
colnames(dfp) <- c('J67783','J81952','J81953','J69863','J69864','J69866','J69867','J69868',
	'J69870','J69871','J69872')
dfa.melt <- reshape::melt(dfa)
dfa.melt <- dfa.melt[,-1]
dfa.melt <- cbind(rep(unique(nose$Percent.Stride),length(MouseIDs)), dfa.melt)
colnames(dfa.melt) <- c('Percent.Stride','MouseID','Displacement')
dfa.melt$Strain <- ifelse(dfa.melt$MouseID %in% c('J67783','J81952','J81953'), 'Pcdh9-/+', 'C57BL/6NJ')
p5 <- ggplot(dfa.melt, aes(x = Percent.Stride, y = Displacement, fill = MouseID, color = Strain)) + geom_line(lwd = 1) + theme_bw(base_size = 22) + 
theme(legend.position = 'top') + labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + 
scale_color_manual(values=c("C57BL/6NJ" = "black","Pcdh9-/+" = "#d94801"))
#dev.print(pdf,'../Temp6/nose-amp.pdf', width = 4, height = 4)

dfp.melt <- reshape::melt(dfp)
dfp.melt <- dfp.melt[,-1]
dfp.melt <- cbind(rep(unique(nose$Percent.Stride),length(MouseIDs)), dfp.melt)
colnames(dfp.melt) <- c('Percent.Stride','MouseID','Displacement')
dfp.melt$Strain <- ifelse(dfp.melt$MouseID %in% c('J67783','J81952','J81953'), 'Pcdh9-/+', 'C57BL/6NJ')
p6 <- ggplot(dfp.melt, aes(x = Percent.Stride, y = Displacement, fill = MouseID, color = Strain)) + geom_line(lwd = 1) + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Nose Displacement', x = 'Percent Stride') + 
labs(y = TeX('$\\gamma$')) + scale_color_manual(values=c("C57BL/6NJ" = "black","Pcdh9-/+" = "#d94801"))

(p1/p3)|(p2/p4)



dev.print(pdf,'../Temp6/Pcdh9.pdf', height = 10, width = 6)

#######################
s <- 84 #'Pcdh9'

cat("Analyzing Strain", paste0(Strains[s]), "\n")
nose <- read.csv(paste0('../Data/Phase_dfs/nose/',Strains[s],'.csv'), header = TRUE)
nose$stride_num <- as.factor(nose$stride_num)

#Attaching MouseID to NetworkID
MouseID <- sapply(seq(dim(nose)[1]), function(x) 
data_per_animal[data_per_animal$label == nose$label[x], 'MouseID'][1])
nose <- cbind(nose,MouseID)
#nose <- nose[nose$MouseID %in% CtrlIDs,]
nose$MouseID <- droplevels(nose$MouseID)
	
#Calculating the Fréchet mean for all strides made by every MouseID
df.list <- list()
invisible(lapply(seq(length(unique(nose$MouseID))), function(x) {
	df <- nose[nose$MouseID == unique(nose$MouseID)[x],]
	df1 <- df[,names(df) %in% c('Percent.Stride','stride_num','Displacement')]
	tmp <- list()
	tmp[[1]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,-1]
	tmp[[2]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,1]
	names(tmp) <- c('f','time')
	tmp$f <- as.matrix(tmp$f)
	df.list[[x]] <<- t(tmp$f)
}))
df.tmp <- do.call(rbind,df.list)
df <- df.tmp
MutantCounts <- sapply(seq(length(unique(nose$MouseID))), function(x) 
	length(unique(subset(nose, MouseID == unique(nose$MouseID)[x])$stride_num)))
MutantIDs <- c('J67783','J81952','J81953')

s <- 31 #'C57B'
CtrlIDs <- c('J69863','J69864','J69866','J69867','J69868','J69870','J69871','J69872',
'J85758','J85759','J85761','J85762','J85763','J85765','J85766')
cat("Analyzing Strain", paste0(Strains[s]), "\n")
nose <- read.csv(paste0('../Data/Phase_dfs/nose/',Strains[s],'.csv'), header = TRUE)
nose$stride_num <- as.factor(nose$stride_num)
labels <- c("LL7-GRB3/2016-10-03/GRB3-1_12-02-54.avi",           
"LL7-GRB3/2016-10-03/GRB3-2_11-25-32.avi",           
"LL7-GRB3/2016-10-03/GRB3-2_12-02-54.avi",           
"LL7-GRB3/2016-10-03/GRB3-3_11-25-32.avi",           
"LL7-GRB3/2016-12-26/GRB3-2_12-26-44.avi",           
"LL7-GRB3/2016-12-26/GRB3-3_09-01-13.avi",           
"LL7-GRB3/2016-12-26/GRB3-4_09-01-13.avi",           
"LL7-GRB3/2016-12-26/GRB3-4_12-26-44.avi",          
"LL8-GRB3/Cropped/2016-10-03/GRB3-10_10-46-04.avi",
"LL8-GRB3/Cropped/2016-10-03/GRB3-10_11-26-45.avi",
"LL8-GRB3/Cropped/2016-10-03/GRB3-12_10-46-04.avi",
"LL8-GRB3/Cropped/2016-10-03/GRB3-9_11-26-45.avi", 
"LL8-GRB3/Cropped/2016-12-26/GRB3-10_09-03-28.avi",
"LL8-GRB3/Cropped/2016-12-26/GRB3-11_09-03-28.avi",
"LL8-GRB3/Cropped/2016-12-26/GRB3-9_12-28-59.avi" )
nose <- nose[nose$label %in% labels,]
nose$label <- droplevels(nose$label)


#Attaching MouseID to NetworkID
MouseID <- sapply(seq(dim(nose)[1]), function(x) 
data_per_animal[data_per_animal$label == nose$label[x], 'MouseID'][1])
nose <- cbind(nose,MouseID)
nose <- nose[nose$MouseID %in% CtrlIDs,]
nose$MouseID <- droplevels(nose$MouseID)
	
#Calculating the Fréchet mean for all strides made by every MouseID
df.list <- list()
invisible(lapply(seq(length(unique(nose$MouseID))), function(x) {
	df <- nose[nose$MouseID == unique(nose$MouseID)[x],]
	df1 <- df[,names(df) %in% c('Percent.Stride','stride_num','Displacement')]
	tmp <- list()
	tmp[[1]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,-1]
	tmp[[2]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,1]
	names(tmp) <- c('f','time')
	tmp$f <- as.matrix(tmp$f)
	df.list[[x]] <<- t(tmp$f)
}))
df.tmp <- do.call(rbind,df.list)
CtrlCounts <- sapply(seq(length(unique(nose$MouseID))), function(x) 
	length(unique(subset(nose, MouseID == unique(nose$MouseID)[x])$stride_num)))

df <- rbind(df,df.tmp)
df <- data.frame(df) 
df <- cbind(MouseID = rep(c(MutantIDs,CtrlIDs), c(MutantCounts, CtrlCounts)), df)
BodyLength <- sapply(seq(length(unique(df$MouseID))), function(x) 
	data_per_animal[data_per_animal$MouseID == unique(data_per_animal$MouseID)[x], 'BodyLength'])
df <- cbind(BodyLength = rep(BodyLength, c(MutantCounts,CtrlCounts)), df)
colnames(df) <- c('BodyLength','MouseID',sapply(seq(ncol(df[,-c(1,2)])), function(x) paste0("nose_",x)))
DF <- list()
DF[['nose_df']] <- df[,-c(1,2)]
DF[['BL']] <- df[,1]
tmp <- fosr(nose_df ~ BL, data = DF)

##################################################################
df_basetail <- read.csv('/Users/sabnig/Documents/Projects/Komp/Data/Phase_dfs/Pcdh9_bt.csv')
MouseID <- sapply(seq(dim(df_basetail)[1]), function(x) 
	data_per_animal[data_per_animal$label == df_basetail$label[x], 'MouseID'][1])
df_basetail <- cbind(df_basetail,MouseID)
p2 <- ggplot(df_basetail, aes(x = Percent.Stride, y = Displacement, color = MouseID)) + stat_smooth() + 
labs(x = 'Percent Stride', y = 'Displacement (Base Tail)') + theme_bw(base_size=22) + theme(legend.position='none')

df_tiptail <- read.csv('/Users/sabnig/Documents/Projects/Komp/Data/Phase_dfs/Pcdh9_tt.csv')
MouseID <- sapply(seq(dim(df_tiptail)[1]), function(x) 
	data_per_animal[data_per_animal$label == df_tiptail$label[x], 'MouseID'][1])
df_tiptail <- cbind(df_tiptail,MouseID)
p3 <- ggplot(df_tiptail, aes(x = Percent.Stride, y = Displacement, color = MouseID)) + stat_smooth() + 
labs(x = 'Percent Stride', y = 'Displacement (Tip Tail)') + theme_bw(base_size=22) + theme(legend.position='top')

df_nose <- read.csv('/Users/sabnig/Documents/Projects/Komp/Data/Phase_dfs/Pcdh9_nose.csv')
MouseID <- sapply(seq(dim(df_nose)[1]), function(x) 
	data_per_animal[data_per_animal$label == df_nose$label[x], 'MouseID'][1])
df_nose <- cbind(df_nose,MouseID)
p1 <- ggplot(df_nose, aes(x = Percent.Stride, y = Displacement, color = MouseID)) + stat_smooth() + 
labs(x = 'Percent Stride', y = 'Displacement (Nose)') + theme_bw(base_size=22) + theme(legend.position='none')

(p1 | p2)/p3

df <- read.csv('/Users/sabnig/Documents/Projects/Komp/Temp/tmp.csv')
MouseID <- sapply(seq(dim(df)[1]), function(x) 
	data_per_animal[data_per_animal$label == df$label[x], 'MouseID'][1])
df <- cbind(df,MouseID)
df$MouseID <- droplevels(df$MouseID)
df$stride_num <- as.factor(df$stride_num)

df.tmp <- list()
lapply(seq(length(unique(df$MouseID))), function(x) {
	df0 <- df[df$MouseID == unique(df$MouseID)[x],]
	df0$MouseID <- droplevels(df0$MouseID)
	df0$stride_num <- droplevels(df0$stride_num)
	df1 <- df0[,names(df0) %in% c('Percent.Stride','stride_num','Displacement')]
	tmp <- list()
	tmp[[1]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,-1]
	tmp[[2]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,1]
	names(tmp) <- c('f','time')
	tmp$f <- as.matrix(tmp$f)
	out <- fdasrvf::time_warping(tmp$f,tmp$time, method = 'median')
	df.tmp[[x]] <<- out$fmean
})

df.tmp2 <- do.call(cbind,df.tmp)
colnames(df.tmp2) <- unique(df$MouseID)
out <- align_fPCA(df.tmp2,unique(df$Percent.Stride))
out2 <- time_warping(df.tmp2, unique(df$Percent.Stride))

df0 <- df[df$MouseID == 'J74086',]
df0$MouseID <- droplevels(df0$MouseID)
df0$stride_num <- droplevels(df0$stride_num)
df1 <- df0[,names(df0) %in% c('Percent.Stride','stride_num','Displacement')]
tmp <- list()
tmp[[1]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,-1]
tmp[[2]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,1]
names(tmp) <- c('f','time')
tmp$f <- as.matrix(tmp$f)
out <- fdasrvf::align_fPCA(tmp$f, tmp$time)
out2 <- fdasrvf::time_warping(tmp$f,tmp$time, method = 'median')
df.tmp <- data.frame(fmean = out2$fmean, time = tmp$time)
ggplot(df.tmp, aes(x = time, y = fmean)) + geom_line()