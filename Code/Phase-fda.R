Phenos.lin <- c("angular_velocity","speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry","base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")

setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
df <- read.csv('../Data/Output/nose.csv', stringsAsFactors = FALSE)
df <- df[,-1]
nose <- read.csv(paste0('../Data/Phase_dfs/nose/C57B.csv'), header = TRUE)

data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)
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
label <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'NetworkFilename'][1])
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
data_per_animal <- cbind(label,Strain, data_per_animal)
data_per_animal$label <- gsub("%2F","/",data_per_animal$label)
Strains <- gsub("*./.*","",unique(Strain))
#Strains[35] <- 'C57BL/6NJ'
Strains <- setdiff(Strains,Strains[table(data_per_animal$Strain) == 1])
StrainsE <- c('Rab11fip3','Rab22a','Hdhd1a','Myo1d','Ehmt2','Adam32','Rhbdf1','Whamm','Stk32c') 
Strains <- setdiff(Strains, StrainsE) #'Mtmr3','Fbxo30'
colnames(df) <- Strains

df.melt <- reshape::melt(df)
df.melt <- cbind(rep(unique(nose$Percent.Stride),length(Strains)), df.melt)
colnames(df.melt) <- c('Percent.Stride','Strain','Displacement')
getPalette <- colorRampPalette(brewer.pal(8, "Spectral"))
p1 <- ggplot(df.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 32) + 
theme(legend.position = 'none') + labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Data') + 
scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-original.pdf', width = 4, height = 4)

df.list <- list()
df.list[["f"]] <- as.matrix(df)
df.list[['time']] <- unique(nose$Percent.Stride) 

out <- fdasrvf::time_warping(df.list$f,df.list$time, method = 'median',MaxItr = 100)
dfa <- as.matrix(out$fn)
colnames(dfa) <- Strains
dfp <- as.matrix(out$gam)
colnames(dfp) <- Strains

dfa.melt <- reshape::melt(dfa)
dfa.melt <- dfa.melt[,-1]
dfa.melt <- cbind(rep(unique(nose$Percent.Stride),length(Strains)), dfa.melt)
colnames(dfa.melt) <- c('Percent.Stride','Strain','Displacement')
p2 <- ggplot(dfa.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 32) + 
theme(legend.position = 'none') + labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Amplitude') + 
labs(y = NULL, x = NULL) + scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-amp.pdf', width = 4, height = 4)

dfp.melt <- reshape::melt(dfp)
dfp.melt <- dfp.melt[,-1]
dfp.melt <- cbind(rep(unique(nose$Percent.Stride),length(Strains)), dfp.melt)
colnames(dfp.melt) <- c('Percent.Stride','Strain','Displacement')
p3 <- ggplot(dfp.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 32) + 
theme(legend.position = 'none') + labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Phase') + 
labs(y = NULL, x = NULL) + scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-phase.pdf', width = 4, height = 4)

amp_box <- fdasrvf::AmplitudeBoxplot(out,showplot=FALSE)
tmp <- data.frame(index = amp_box$outlier_index,Strain = Strains[amp_box$outlier_index])
dfbox <- data.frame(min = amp_box$minn, Q1 = amp_box$Q1, 
	M = amp_box$median_y, Q3 = amp_box$Q3, max = amp_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
palette_Dark2 <- colorRampPalette(brewer.pal(14, "Dark2"))
p4 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Summary)) + geom_line(lwd=1) + theme_bw(base_size = 32) + 
labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Amplitude Boxplot') + 
theme(legend.position = 'none',legend.text=element_text(size=12)) + 
discrete_scale("color", "manual", palette_Dark2) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/amp-box.pdf', width = 6, height = 6)

#phase_box <- fdasrvf::PhaseBoxplot(out,showplot=FALSE)
#tmp <- data.frame(index = phase_box$outlier_index,Strain = Strains[phase_box$outlier_index])
dfbox <- data.frame(Pcdh9 = out$fn[,79], Cdkl4 = out$fn[,199], 
	Sipa1l1 = out$fn[,98], Slc43a3 = out$fn[,42], Tmem222 = out$fn[,6])
dfsumm <- data.frame(min = amp_box$minn, Q1 = amp_box$Q1, 
	M = amp_box$median_y, Q3 = amp_box$Q3, max = amp_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
dfbox.melt <- cbind(dfbox.melt, min = dfsumm$min, max = dfsumm$max)
colnames(dfbox.melt)[colnames(dfbox.melt)=='variable'] <- 'Strain'
#colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
#dfbox.melt$Type <- as.factor(ifelse(dfbox.melt$Summary %in% tmp$Strain, paste0(dfbox.melt$Summary), 'Summary'))
p6 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Strain)) + geom_line(lwd=1) + theme_bw(base_size = 32) + 
labs(y = 'Nose Displacement', x = 'Percent Stride') + geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2, color ='grey') + 
ggtitle('Amplitude:Outliers') + 
scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#6a3d9a','#ff7f00')) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/amp-box-out.pdf', width = 6, height = 6)

phase_box <- fdasrvf::PhaseBoxplot(out,showplot=FALSE)
dfbox <- data.frame(min = phase_box$minn, Q1 = phase_box$Q1, 
	M = phase_box$median_x, Q3 = phase_box$Q3, max = phase_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
p5 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Summary)) + geom_line(lwd=1) + theme_bw(base_size = 32) + 
labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Phase Boxplot') + theme(legend.position = 'none') + 
discrete_scale("color", "manual", palette_Dark2) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/phase-box.pdf', width = 6, height = 6)

tmp <- data.frame(index = phase_box$outlier_index,Strain = Strains[phase_box$outlier_index])
dfbox <- data.frame(Pcdh9 = out$gam[,79], Cdkl4 = out$gam[,199], 
	Sipa1l1 = out$gam[,98], Slc43a3 = out$gam[,42], Tmem222 = out$gam[,6],Epb41l1 = out$gam[,210],
	Wdr24 = out$gam[,130], Poc1a = out$gam[,8])
dfsumm <- data.frame(min = phase_box$minn, Q1 = phase_box$Q1, 
	M = phase_box$median_x, Q3 = phase_box$Q3, max = phase_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
dfbox.melt <- cbind(dfbox.melt, min = dfsumm$min, max = dfsumm$max)
colnames(dfbox.melt)[colnames(dfbox.melt)=='variable'] <- 'Strain'
#colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
#dfbox.melt$Type <- as.factor(ifelse(dfbox.melt$Summary %in% tmp$Strain, paste0(dfbox.melt$Summary), 'Summary'))
p7 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Strain)) + geom_line(lwd=1) + 
theme_bw(base_size = 32) + geom_ribbon(aes(ymin=min,ymax=max),alpha = 0.2, color = 'grey') + 
labs(y = 'Nose Displacement', x = 'Percent Stride') + ggtitle('Phase:Outliers') + 
scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#6a3d9a','#ff7f00','#b15928','#f781bf',
'#a65628')) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/phase-box-out.pdf', width = 6, height = 6)

(p2/p3) | (p4/p5) | (p6/p7)
dev.print(pdf,'../Temp7/nose-fda-2.pdf', width = 22, height = 12)
dev.print(pdf,'../Temp7/nose-data.pdf', width = 9, height = 9)

####Base Tail 
df <- read.csv('../Data/base_tail.csv', stringsAsFactors = FALSE)
df <- df[,-1]
base_tail <- read.csv(paste0('../Data/Phase_dfs/base_tail/C57B.csv'), header = TRUE)

data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)
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
label <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'NetworkFilename'][1])
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
data_per_animal <- cbind(label,Strain, data_per_animal)
data_per_animal$label <- gsub("%2F","/",data_per_animal$label)
Strains <- gsub("*./.*","",unique(Strain))
#Strains[35] <- 'C57BL/6NJ'
Strains <- setdiff(Strains,Strains[table(data_per_animal$Strain) == 1])
StrainsE <- c('Rab11fip3','Rab22a','Hdhd1a','Myo1d','Ehmt2','Adam32','Rhbdf1','Whamm','Stk32c') 
Strains <- setdiff(Strains, StrainsE) #'Mtmr3','Fbxo30'
colnames(df) <- Strains

df.melt <- reshape::melt(df)
df.melt <- cbind(rep(unique(nose$Percent.Stride),length(Strains)), df.melt)
colnames(df.melt) <- c('Percent.Stride','Strain','Displacement')
getPalette <- colorRampPalette(brewer.pal(8, "Spectral"))
p1 <- ggplot(df.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Data') + 
scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-original.pdf', width = 4, height = 4)

df.list <- list()
df.list[["f"]] <- as.matrix(df)
df.list[['time']] <- unique(base_tail$Percent.Stride) 

out <- fdasrvf::time_warping(df.list$f,df.list$time, method = 'median',MaxItr = 100)
dfa <- as.matrix(out$fn)
colnames(dfa) <- Strains
dfp <- as.matrix(out$gam)
colnames(dfp) <- Strains

dfa.melt <- reshape::melt(dfa)
dfa.melt <- dfa.melt[,-1]
dfa.melt <- cbind(rep(unique(base_tail$Percent.Stride),length(Strains)), dfa.melt)
colnames(dfa.melt) <- c('Percent.Stride','Strain','Displacement')
p2 <- ggplot(dfa.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Amplitude') + 
labs(y = NULL, x = NULL) + scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-amp.pdf', width = 4, height = 4)

dfp.melt <- reshape::melt(dfp)
dfp.melt <- dfp.melt[,-1]
dfp.melt <- cbind(rep(unique(base_tail$Percent.Stride),length(Strains)), dfp.melt)
colnames(dfp.melt) <- c('Percent.Stride','Strain','Displacement')
p3 <- ggplot(dfp.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Phase') + 
labs(y = NULL, x = NULL) + scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-phase.pdf', width = 4, height = 4)

amp_box <- fdasrvf::AmplitudeBoxplot(out,showplot=FALSE)
tmp <- data.frame(index = amp_box$outlier_index,Strain = Strains[amp_box$outlier_index])
dfbox <- data.frame(min = amp_box$minn, Q1 = amp_box$Q1, 
	M = amp_box$median_y, Q3 = amp_box$Q3, max = amp_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(base_tail$Percent.Stride), dfbox.melt)
colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
p4 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Summary)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Amplitude Boxplot') + theme(legend.position = 'none') + 
discrete_scale("color", "manual", palette_Dark2) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/amp-box.pdf', width = 6, height = 6)

tmp <- data.frame(index = amp_box$outlier_index,Strain = Strains[amp_box$outlier_index])
dfbox <- data.frame(Pcdh9 = out$fn[,79], Gpr87 = out$fn[,54])
dfsumm <- data.frame(min = amp_box$minn, Q1 = amp_box$Q1, 
	M = amp_box$median_y, Q3 = amp_box$Q3, max = amp_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
dfbox.melt <- cbind(dfbox.melt, min = dfsumm$min, max = dfsumm$max)
colnames(dfbox.melt)[colnames(dfbox.melt)=='variable'] <- 'Strain'
#dfbox.melt <- reshape::melt(dfbox)
#dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
#colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
#dfbox.melt$Type <- as.factor(ifelse(dfbox.melt$Summary %in% tmp$Strain, paste0(dfbox.melt$Summary), 'Summary'))
p6 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Strain)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Amplitude:Outliers') + 
geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2, color ='grey') + 
scale_color_manual(values = c('#e41a1c','#377eb8','#999999')) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/amp-box-out.pdf', width = 6, height = 6)

phase_box <- fdasrvf::PhaseBoxplot(out,showplot=FALSE)
dfbox <- data.frame(min = phase_box$minn, Q1 = phase_box$Q1, 
	M = phase_box$median_x, Q3 = phase_box$Q3, max = phase_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
p5 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Summary)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Phase Boxplot') + theme(legend.position = 'none') + 
discrete_scale("color", "manual", palette_Dark2) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/phase-box.pdf', width = 6, height = 6)

tmp <- data.frame(index = phase_box$outlier_index,Strain = Strains[phase_box$outlier_index])
dfbox <- data.frame(Poc1a = out$gam[,8], Hoxc12 = out$gam[,141],
	Epb41l1 = out$gam[,106])
dfsumm <- data.frame(min = phase_box$minn, Q1 = phase_box$Q1, 
	M = phase_box$median_x, Q3 = phase_box$Q3, max = phase_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
dfbox.melt <- cbind(dfbox.melt, min = dfsumm$min, max = dfsumm$max)
colnames(dfbox.melt)[colnames(dfbox.melt)=='variable'] <- 'Strain'
#dfbox.melt <- reshape::melt(dfbox)
#dfbox.melt <- cbind('Percent.Stride' = unique(nose$Percent.Stride), dfbox.melt)
#colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
#dfbox.melt$Type <- as.factor(ifelse(dfbox.melt$Summary %in% tmp$Strain, paste0(dfbox.melt$Summary), 'Summary'))
p7 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
labs(y = 'Base Tail Displacement', x = 'Percent Stride') + ggtitle('Phase:Outliers') + 
geom_ribbon(aes(ymin=min,ymax=max),alpha = 0.2, color = 'grey') + 
scale_color_manual(values = c('#33a02c','#ff7f00','#6a3d9a')) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/phase-box-out.pdf', width = 6, height = 6)

(p2/p3) | (p4/p5) | (p6/p7)
dev.print(pdf,'../Temp6/bt-fda-2.pdf', width = 16, height = 9)
dev.print(pdf,'../Temp6/bt-data.pdf', width = 9, height = 9)


####Tip Tail 
df <- read.csv('../Data/tip_tail.csv', stringsAsFactors = FALSE)
df <- df[,-1]
tip_tail <- read.csv(paste0('../Data/Phase_dfs/tip_tail/C57B.csv'), header = TRUE)

data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)
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
label <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'NetworkFilename'][1])
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
data_per_animal <- cbind(label,Strain, data_per_animal)
data_per_animal$label <- gsub("%2F","/",data_per_animal$label)
Strains <- gsub("*./.*","",unique(Strain))
#Strains[35] <- 'C57BL/6NJ'
Strains <- setdiff(Strains,Strains[table(data_per_animal$Strain) == 1])
StrainsE <- c('Rab11fip3','Rab22a','Hdhd1a','Myo1d','Ehmt2','Adam32','Rhbdf1','Whamm','Stk32c') 
Strains <- setdiff(Strains, StrainsE) #'Mtmr3','Fbxo30'
colnames(df) <- Strains

df.melt <- reshape::melt(df)
df.melt <- cbind(rep(unique(tip_tail$Percent.Stride),length(Strains)), df.melt)
colnames(df.melt) <- c('Percent.Stride','Strain','Displacement')
getPalette <- colorRampPalette(brewer.pal(8, "Spectral"))
p1 <- ggplot(df.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Data') + 
scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-original.pdf', width = 4, height = 4)

df.list <- list()
df.list[["f"]] <- as.matrix(df)
df.list[['time']] <- unique(tip_tail$Percent.Stride) 

#out <- fdasrvf::time_warping(df.list$f,df.list$time, method = 'median',MaxItr = 200)
load('../Data/Output/tip_tail_strain.RData')
dfa <- as.matrix(out$fn)
colnames(dfa) <- Strains
dfp <- as.matrix(out$gam)
colnames(dfp) <- Strains

dfa.melt <- reshape::melt(dfa)
dfa.melt <- dfa.melt[,-1]
dfa.melt <- cbind(rep(unique(tip_tail$Percent.Stride),length(Strains)), dfa.melt)
colnames(dfa.melt) <- c('Percent.Stride','Strain','Displacement')
p2 <- ggplot(dfa.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Amplitude') + 
labs(y = NULL, x = NULL) + scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-amp.pdf', width = 4, height = 4)

dfp.melt <- reshape::melt(dfp)
dfp.melt <- dfp.melt[,-1]
dfp.melt <- cbind(rep(unique(tip_tail$Percent.Stride),length(Strains)), dfp.melt)
colnames(dfp.melt) <- c('Percent.Stride','Strain','Displacement')
p3 <- ggplot(dfp.melt, aes(x = Percent.Stride, y = Displacement, color = Strain)) + geom_line() + theme_bw(base_size = 22) + 
theme(legend.position = 'none') + labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Phase') + 
labs(y = NULL, x = NULL) + scale_color_manual(values=getPalette(length(Strains)))
#dev.print(pdf,'../Temp6/nose-phase.pdf', width = 4, height = 4)

amp_box <- fdasrvf::AmplitudeBoxplot(out,showplot=FALSE)
tmp <- data.frame(index = amp_box$outlier_index,Strain = Strains[amp_box$outlier_index])
dfbox <- data.frame(min = amp_box$minn, Q1 = amp_box$Q1, 
	M = amp_box$median_y, Q3 = amp_box$Q3, max = amp_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(tip_tail$Percent.Stride), dfbox.melt)
colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
palette_Dark2 <- colorRampPalette(brewer.pal(14, "Dark2"))
p4 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Summary)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Amplitude Boxplot') + theme(legend.position = 'none') + 
discrete_scale("color", "manual", palette_Dark2) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/amp-box.pdf', width = 6, height = 6)

tmp <- data.frame(index = amp_box$outlier_index,Strain = Strains[amp_box$outlier_index])
dfbox <- data.frame(Pcdh9 = out$fn[,79], Gpr87 = out$fn[,54],
	Slc43a3 = out$fn[,42], Adgra1 = out$fn[,9], Mrgprd = out$fn[,55], Lamp5 = out$fn[,180])
dfsumm <- data.frame(min = amp_box$minn, Q1 = amp_box$Q1, 
	M = amp_box$median_y, Q3 = amp_box$Q3, max = amp_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(tip_tail$Percent.Stride), dfbox.melt)
dfbox.melt <- cbind(dfbox.melt, min = dfsumm$min, max = dfsumm$max)
colnames(dfbox.melt)[colnames(dfbox.melt)=='variable'] <- 'Strain'
p6 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Strain)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Amplitude:Outliers') + 
geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2, color ='grey') + 
scale_color_manual(values = c('#e41a1c','#377eb8','#f781bf','#a65628','#4daf4a','#984ea3','#999999')) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/amp-box-out.pdf', width = 6, height = 6)

phase_box <- fdasrvf::PhaseBoxplot(out,showplot=FALSE)
dfbox <- data.frame(min = phase_box$minn, Q1 = phase_box$Q1, 
	M = phase_box$median_x, Q3 = phase_box$Q3, max = phase_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(tip_tail$Percent.Stride), dfbox.melt)
colnames(dfbox.melt) <- c('Percent.Stride', 'Summary', 'value')
p5 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Summary)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Phase Boxplot') + theme(legend.position = 'none') + 
discrete_scale("color", "manual", palette_Dark2) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/phase-box.pdf', width = 6, height = 6)

tmp <- data.frame(index = phase_box$outlier_index,Strain = Strains[phase_box$outlier_index])
dfbox <- data.frame(Pcdh9 = out$gam[,79], Hoxc12 = out$gam[,141],
	Slc43a3 = out$gam[,42], Psmc6 = out$gam[,143])
dfsumm <- data.frame(min = phase_box$minn, Q1 = phase_box$Q1, 
	M = phase_box$median_x, Q3 = phase_box$Q3, max = phase_box$maxx)
dfbox.melt <- reshape::melt(dfbox)
dfbox.melt <- cbind('Percent.Stride' = unique(tip_tail$Percent.Stride), dfbox.melt)
dfbox.melt <- cbind(dfbox.melt, min = dfsumm$min, max = dfsumm$max)
colnames(dfbox.melt)[colnames(dfbox.melt)=='variable'] <- 'Strain'
p7 <- ggplot(dfbox.melt, aes(x = Percent.Stride, y = value, color = Strain)) + geom_line(lwd=1) + theme_bw(base_size = 22) + 
labs(y = 'Tip Tail Displacement', x = 'Percent Stride') + ggtitle('Phase:Outliers') + 
geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2, color ='grey') + 
scale_color_manual(values = c('#e41a1c','#fdb462','#f781bf','#8dd3c7')) + labs(y = NULL, x = NULL)
#dev.print(pdf,'../Temp6/phase-box-out.pdf', width = 6, height = 6)

(p2/p3) | (p4/p5) | (p6/p7)

dev.print(pdf,'../Temp6/tt-fda-2.pdf', width = 16, height = 9)
dev.print(pdf,'../Temp6/tt-data.pdf', width = 9, height = 9)