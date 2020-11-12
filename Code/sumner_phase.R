setwd("/projects/kumar-lab/sabnig/Projects/KOMP-Analysis")
data_per_stride <- read.delim('/Data/kompdf-corr', stringsAsFactors = FALSE)

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
label <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'NetworkFilename'][1])
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
data_per_animal <- cbind(label,Strain, data_per_animal)
data_per_animal$label <- gsub("%2F","/",data_per_animal$label)
Strains <- gsub("*./.*","",unique(Strain))
Strains[35] <- 'C57BL/6NJ'
Strains <- setdiff(Strains,Strains[table(data_per_animal$Strain) == 1])

allnose.list <- list()
allbt.list <- list()
alltt.list <- list()

for (s in length(Strains)){

	cat("Analyzing Strain", paste0(Strains[s]), "\n")

	#Nose 
	nose <- read.csv(paste0('../Data/Phase_dfs/nose/',Strains[s],'.csv'), header = TRUE)
	nose$stride_num <- as.factor(nose$stride_num)

	#Attaching MouseID to NetworkID
	MouseID <- sapply(seq(dim(nose)[1]), function(x) 
	data_per_animal[data_per_animal$label == nose$label[x], 'MouseID'][1])
	nose <- cbind(nose,MouseID)
	nose$MouseID <- droplevels(nose$MouseID)
	
	#Calculating the Fréchet mean for all strides made by every MouseID
	df.list <- list()
	invisible(lapply(seq(length(unique(nose$MouseID))), function(x) {
		df <- nose[nose$MouseID == unique(nose$MouseID)[x],]
		df1 <- df[,names(df) %in% c('Percent.Stride','stride_num','Displacement')]
		df1$stride_num <- droplevels(df1$stride_num)
		tmp <- list()
		tmp[[1]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,-1]
		tmp[[2]] <- reshape2::dcast(data = df1, formula = Percent.Stride ~ stride_num, value.var = 'Displacement')[,1]
		names(tmp) <- c('f','time')
		tmp$f <- as.matrix(tmp$f)
		tryCatch({
			invisible(capture.output(out <- time_warping(tmp$f,tmp$time, method = 'median',showplot=FALSE)))
		}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
		df.list[[x]] <<- out$fmean
	}))
	df.tmp <- do.call(cbind,df.list)

	#Calculating the Fréchet mean for the strain line
	invisible(capture.output(out <- time_warping(df.tmp, unique(nose$Percent.Stride), showplot=FALSE)))
	allnose.list[[s]] <- out$fmean


	#Tip Tail
	tip_tail <- read.csv(paste0('../Data/Phase_dfs/tip_tail/',Strains[s],'.csv'), header = TRUE)
	tip_tail$stride_num <- as.factor(tip_tail$stride_num)

	#Attaching MouseID to NetworkID
	MouseID <- sapply(seq(dim(tip_tail)[1]), function(x) 
	data_per_animal[data_per_animal$label == tip_tail$label[x], 'MouseID'][1])
	tip_tail <- cbind(tip_tail,MouseID)
	tip_tail$MouseID <- droplevels(tip_tail$MouseID)
	
	#Calculating the Fréchet mean for all strides made by every MouseID
	df.list <- list()
	invisible(lapply(seq(length(unique(tip_tail$MouseID))), function(x) {
		df <- tip_tail[tip_tail$MouseID == unique(tip_tail$MouseID)[x],]
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
	alltt.list[[s]] <- out$fmean


	#Base Tail
	base_tail <- read.csv(paste0('../Data/Phase_dfs/base_tail/',Strains[s],'.csv'), header = TRUE)
	base_tail$stride_num <- as.factor(base_tail$stride_num)

	#Attaching MouseID to NetworkID
	MouseID <- sapply(seq(dim(base_tail)[1]), function(x) 
	data_per_animal[data_per_animal$label == base_tail$label[x], 'MouseID'][1])
	base_tail <- cbind(base_tail,MouseID)
	base_tail$MouseID <- droplevels(base_tail$MouseID)
	
	#Calculating the Fréchet mean for all strides made by every MouseID
	df.list <- list()
	invisible(lapply(seq(length(unique(base_tail$MouseID))), function(x) {
		df <- base_tail[base_tail$MouseID == unique(base_tail$MouseID)[x],]
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
	allbt.list[[s]] <- out$fmean

}

df.nose <- data.frame(do.call(cbind,allnose.list))
df.bt <- data.frame(do.call(cbind,allbt.list)) 
df.tt <- data.frame(do.call(cbind,alltt.list)) 

write.csv(x = df.nose, file = 'Output/nose.csv')
write.csv(x = df.bt, file = 'Output/basetail.csv')
write.csv(x = df.tt, file = 'Output/tiptail.csv')
