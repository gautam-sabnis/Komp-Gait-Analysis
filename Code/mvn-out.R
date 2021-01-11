setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)

Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry","base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Speed","LDF","Step Length","Step Width","Stride Length",
"TS","Base Tail LD","Tip Tail LD","Nose LD")
Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
"nose_lateral_displacement_phase")

names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Genotype'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
names(data_per_stride)[names(data_per_stride) == "OFA_Arena.ID"] <- 'Arena'
names(data_per_stride)[names(data_per_stride) == "OFA_Experimenter.ID"] <- 'Experimenter'
data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))
levels(data_per_stride$Strain)[1] <- "C57BL/6NJ"
levels(data_per_stride$Strain)[3] <- "Rik1<em1J> -/-"
levels(data_per_stride$Strain)[4] <- "Rik2<tm1.1(KOMP)Vlcg> -/-"
levels(data_per_stride$Strain)[119] <- "Mrps22<tm1.1(KOMP)Vlcg> -/+"

#Remove Strains
toMatch <- c("B6.Cg-Esrrb<tm1(cre)Yba>/J", "<em2J>/J COIN","Tex2")
matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
Strains <- setdiff(unique(data_per_stride$Strain), matches)
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
	'speed_25_ang_vel_neg20'),]
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
data_per_animal <- cbind(Strain, TestDate, data_per_animal)

#Filter Strains for which at least 8 animals were tested 
Strains8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 5]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strains8,] 
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains8,]
data_per_animal$Strain <- droplevels(data_per_animal$Strain)

data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("Strain")], 
	FUN = function(x) mean(x, na.rm=TRUE))
levels(data_per_strain$Strain)[1] <- "C57BL/6NJ"
data_per_strain[1,'Strain'] <- 'C57BL/6NJ'
#<em1J> - em1J, (KOMP) - Mbp, Wtsi, Vlcg, (EUCOMM) - Hmgu 
data_per_animal$BG <- as.factor(ifelse(grepl("Mbp", data_per_animal$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_animal$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_animal$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_animal$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_animal$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))

#C57BL/6NJ      em1J      Hmgu       Mbp      Vlcg      Wtsi
df <- data_per_animal[data_per_animal$BG %in% c('C57BL/6NJ'),]
tmp <- rrcov::CovRobust(df[,names(df) %in% Phenos.lin])
plot(tmp, which='dd')

df <- data_per_strain
tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])

df.out <- data.frame(rodist = tmp$md.rob, mahadist = tmp$md.cla, Strain = df[,'Strain'], 
	Outlier = tmp$outliers)
#(ifelse(mvoutlier::aq.plot(df[,names(df) %in% Phenos.lin])$outliers==TRUE,1,0))
df.out[df.out$Strain == 'C57BL/6NJ', 'Outlier'] <- -1
df.out$Outlier <- as.factor(df.out$Outlier)
df.out$Genotype <- ifelse(df.out$Strain == 'C57BL/6NJ', "C57BL/6NJ", "Mutant")
textdf <- df.out[df.out$Strain == 'C57BL/6NJ', ]
mycolors <- c("C57BL/6NJ" = "red", "Mutant" = "grey50")
ggplot(df.out, aes(x = mahadist, y=rodist)) + geom_point(alpha=0.8, aes(color=Outlier), size=4) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1,-1),as.character(Strain),'')),size=6,box.padding=2) + 
labs(x = 'Mahalanobis Distance', y = 'Robust Distance') + ggtitle('KOMP Outliers') + scale_color_manual(values=c("black","grey50","red")) + 
theme_bw(base_size = 16) + theme(legend.position='none')
ggsave('../Temp3/Mvout-mean-greaterthan8animals.pdf', width=9, height=9)

Mutants.out <- unique(df.out[df.out$Outlier==1,'Strain'])

df <- df[,-1]
names(df) <- Phenos.lin.Nomen
df <- cbind('Strain' = data_per_strain$Strain,df)
df <- cbind(id = 1:dim(df)[1], df)
df.melt <- reshape::melt(df[,-2], id.vars = 'id')
df.melt <- cbind(df.melt, Strain = rep(df$Strain, length(Phenos.lin)), Outlier = as.factor(rep(df.out$Outlier, length(Phenos.lin))))
ggplot(df.melt, aes(y = value, color = Outlier, x = id)) + geom_jitter(alpha=0.7, size = 3) + 
facet_wrap(~variable, scales='free') + scale_color_manual(values = c('black','grey50','red')) + theme_bw(base_size = 16) + 
theme(legend.position = 'none') + labs(x = 'Index', y='Phenotype') 
ggsave('../Temp2/Mvout2.pdf', width=9, height=9)


tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin])
df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, Strain = df[,'Strain'], 
	Outlier = as.numeric(!tmp$wfinal01))
df.out$Genotype <- ifelse(df.out$Strain == 'C57BL/6NJ', "C57BL/6NJ", "Mutant")
textdf <- df.out[df.out$Strain == 'C57BL/6NJ', ]
mycolors <- c("C57BL/6NJ" = "red", "Mutant" = "grey50")
ggplot(df.out, aes(x = Distance1, y=Distance2)) + geom_point(alpha=0.8, aes(color=Genotype), size=4) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier==1,as.character(Strain),'')),size=4,box.padding=2) + 
labs(x = 'Distance1', y = 'Distance2') + ggtitle('KOMP Outliers') + scale_color_manual(values=mycolors) + 
theme_bw(base_size = 16)


gBG <- c('em1J','Hmgu','Mbp','Vlcg','Wtsi')
lapply(seq(gBG), function(x) {
	df <- data_per_animal[data_per_animal$BG %in% c(gBG[x]),]
	tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])
	df.out <- data.frame(rodist = tmp$md.rob, mahadist = tmp$md.cla, MouseID = df[,'MouseID'], 
		Outlier = ifelse(tmp$outliers==TRUE,1,0))
	p1 <- ggplot(df.out, aes(x = mahadist, y=rodist)) + geom_point(alpha=0.8) + 
	ggrepel::geom_text_repel(aes(label=ifelse(Outlier==1,as.character(MouseID),'')),hjust=0,vjust=0) + 
	labs(x = 'Mahalanobis Distance', y = 'Robust Distance') + ggtitle(paste0('Background - ',gBG[x]))
	ggsave(paste0('Mvout-',gBG[x],'.pdf'),p1)
})

Mvoutliers_komp <- function(gBG){
	df <- data_per_animal[data_per_animal$BG %in% c(gBG),]
	tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])
	df.out <- data.frame(rodist = tmp$md.rob, mahadist = tmp$md.cla, MouseID = df[,'MouseID'], 
		Outlier = ifelse(tmp$outliers==TRUE,1,0))
	ggplot(df.out, aes(x = mahadist, y=rodist)) + geom_point(alpha=0.8) + 
	ggrepel::geom_text_repel(aes(label=ifelse(Outlier==1,as.character(MouseID),'')),hjust=0,vjust=0) + 
	labs(x = 'Mahalanobis Distance', y = 'Robust Distance') + ggtitle(paste0('Background - ',gBG))
}


df <- data_per_animal
tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])

df.out <- data.frame(rodist = tmp$md.rob, mahadist = tmp$md.cla, MouseID = df[,'MouseID'], 
	Outlier = tmp$outliers)
df.out[df.out$rodist > 5.5, 'Outlier'] <- -1
#(ifelse(mvoutlier::aq.plot(df[,names(df) %in% Phenos.lin])$outliers==TRUE,1,0))
df.out$Outlier <- as.factor(df.out$Outlier)
ggplot(df.out, aes(x=mahadist, y=rodist)) + geom_point(alpha=0.6, aes(color=Outlier), size=4) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(-1),as.character(MouseID),'')),size=4,box.padding=2) + 
labs(x = 'Mahalanobis Distance', y = 'Robust Distance') + ggtitle('KOMP animal outliers') + scale_color_manual(values=c("red","black","grey50")) + 
theme_bw(base_size = 16) + theme(legend.position='none')

df <- data.frame(Strain = data_per_animal$Strain, MouseID = data_per_animal$MouseID, Sex = data_per_animal$Sex)
df <- cbind(df, apply(data_per_animal[,names(data_per_animal) %in% Phenos.lin],2,function(x) (x - mean(x))/sd(x)))
tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin],makeplot=TRUE)
df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, 
	Label = paste0(df[,'MouseID']," (", df[,'Strain'], ")"), MouseID = df[,'MouseID'], Strain = df[,'Strain'],
	Outlier = as.numeric(!tmp$wfinal01))
df.out[(df.out$Distance2 >= 5.9 | df.out$Distance1 > 10), 'Outlier'] <- -1
df.out[df.out$MouseID == 'J80962', 'Outlier'] <- -1
df.out[df.out$MouseID == 'J76119', 'Outlier'] <- -1
df.out$Outlier <- as.factor(df.out$Outlier)
ggplot(df.out, aes(x = Distance1, y=Distance2)) + geom_point(alpha=0.8, aes(color=Outlier), size=4) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier==-1,as.character(Label),'')),size=5,box.padding=2) + 
labs(x = 'Distance1', y = 'Distance2') + ggtitle('KOMP Outliers') + scale_color_manual(values=c("red","black","grey50")) + 
theme_bw(base_size = 28) + theme(legend.position='none') + ggtitle('KOMP animal outliers')
ggsave('../Temp7/komp-animal-outliers2.pdf', width=36, height=36, limitsize=FALSE) #base_size=98 for plots

#df.out[df.out$Outlier %in% c(1,-1),'Outlier'] <- -1
#df.out[df.out$Strain == 'Rab6b-/-','Outlier'] <- 1
#df.out[df.out$Strain == 'Lamp5-/-','Outlier'] <- 1
#df.out$Outlier <- as.factor(df.out$Outlier)
ggplot(df.out, aes(x = Distance1, y=Distance2)) + geom_point(alpha=0.8, aes(color=Outlier), size=4) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier==1,as.character(Label),'')),size=5,box.padding=2) + 
labs(x = 'Distance1', y = 'Distance2') + ggtitle('KOMP Outliers') + scale_color_manual(values=c("grey50","black","red")) + 
theme_bw(base_size = 28) + theme(legend.position='none') + ggtitle('KOMP animal outliers')


out.df <- df.out[df.out$Outlier %in% c(1,-1), names(df.out) %in% c('Strain','MouseID')]
#out.df <- out.df[out.df$Strain %in% names(table(data_per_animal$Strain)[table(data_per_animal$Strain) < 5]),]
Strains.out <- out.df$Strain
#out.df <- out.df[out.df$Strain %in% setdiff(Strains.out,Strains8),]
#out.df$Strain <- droplevels(out.df$Strain)
out.df$MouseID <- droplevels(out.df$MouseID)
total_animals <- sapply(seq(length(unique(out.df$Strain))), function(x) table(data_per_animal$Strain)[paste0(unique(out.df$Strain)[x])])
out_animals <- sapply(seq(length(unique(out.df$Strain))), function(x) table(out.df$Strain)[paste0(unique(out.df$Strain)[x])])
df <- data.frame(Proportion = out_animals/total_animals,
	Animals = total_animals)
df <- cbind(Strain = rownames(df), df)
df <- df[with(df,order(-Proportion)),]
#df$color <- ifelse(df$Strain %in% c('Pcdh9-/+','Alg11-/+','Sfxn5-/+','Nsun5-/-','Whamm-/-'), "red", "black")
#df$Strain <- glue("<i style='color:{color}'>{Strain}</i>")
#sapply(seq(nrow(out.df)), function(x) length(unique(data_per_animal[data_per_animal$Strain %in% out.df$Strain[x],'MouseID'])))
ggplot(df, aes(x=Strain,y=Proportion)) + geom_bar(stat='identity') + theme_minimal(base_size=22) +  
theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5), 
	axis.text.x.top = element_text(vjust=0.5)) + labs(y = 'Proportion') + 
geom_text(label = paste0(df$Animals),vjust = -0.1, hjust = 0.5,size=5.5) + coord_cartesian(clip = "off")
ggsave('../Temp7/animal-by-strain-prop-outliers.pdf', width=11.75, height=3.8)

mouseID.out <- sapply(seq(length(unique(out.df$Strain))), function(x) 
	out.df[out.df$Strain == unique(out.df$Strain)[x],'MouseID'])

df <- data_per_stride[data_per_stride$MouseID %in% unlist(mouseID.out[1:length(mouseID.out)]), ]
df$MouseID <- droplevels(df$MouseID)
df$Strain <- droplevels(df$Strain)



df <- data_per_stride[data_per_stride$Strain == 'Sfxn5-/+',c('MouseID',Phenos.lin)]
#colnames(df)[which(names(df) %in% Phenos.lin)] <- Phenos.lin.Nomen
df$MouseID <- droplevels(df$MouseID)
#df$Strain <- droplevels(df$Strain)
df$Index <- ave(df$stride_length, df$MouseID, FUN = seq_along)
df$Outlier <- as.factor(ifelse(df$MouseID %in% out.df$MouseID, 1, 0)) 


invisible(lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), 
	ggplot(df, aes_string(x = 'Index', y = Phenos.lin[x])) + geom_line(aes(color=MouseID)) + 
	scale_color_manual(values = c(ifelse(sapply(seq(length(unique(df$MouseID))), function(x) 
		unique(df[df$MouseID==unique(df$MouseID)[x],'Outlier'])) == 1, 'red','grey50'))) + 
	labs(y = paste0(Phenos.lin.Nomen[x])), inherits=TRUE)))
legend <- get_legend(p1)
p <- plot_grid(p1+labs(x=NULL)+theme(legend.position='none'),p2+labs(x=NULL)+theme(legend.position='none'),p3+labs(x=NULL)+theme(legend.position='none'),
	p4+labs(x=NULL)+theme(legend.position='none'),p5+labs(x=NULL)+theme(legend.position='none'),
	p6+labs(x=NULL)+theme(legend.position='none'),p7+labs(x=NULL)+theme(legend.position='none'),
	p8+labs(x=NULL)+theme(legend.position='none'),p9+labs(x=NULL)+theme(legend.position='none'),nrow=9)
plot_grid(p,legend,rel_widths = c(3, .4))
dev.print(pdf,'../Temp4/Sfxn5.pdf')

#Remove outliers 
invisible(sapply(seq(length(unique(data_per_animal$Strain))), function(s) Map(function(p) {
	vals <- data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][[Phenos.lin[p]]];
	outliers <- boxplot.stats(vals)$out
	ids <- match(outliers, vals) 
	data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][paste0(Phenos.lin[p])][ids,] <<- NA

}, seq(length(Phenos.lin)))))


komp_select_controls <- function(CtrlStrain="C57BL/6NJ", tw){
	control.list <- subset(data_per_stride, Strain == CtrlStrain)[,c('MouseID','Strain','TestDate')]
	mutant.list <- subset(data_per_stride, Strain != CtrlStrain)[,c('MouseID','Strain','TestDate')]
	mutant.list$Type <- 'Mutant'
	all.control.dates.df <- NULL
	control.dates.df.tmp <- NULL
	all.control.mouseid.df <- NULL

	for (s in unique(mutant.list$Strain)){
    	mutant.dates.list <- unique(subset(mutant.list,Strain == s)$TestDate)
    	control.dates.list <- NULL
    	mouse.id.list <- NULL
    	for (d in mutant.dates.list){
        	start.d <- as.Date(d) - tw
        	end.d <- as.Date(d) + tw
        	control.dates <- unique(subset(control.list, as.Date(TestDate) >= start.d & as.Date(TestDate) <= end.d)$TestDate)
        	control.dates.list <- c(control.dates.list, format(control.dates, format = '%Y-%m-%d'))
        	mouse.id <- unique(subset(control.list, as.Date(TestDate) >= start.d & as.Date(TestDate) <= end.d)$MouseID)
    		mouse.id.list <- c(mouse.id.list, as.character(mouse.id))
    	}
    	if (identical(mouse.id.list,character(0))) next
    	control.dates.df <- data.frame(Strain = s, TestDate = format(control.dates.list, format = '%Y-%m-%d'), Type = 'Control')
    	all.control.dates.df <- rbind(all.control.dates.df, control.dates.df[,c('Strain','Type')])
    	control.dates.df.tmp <- c(control.dates.df.tmp, format(control.dates.list, format = '%Y-%m-%d'))

    	control.mouseid.df <- data.frame(Strain = s, MouseID = mouse.id.list, Type = 'Control')
    	all.control.mouseid.df <- rbind(all.control.mouseid.df,control.mouseid.df[,c('Strain','MouseID','Type')])
	}

	all.control.dates.df <- cbind(all.control.dates.df, TestDate = format(control.dates.df.tmp, format = '%Y-%m-%d'))
	dates.df <- rbind(all.control.dates.df,mutant.list[,c('Strain','TestDate','Type')])
	dates.df$TestDate <- as.Date(dates.df$TestDate, format = "%Y-%m-%d")
	mouseid.df <- rbind(all.control.mouseid.df, mutant.list[,c('Strain','MouseID','Type')])
	return(mouseid.df)
}

controlids.df <- komp_select_controls(CtrlStrain="C57BL/6NJ", tw=21) #tw: time window
Mutants.out <- c('Pcdh9-/+','Alg11-/+','Nsun5-/-','Sfxn5-/+','Whamm-/-','Hoxc12-/-','Mettl10-/-','Tusc3-/+',
	'Rmnd5b-/-','Msn-/-','Prss8-/+','Zbtb43-/-','Ccdc28a-/-','Pigc-/+')
#'Whamm-/-'

lapply(seq(length(Mutants)), function(m) {
	CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
	dfa <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs, ]
	dfa$Genotype <- ifelse(dfa$Strain == 'C57BL/6NJ', 'Control', 'Mutant')
	dfa$MouseID <- droplevels(dfa$MouseID)
	#dfa <- dfa[dfa$TestDate %in% '2016-09-06', ] #only for Alg11 for plotting purposes
	dfa <- dfa[dfa$TestDate %in% names(which(table(dfa$TestDate, dfa$Genotype)[,2] >= 1)), ] #for everything else
	dfa$TestDate <- droplevels(dfa$TestDate)
	#CtrlIDs <- sample(unique(dfa[dfa$Genotype=='Control','MouseID']), length(unique(dfa[dfa$Genotype=='Mutant','MouseID'])))
	CtrlIDs <- setdiff(unique(dfa[dfa$Genotype == 'Control', 'MouseID']), unique(out.df$MouseID))
	df <- data_per_stride[data_per_stride$Strain == Mutants[m],c('MouseID','BodyLength','Sex',Phenos.lin)]
	df <- rbind(df,data_per_stride[data_per_stride$MouseID %in% CtrlIDs, c('MouseID', 'BodyLength','Sex',Phenos.lin)])
	df$MouseID <- droplevels(df$MouseID)
	df$Outlier <- as.factor(ifelse(df$MouseID %in% out.df$MouseID, 1, 
		ifelse(df$MouseID %in% CtrlIDs, -1, 0))) 
	if (!("0" %in% levels(df$Outlier))){df$Outlier <- factor(df$Outlier, levels = c("0",levels(df$Outlier)))}
	df2 <- df[,names(df) %in% c(Phenos.lin)]
	#df2 <- data.frame(do.call(cbind,lapply(seq(length(Phenos.lin)), function(p) 
	#	as.numeric(resid(lm(df[[Phenos.lin[p]]] ~ BodyLength, data = df))))))
	names(df2) <- Phenos.lin.Nomen
	df2 <- cbind(id = 1:dim(df)[1], df2)
	df.melt <- reshape::melt(df2, id.vars = 'id')
	df.melt <- cbind(df.melt, MouseID = rep(rep(names(table(df$MouseID)), as.numeric(table(df$MouseID))),length(Phenos.lin)))
	df.melt$Outlier <- as.factor(ifelse(df.melt$MouseID %in% out.df$MouseID, 1, 
		ifelse(df.melt$MouseID %in% CtrlIDs, -1, 0))) 
	if (!("0" %in% levels(df.melt$Outlier))){df.melt$Outlier <- factor(df.melt$Outlier, levels = c("0",levels(df.melt$Outlier)))}
	p2 <- ggplot(df.melt[df.melt$variable %in% c('Base Tail LD','Tip Tail LD','Nose LD'),], aes(x=MouseID,y=value)) + 
	geom_boxplot(outlier.shape=NA,aes(fill = Outlier),alpha = 0.4) + 
	facet_wrap(~ variable, scales = 'free') + ggtitle(paste0(Mutants[m])) + 
	scale_fill_manual(name = 'Genotype', values =c("0" = "black","1" = "#d94801", "-1" = "#6a51a3"), 
		labels = c("0"='Mutant',"-1"='Control',"1"="Mutant (Outlier)"), drop = FALSE) + theme_bw(base_size=22) + 
	theme(axis.text.x=element_text(angle=90,size=16), legend.position='none') + labs(y = 'Residuals') 
	#ggsave(paste0('../Temp5/',gsub("*./.*","",Mutants[m]),'.pdf'), width=16,height=16)
	ggsave(paste0('../Temp5/',gsub("*./.*","",Mutants[m]),'.pdf'), width=25,height=9)
	
#geom_jitter(color='grey50',alpha=0.3,width=0.01, shape = 1, stroke=1) +
})

ggsave('../Temp7/sfxn-alg11-2.pdf', width = 12, height = 10)

dfa <- data_per_animal[data_per_animal$Strain == 'Nsun5-/-',]
dfa$Outlier <- as.factor(ifelse(dfa$MouseID %in% out.df$MouseID, 1, 0)) 
dfa$MouseID <- droplevels(dfa$MouseID)
dfa$Strain <- droplevels(dfa$Strain)
t <- list(
  size = 14,
  color = 'black')
p <- plot_ly(dfa,x=~`stride_length`, y=~`step_width`, z=~`step_length1`, type='scatter3d', color=~`Outlier`,colors=c("black","red")) %>% 
layout(scene = list(
      xaxis = list(title = "Stride Length"),
      yaxis = list(title = "Step Width"),
      zaxis = list(title = "Step Length")
    ), font=t)

p <- plot_ly(dfa,x=~`nose_lateral_displacement`, y=~`tip_tail_lateral_displacement`, z=~`base_tail_lateral_displacement`, type='scatter3d', color=~`Outlier`,colors=c("black","red")) %>% 
layout(scene = list(
      xaxis = list(title = "Nose LD"),
      yaxis = list(title = "Tip Tail LD"),
      zaxis = list(title = "Base Tail LD")
    ), font=t)


ggplot(df,aes(x=Index,y=stride_length)) + geom_line(aes(color=MouseID)) + geom_point(size=0.1)+ 
scale_color_manual(values=c("red","grey50","grey50","grey50","grey50","grey50","grey50","grey50"))

#df.out[df.out$MouseID %in% c('J79719','J82869','J86327'),'Outlier'] <- 1
#df.out[df.out$MouseID %in% c('J67783','J81952','J81953'),'Outlier'] <- 1
#layout_matrix <- rbind(c(1,1,1,2,2),c(1,1,1,3,3))
#p <- gridExtra::grid.arrange(C,CX,CY,layout_matrix=layout_matrix)
#p1|p2

df <- data_per_strain
tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin])
df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, Strain = df[,'Strain'],
	Outlier = as.numeric(!tmp$wfinal01),WeightL = tmp$wloc, WeightS = tmp$wscat)
df.out[df.out$Strain == 'C57BL/6NJ', 'Outlier'] <- -1
df.out$Outlier <- as.factor(df.out$Outlier)
ggplot(df.out, aes(x = Distance1, y=Distance2)) + geom_point(alpha=0.8, aes(color=Outlier), size=11) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(-1,1),as.character(Strain),'')),size=24,box.padding=6) + 
labs(x = 'Distance1', y = 'Distance2') + ggtitle('KOMP strain outliers') + scale_color_manual(values=c("black","grey50","red")) + 
theme_bw(base_size = 60) + theme(legend.position='none')
ggsave('../Temp7/komp-strain-outliers.pdf', width=31.9, height=31.9)

Mutants.out <- c("Cldn13-/-","Kcnd3-/-","Keg1-/-","Rab11a-/+","Rapgef1-/+","Sema6a-/+","Sfxn5-/+")
tmp <- gtools::permutations(n=7,r=2,v=1:7)
lapply(seq(nrow(tmp)), function(x) {cat("Permutation = ", x, "\n");
	Mutants <- Mutants.out[tmp[x,]];
	komp_lda(Mutants) 
}
)

#Mutants <- c('Arpc5l-/+','Fam120b-/+')
#Mutants <- c('Tlk1-/-','Fam120b-/+')
Mutants <- c('Ndufs8-/+','Fam120b-/+')
Mutants <- c('Rusc1-/-','Spin1-/+')
Mutants <- c('Dgkh-/-','Tmem222-/+')
komp_lda <- function(Mutants){
	df <- data.frame()
	CtrlStrain <- "C57BL/6NJ"
	for (m in seq(Mutants)){
		#cat("Mutant", paste0(Mutants[m]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
		df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
		df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
		df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 1)), ]
		df1$Strain <- droplevels(df1$Strain)
		df <- rbind(df,df1)
	}
	#df <- unique(df)
	df <- df[,-which(names(df) %in% c('BodyLength','TestDate','TestAge','Sex'))] #Remove BodyLength
	#df$Outlier <- as.factor(sapply(seq(nrow(df)), function(x) df.out[df.out$MouseID == df$MouseID[x], 'Outlier']))
	#df[df$Strain == CtrlStrain, 'Outlier'] <- 0 
	MutStrain <- setdiff(levels(df$Strain),CtrlStrain)
	df$Strain <- factor(df$Strain, levels = c(CtrlStrain,MutStrain[1],MutStrain[2]), ordered=TRUE)
	#df <- df[,-which(names(df) %in% c('BodyLength'))] #Remove BodyLength
	df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
	#FReffects <- 'BodyLength'
	#formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
   	#fits <- lapply(formulas, lm, data = df)
   	#df_resid <- data.frame(sapply(seq(Phenos.lin), function(x) resid(fits[[x]])))
   	#colnames(df_resid) <- Phenos.lin
   	#df_resid <- cbind(Strain = df$Strain, df_resid)
   	#df_lda <- df_resid
   	#df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])

   	
   	df_resid <- df
	df_svd <- svd(df_resid[,sapply(df_resid,is.numeric)])
	df_pca <- df_svd$u %*% diag(df_svd$d)
	tmp <- df_svd$d^2/(nrow(df)-1)
	df_lda <- data.frame(Strain = df$Strain, df_pca[,1:6])
	#colnames(df_lda)[2:ncol(df_lda)] <- Phenos.lin
	#df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])
	fit_lda <- lda(Strain ~ ., data = df_lda, prior = c(1,1,1)/3)
	lda_values <- predict(fit_lda, df_lda[,-1])
	C <- ggplot(data = data.frame(Strain = df$Strain, lda_values$x), aes(x=LD1,y=LD2,shape=Strain,color=Strain,fill=Strain)) + geom_point(size = 2,aes(color=Strain)) + 
	stat_ellipse(geom = "polygon", alpha = 1/3, aes(fill=Strain)) + theme_bw(base_size = 16) + theme(legend.position = 'top') + 
	scale_color_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) + 
	scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) 
	#ggrepel::geom_text_repel(aes(label=ifelse(df$Outlier ==1,as.character(df$MouseID),'')),size=8,box.padding=2,show.legend=FALSE) 


	#tmp <- plsda(df[,names(df) %in% c(Phenos.lin)], df$Strain)
	#ggplot(data = data.frame(Strain = df$Strain, tmp$variates$X), aes(x=comp1,y=comp2,shape=Strain)) + geom_point(aes(color=Strain)) + stat_ellipse(aes(color = Strain))

	#ggord::ggord(fit_lda,df$Strain,veclsz=NA,labcol=NA) + theme_bw(base_size=16) + theme(legend.position = 'top')
	#lda_result <- data.frame(Strain = df$Strain, lda = predict(fit_lda)$x)
	PCLD_df <- as.data.frame(fit_lda$scaling)
	rownames(PCLD_df) <- Phenos.lin
	PCLD1 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,1]))
	PCLD1$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
		'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
	CX <- ggplot(PCLD1, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = 'Loadings') + ggtitle('LD1')
	PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
	PCLD2$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
		'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
	CY <- ggplot(PCLD2, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = NULL) + ggtitle('LD2')
	#PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
	#acc <- sum(df$Strain == predict(fit_lda)$class)/nrow(df)
	#C <- ggord::ggord(fit_lda,df$Strain,veclsz=NA,labcol=NA) + theme_classic(base_size=16) + theme(legend.position = 'top') 
	#+ ggtitle(paste0('Accuracy: ', round(acc,2)))
	#blankPlot <- ggplot() + geom_blank(aes(1,1)) + theme_void()
	#CC <- gridExtra::grid.arrange(CY,C,blankPlot,CX, ncol=2,nrow=2,widths = c(1,4), heights=c(4,1))
	#p <- plot(fit_lda)
	#p <- plot_grid(C,CX,CY,ncol=3)

	phenotype <- PCLD1$Pheno[which.max(PCLD1$value)]
	phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
	p1 <- ggplot(df_lda, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + 
			theme(legend.position = 'none', axis.text.x=element_blank(), 
				axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + 
			scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a")))


	phenotype <- PCLD2$Pheno[which.max(PCLD2$value)]
	phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
	p2 <- ggplot(df_lda, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + 
			theme(legend.position = 'none', axis.text.x=element_blank(), 
				axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + 
				scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a")))
	

	layout_matrix <- rbind(c(1,1,1,2,3),c(1,1,1,4,5))
	p <- gridExtra::grid.arrange(C,CX,CY,p1,p2,layout_matrix = layout_matrix)
	return(p)
}

#Mutants <- c("Ndufs8-/+","Fam120b-/+")
Mutants <- c('Nsun5-/-','Fam120b-/+')
Mutants <- c('Rusc1-/-','Spin1-/+')
Mutants <- c('Arpc5l-/+','Fam120b-/+')
df <- data.frame()
CtrlStrain <- "C57BL/6NJ"
for (m in seq(Mutants)){
	#cat("Mutant", paste0(Mutants[m]), "\n")
	CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
	df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
	df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
   	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
	df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 1)), ]
	df1$Strain <- droplevels(df1$Strain)
	df <- rbind(df,df1)
}
#df <- unique(df)
df <- df[,-which(names(df) %in% c('BodyLength','TestDate','TestAge','Sex'))] #Remove BodyLength
#df$Outlier <- as.factor(sapply(seq(nrow(df)), function(x) df.out[df.out$MouseID == df$MouseID[x], 'Outlier']))
#df[df$Strain == CtrlStrain, 'Outlier'] <- 0 
MutStrain <- setdiff(levels(df$Strain),CtrlStrain)
df$Strain <- factor(df$Strain, levels = c(CtrlStrain,MutStrain[1],MutStrain[2]), ordered=TRUE)
#df <- df[,-which(names(df) %in% c('BodyLength'))] #Remove BodyLength


#Step 1
df0 <- df[df$Strain %in% CtrlStrain,]
df1 <- df[df$Strain %in% Mutants[1],]
df2 <- df[df$Strain %in% Mutants[2],]

mean_vec <- matrix(0,nrow=length(Phenos.lin),3)
mean_vec[,1] <- apply(df0[,sapply(df0,is.numeric)],2,mean)
mean_vec[,2] <- apply(df1[,sapply(df1,is.numeric)],2,mean)
mean_vec[,3] <- apply(df2[,sapply(df2,is.numeric)],2,mean)

#Computing the scatter matrices using Ledoit-Wolf Shrinkage estimator
cov0 <- nlshrink::linshrink_cov(as.matrix(df0[,-which(names(df0) %in% c('Strain','MouseID','Genotype'))])) 
cov1 <- nlshrink::linshrink_cov(as.matrix(df1[,-which(names(df1) %in% c('Strain','MouseID','Genotype'))]))
cov2 <- nlshrink::linshrink_cov(as.matrix(df2[,-which(names(df2) %in% c('Strain','MouseID','Genotype'))]))

#cov0 <- cov(as.matrix(df0[,-which(names(df0) %in% c('Strain','MouseID','Genotype'))]))
#cov1 <- cov(as.matrix(df1[,-which(names(df1) %in% c('Strain','MouseID','Genotype'))]))
#cov2 <- cov(as.matrix(df2[,-which(names(df2) %in% c('Strain','MouseID','Genotype'))])) 
covw <- cov0 + cov1 + cov2
#Computing Between-class scatter matrix 
overall_mean <- apply(mean_vec,1,mean)
N <- c(nrow(df0),nrow(df1),nrow(df2))
sum <- 0 
for (i in 1:3){
	sum <- sum + N[i]*(mean_vec[,i] - overall_mean)%*%t(mean_vec[,i] - overall_mean)
}
covb <- sum
tmp <- eigen(solve(covw)%*%covb)
tmp2 <- tmp$vectors[,1:2]
Y <- as.matrix(df[,-which(names(df) %in% c('Strain','MouseID','Genotype'))])%*%tmp2
Y <- as.data.frame(Y)
Y$Strain <- df$Strain
C <- ggplot(data = Y, aes(x=V1,y=V2,shape=Strain,color=Strain,fill=Strain)) + geom_point(size = 2,aes(color=Strain)) + 
	stat_ellipse(geom = "polygon", alpha = 1/3, aes(fill=Strain)) + theme_bw(base_size = 16) + theme(legend.position = 'top') + 
	scale_color_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) + 
	scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) + 
	labs(x = 'LD1', y = 'LD2')

PCLD_df <- as.data.frame(tmp2)
rownames(PCLD_df) <- Phenos.lin
PCLD1 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,1]))
PCLD1$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
	'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
CX <- ggplot(PCLD1, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = 'Loadings') + ggtitle('LD1')
PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
PCLD2$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
	'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
CY <- ggplot(PCLD2, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = NULL) + ggtitle('LD2')
phenotype <- PCLD1$Pheno[which.max(PCLD1$value)]
phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
p1 <- ggplot(df, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + 
		theme(legend.position = 'none', axis.text.x=element_blank(), 
			axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + 
		scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a")))


phenotype <- PCLD2$Pheno[which.max(PCLD2$value)]
phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
p2 <- ggplot(df, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + 
		theme(legend.position = 'none', axis.text.x=element_blank(), 
			axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + 
			scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a")))
	

layout_matrix <- rbind(c(1,1,1,2,3),c(1,1,1,4,5))
p <- gridExtra::grid.arrange(C,CX,CY,p1,p2,layout_matrix = layout_matrix)



#sort(setdiff(unique(data_per_animal$Strain), tmp3v[[4]]))
komp_lda_supp <- function(Mutants){
	df <- data.frame()
	CtrlStrain <- "C57BL/6NJ"
	for (m in seq(Mutants)){
		#cat("Mutant", paste0(Mutants[m]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
		df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
		df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
		df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 1)), ]
		df1$Strain <- droplevels(df1$Strain)
		df <- rbind(df,df1)
	}
	df <- unique(df)
	df <- df[complete.cases(df),]
	df$Outlier <- as.factor(sapply(seq(nrow(df)), function(x) df.out[df.out$MouseID == df$MouseID[x], 'Outlier']))
	df[df$Strain == CtrlStrain, 'Outlier'] <- 0 
	MutStrain <- setdiff(levels(df$Strain),CtrlStrain)
	df$Strain <- factor(df$Strain, levels = c(CtrlStrain,MutStrain[1],MutStrain[2]), ordered=TRUE)
	#FReffects <- 'BodyLength'
	#formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
   	#fits <- lapply(formulas, lm, data = df)
   	#df_resid <- data.frame(sapply(seq(Phenos.lin), function(x) resid(fits[[x]])))
   	#colnames(df_resid) <- Phenos.lin
   	#df_resid <- cbind(Strain = df$Strain, df_resid)
   	#df_lda <- df_resid
   	#df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])
   	df_resid <- df
	df_svd <- svd(df_resid[,sapply(df_resid,is.numeric)])
	df_pca <- df_svd$u %*% diag(df_svd$d)
	df_lda <- data.frame(Strain = df$Strain, df_pca[,1:6])
	#colnames(df_lda)[2:ncol(df_lda)] <- Phenos.lin
	#df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])
	fit_lda <- lda(Strain ~ ., data = df_lda)
	lda_values <- predict(fit_lda, df_lda[,-1])
	C <- ggplot(data = data.frame(Strain = df$Strain, lda_values$x), aes(x=LD1,y=LD2,shape=Strain,color=Strain,fill=Strain)) + 
	geom_point(size = 5,aes(color=Strain)) + stat_ellipse(geom = "polygon", alpha = 1/3, aes(fill=Strain)) + 
	theme_bw(base_size = 25) + theme(legend.position = 'none') + 
	scale_color_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) + 
	scale_fill_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) 
	#+ ggrepel::geom_text_repel(aes(label=ifelse(df$Outlier ==1,as.character(df$MouseID),'')),size=8,box.padding=2,show.legend=FALSE)
	return(C)
}

C1 <- komp_lda_supp(c("Adgre4-/-","Zfp422-/-"))
C2 <- komp_lda_supp(c("Xpnpep3-/+","Tomm22-/+"))
C3 <- komp_lda_supp(c("Hap1-/+","Sema6a-/+"))
C4 <- komp_lda_supp(c("Cmtr2-/+","Mrps22-/+"))
C5 <- komp_lda_supp(c("Rapgef1-/+","Elavl1-/+"))
C6 <- komp_lda_supp(c("Ofcc1-/-","Il1a-/-"))
(C1|C2|C3)/(C4|C5|C6)
dev.print(pdf,'../Temp7/LDA-supp-LabM.pdf',width=17.8, height = 11.1) 
#"#e41a1c","#377eb8"

#####Legends#####
Mutants <- c("Ofcc1-/-","Il1a-/-")
df.tmp <- data.frame(value = c(0,1,-1), Strain = c(CtrlStrain,Mutants[1],Mutants[2]))
df.tmp$Strain <- factor(df.tmp$Strain, levels = c(CtrlStrain,Mutants[1],Mutants[2]), ordered=TRUE)
p0 <- ggplot(df.tmp,aes(y=value,x=Strain,color=Strain,shape=Strain,fill=Strain)) + geom_point() + 
scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a")) + theme_bw(base_size=22) +  
theme(legend.position='top') + 
guides(color = guide_legend(override.aes = list(size = 5))) 
legend <- cowplot::get_legend(p0)
grid.newpage()
grid.draw(legend)
dev.print(pdf,'../Temp7/LDA-supp-legend-F.pdf', width=6.26, height=0.66)


df <- data.frame(Strain = data_per_animal$Strain, MouseID = data_per_animal$MouseID, Sex = data_per_animal$Sex)
df <- cbind(df, apply(data_per_animal[,names(data_per_animal) %in% Phenos.lin],2,function(x) (x - mean(x))/sd(x)))
tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin])
df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, 
	Label = paste0(df[,'MouseID']," (", df[,'Strain'], ")"), MouseID = df[,'MouseID'], Strain = df[,'Strain'],
	Outlier = as.numeric(!tmp$wfinal01))

komp_lda_outliers <- function(Mutants){
	df <- data.frame()
	CtrlStrain <- "C57BL/6NJ"
	for (m in seq(Mutants)){
		#cat("Mutant", paste0(Mutants[m]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
		df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
		df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
		df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 1)), ]
		df1$Strain <- droplevels(df1$Strain)
		df <- rbind(df,df1)
	}
	df <- unique(df)
	df$Outlier <- as.factor(sapply(seq(nrow(df)), function(x) df.out[df.out$MouseID == df$MouseID[x], 'Outlier']))
	df[df$Strain == CtrlStrain, 'Outlier'] <- 0 
	#df <- cbind(df,Outlier = as.factor(df.out[df.out$MouseID %in% df$MouseID, 'Outlier']))
	MutStrain <- setdiff(levels(df$Strain),CtrlStrain)
	#df <- df[,-which(names(df) %in% c('BodyLength'))] #Remove BodyLength
	df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
	FReffects <- 'BodyLength'
	formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
   	fits <- lapply(formulas, lm, data = df)
   	df_resid <- data.frame(sapply(seq(Phenos.lin), function(x) resid(fits[[x]])))
   	colnames(df_resid) <- Phenos.lin
   	df_resid <- cbind(Strain = df$Strain, df_resid)
   	#df_lda <- df_resid
    #df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])

	df_svd <- svd(df_resid[,sapply(df_resid,is.numeric)])
	df_pca <- df_svd$u %*% diag(df_svd$d)
	df_lda <- data.frame(Strain = df$Strain, df_pca[,1:2])
	df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])
	fit_lda <- rda::rda(Strain ~ ., data = df_lda)
	lda_values <- predict(fit_lda, df_lda[,-1])
	C <- ggplot(data = data.frame(Strain = df$Strain, lda_values$x), aes(x=LD1,y=LD2,shape=Strain)) + geom_point(size = 2,aes(color=Strain)) + 
	stat_ellipse(geom = "polygon", alpha = 1/3, aes(fill = Strain)) + theme_bw(base_size = 16) + theme(legend.position = 'top') + 
	scale_color_manual(values = c("C57BL/6NJ" = "#F8766D",assign(MutStrain[1],"#619CFF"),assign(MutStrain[2],"#00BA38"))) + 
	scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#619CFF"),assign(MutStrain[2],"#00BA38"))) + 
	ggrepel::geom_text_repel(aes(label=ifelse(df$Outlier ==1,as.character(df$MouseID),'')),size=8,box.padding=2,show.legend=FALSE)


	#tmp <- plsda(df[,names(df) %in% c(Phenos.lin)], df$Strain)
	#ggplot(data = data.frame(Strain = df$Strain, tmp$variates$X), aes(x=comp1,y=comp2,shape=Strain)) + geom_point(aes(color=Strain)) + stat_ellipse(aes(color = Strain))

	#ggord::ggord(fit_lda,df$Strain,veclsz=NA,labcol=NA) + theme_bw(base_size=16) + theme(legend.position = 'top')
	#lda_result <- data.frame(Strain = df$Strain, lda = predict(fit_lda)$x)
	PCLD_df <- as.data.frame(fit_lda$scaling)
	rownames(PCLD_df) <- Phenos.lin
	PCLD1 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,1]))
	PCLD1$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
		'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
	CX <- ggplot(PCLD1, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = 'Loadings') + ggtitle('LD1')
	PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
	PCLD2$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
		'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
	CY <- ggplot(PCLD2, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = NULL) + ggtitle('LD2')
	#PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
	#acc <- sum(df$Strain == predict(fit_lda)$class)/nrow(df)
	#C <- ggord::ggord(fit_lda,df$Strain,veclsz=NA,labcol=NA) + theme_classic(base_size=16) + theme(legend.position = 'top') 
	#+ ggtitle(paste0('Accuracy: ', round(acc,2)))
	#blankPlot <- ggplot() + geom_blank(aes(1,1)) + theme_void()
	#CC <- gridExtra::grid.arrange(CY,C,blankPlot,CX, ncol=2,nrow=2,widths = c(1,4), heights=c(4,1))
	#p <- plot(fit_lda)
	#p <- plot_grid(C,CX,CY,ncol=3)

	#phenotype <- PCLD1$Pheno[which.max(PCLD1$value)]
	#phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
	#p1 <- ggplot(df_lda, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + 
	#		theme(legend.position = 'none', axis.text.x=element_blank(), 
	#			axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + 
	#		scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#619CFF"),assign(MutStrain[2],"#00BA38")))


	#phenotype <- PCLD2$Pheno[which.max(PCLD2$value)]
	#phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
	#p2 <- ggplot(df_lda, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + 
	#		theme(legend.position = 'none', axis.text.x=element_blank(), 
	#			axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + 
	#			scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#619CFF"),assign(MutStrain[2],"#00BA38")))
	

	layout_matrix <- rbind(c(1,2),c(1,3))
	p <- gridExtra::grid.arrange(C,CX,CY,layout_matrix = layout_matrix)
	return(p)
}

pairs <- data.frame(expand.grid(Mutants.out,Mutants.out))
apply(pairs, 1, function(x) if (unname(pairs[x,])[1] == unname(pairs[x,])[2]) {next} else {komp_lda_outliers(as.vector(unname(pairs[x,])))})



p1 <- plot_grid(C,CX,CY,ncol = 3)
p2 <- plot_grid(C,CX,CY,ncol = 3)
p3 <- plot_grid(C,CX,CY, ncol = 3)
p4 <- plot_grid(C,CX,CY, ncol = 3)

plot_grid(p1,p2,p3,p4,nrow = 4)
dev.print(pdf,'../Temp5/lda_new_plot.pdf', width = 16, height = 22)


#Specific to Pcdh and Sfxn5 (mv outlier animal lines)
lda_result$MouseID <- rep(0,30)
lda_result$MouseID[c(1,10,11,22)] <- 1
lda_result$MouseID <- as.factor(lda_result$MouseID) 
C <- C + ggrepel::geom_text_repel(aes(label=ifelse(lda_result$MouseID %in% c(1),as.character(df$MouseID),'')),size=8,box.padding=2)

Mutants.out <- df.out[df.out$Outlier==1,'Strain']
Mutants <- Mutants.out[c(1,6)]
komp_lda(Mutants)[[2]]
dev.print(pdf,paste0('../Temp5/',paste(c(gsub("./.","",Mutants)),collapse="-"),"-lda.pdf"),width=11,height=11)

all_comb <- combn(Mutants.out,2)
acc.mat <- matrix(0,nrow=length(Mutants.out),ncol=length(Mutants.out))
rownames(acc.mat) <- Mutants.out; colnames(acc.mat) <- Mutants.out; diag(acc.mat) <- rep(1,length(Mutants.out))
apply(all_comb,2,function(x) {komp_lda(Mutants = all_comb[,x])[[2]];})
acc <- sapply(seq(ncol(all_comb)), function(x) komp_lda(Mutants=all_comb[,x])[[2]])
acc.mat[lower.tri(acc.mat,diag=FALSE)] <- acc
acc.mat <- forceSymmetric(acc.mat,uplo="L")
col_fun <- circlize::colorRamp2(c(0.70,0.80,0.90,1),c("#e6f598",
"#99d594","#3288bd","#5e4fa2"))

########
Mutants <- Mutants.out[c(12,16)]
df <- data.frame()
for (m in seq(Mutants)){
	#cat("Mutant", paste0(Mutants[m]), "\n")
	CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
	df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
	df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
   	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
	df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 2)), ]
	df1$Strain <- droplevels(df1$Strain)
	df <- rbind(df,df1)
}
df <- unique(df)
df <- df[,-which(names(df) %in% c('BodyLength'))] #Remove BodyLength
df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
df_svd <- svd(df[,sapply(df,is.numeric)])
df_pca <- df_svd$u %*% diag(df_svd$d)
df1 <- df[,Phenos.lin]
rownames(df1) <- paste0(df$MouseID,"(",df$Strain,")") 
kmu.list <- list(); tmp <- numeric(); nclusters <- 3
invisible(lapply(1:choose(62,nclusters), function(c) {kmu.list[[c]] <<- kmeans(df1, centers = df1[sample(nrow(df1),nclusters),], 
nstart = 25);tmp[c] <<- as.numeric(kmu.list[[c]]['tot.withinss']);}))
mykmeans <- kmu.list[[which.min(tmp)]]
factoextra::fviz_cluster(mykmeans, data = df1,
   star.plot = TRUE, ggtheme = theme_minimal(), repel = TRUE, axes = c(1,2),
   palette = c("#999999", "#E69F00", "#56B4E9"),
   geom = 'text',ellipse = TRUE, ellipse.type = 'convex', ellipse.level = 0.80) + 
theme_minimal(base_size = 18)



#"#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026"
Heatmap(as.matrix((acc.mat)), cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = 'left',
	column_names_side = 'bottom',col = col_fun, row_names_gp = gpar(fontsize = 16, fontface="italic"),
	column_names_gp = gpar(fontsize = 16, fontface="italic"), 
	heatmap_legend_param = list(at = c(0.70,0.80,0.90,1),title = "Accuracy", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(8, "cm"), just = c("right", "top")), border=TRUE,
	cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", acc.mat[i, j]), x, y, gp = gpar(fontsize = 15))
})
dev.print(pdf,'../Temp4/lda-accuracy.pdf',width=9,height=9)

#Strains: em1J 
df <- data_per_animal[data_per_animal$Strain %in% c(as.character(df.out[df.out$Outlier==1,'Strain']), 'C57BL/6NJ'),]
df$Strain <- droplevels(df$Strain)


df <- data_per_animal[data_per_animal$Strain %in% c('C57BL/6NJ','C57BL/6NJ-Zfp422<em1J>/J',
	'B6N(Cg)-Coq8b<tm1b(EUCOMM)Hmgu>/J','B6N(Cg)-Tox4<tm1b(KOMP)Mbp>/J',
	'B6N(Cg)-Nsun5<tm2b(EUCOMM)Wtsi>/J','B6N(Cg)-Ptcd3<tm1.1(KOMP)Vlcg>/J'),c('Strain','BG',Phenos.lin)]
df$Strain <- droplevels(df$Strain)
df$BG <- droplevels(df$BG)

df[sample(df$Strain=='C57BL6/NJ',10),]


df <- data_per_animal[data_per_animal$BG %in% c('em1J','C57BL/6NJ','Wtsi'),names(df) %in% c('BG',Phenos.lin)]
df$BG <- droplevels(df$BG) 
df[,names(df) %in% c(Phenos.lin)] <- apply(df[,names(df) %in% c(Phenos.lin)], 2, function(x) (x - mean(x))/sd(x))
df_svd <- svd(df[,sapply(df,is.numeric)])
df_eigvals <- (df_svd$d^2)/(dim(df)[1] - 1)
B <- ggplot(data.frame(evals = df_eigvals/sum(df_eigvals)), aes(y = evals, x = seq(1,9))) + geom_bar(stat="identity") + 
geom_line(linetype = 'solid', color = 'black', size = 1) + 
geom_text(label = paste0(round(df_eigvals/sum(df_eigvals),2), '%'),vjust = -.2, hjust = 0) + 
labs(x = 'Dimension', y = 'Percent of Variance') + theme_bw() + scale_x_discrete(limits=seq(1,9))
df_eigvecs <- df_svd$v
df_pca <- df_svd$u %*% diag(df_svd$d)
df_lda <- data.frame(Strain = df$Strain, df_pca)
fit_lda <- lda(Strain ~ ., data = df_lda)
fit_lda <- rrcov::Linda(BG ~ ., data = df_lda, method='fsa')
lda_result <- data.frame(Strain = df$Strain, lda = predict(fit_lda)$x)

C <- ggplot(lda_result, aes(lda.LD1,lda.LD2, color=Strain)) + geom_point(alpha=0.5,size=3) + labs(x='LD1',y='LD2')

PCLD_df <- as.data.frame(df_eigvecs %*% fit_lda$scaling)
rownames(PCLD_df) <- Phenos.lin
PCLD1 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,1]))
CX <- ggplot(PCLD1, aes(x = reorder(Phenos, value), y = value)) + geom_bar(stat = 'identity', color = 'black') + 
theme(axis.text.x = element_text(angle=45,vjust=.5)) + labs(x = 'Phenotypes', y = 'Loadings')
PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
CY <- ggplot(PCLD2, aes(x = reorder(Phenos, value), y = value)) + geom_bar(stat = 'identity', color = 'black') + 
theme(axis.text.y = element_text(angle=45)) + labs(x = 'Phenotypes', y = 'Loadings') + 
coord_flip()

#element_text(angle = 45, vjust = 0.5, hjust=1) #element_blank
blankPlot <- ggplot() + geom_blank(aes(1,1)) + theme_void()
CC <- gridExtra::grid.arrange(CY,C,blankPlot,CX, ncol=2,nrow=2,widths = c(1,4), heights=c(4,1))

cov_df <- as.data.frame(cov(df[,sapply(df, is.numeric)])) 
cov_df$Phenos <- Phenos.lin.Nomen
cov_melt <- reshape::melt(cov_df, id = "Phenos")  
cov_melt$Phenos <- factor(cov_melt$Phenos,levels = unique(cov_melt$Phenos))



#Analysis at the individual animal level
#1
Mutants <- unique(data_per_animal$Strain)
df <- data_per_animal[data_per_animal$Strain %in% Mutants[1],
names(data_per_animal) %in% c('MouseID','Sex','speed','step_length1','step_width','stride_length')]
tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])

df.out <- data.frame(rodist = tmp$md.rob, mahadist = tmp$md.cla, MouseID = df[,'MouseID'], 
	Outlier = (ifelse(mvoutlier::aq.plot(df[,names(df) %in% Phenos.lin])$outliers==TRUE,1,0)))
df.out$Outlier <- as.factor(df.out$Outlier)
ggplot(df.out, aes(x = mahadist, y=rodist)) + geom_point(alpha=0.8, aes(color=Outlier), size=4) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1,-1),as.character(MouseID),'')),size=6,box.padding=2) + 
labs(x = 'Mahalanobis Distance', y = 'Robust Distance') + ggtitle('KOMP Outliers') + scale_color_manual(values = c('grey50','red'))
theme_bw(base_size = 16) + theme(legend.position='none')


df <- df[,-1]
df <- cbind(id = 1:dim(df)[1], df)
df.melt <- reshape::melt(df[,-2], id.vars = 'id')
df.melt <- cbind(df.melt, Strain = rep(df$MouseID, 4), Outlier = as.factor(rep(df.out$Outlier, 4)))
ggplot(df.melt, aes(y = value, color = Outlier, x = id)) + geom_jitter(alpha=0.7, size = 3) + 
facet_wrap(~variable, scales='free') + scale_color_manual(values = c('grey50','red')) + theme_bw(base_size = 16) + 
theme(legend.position = 'none') + labs(x = 'Index', y='Phenotype') 


#2
Mutants <- unique(data_per_animal$Strain)
df <- data_per_animal[data_per_animal$Strain %in% Mutants[1],names(data_per_animal) %in% c('MouseID',Phenos.lin)]
df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/(sd(x)))





#3
Mutants <- unique(data_per_animal$Strain)
df <- data_per_animal[data_per_animal$Strain %in% Mutants[1], names(data_per_animal) %in% c('MouseID',Phenos.lin)]
X <- as.matrix(df)
lambda <- 1/sqrt(max(dim(X)[1],dim(X)[2]))
i <- 1
tau <- 0.5
eps <- 0.2
while((tau > 0.3) & (eps > 0.1)){
	X.svd <- svd(X)
	d <- thresh.l1(X.svd$d,lambda)
	L <- X.svd$u %*% diag(d) %*% X.svd$v
	S <- X - L 
	eps <- sum(abs(as.numeric(X - L - S)))
	X <- L + S
	lambda <- lambda - 0.01
	tau <- sum(S!=0)/(dim(S)[1]*dim(S)[2])
	cat("Iteration = ", i + 1, "\n")
}

Mutants <- unique(data_per_animal$Strain)
tmp <- list()
lapply(seq(length(Mutants)), function(m) {
	df <- data_per_animal[data_per_animal$Strain %in% Mutants[m], names(data_per_animal) %in% c('MouseID',Phenos.lin)];
	X <- as.matrix(df[,-1]);
	rpcaX <- rpca::rpca(X,max.iter=30000,lambda=0.3);
	od1 <- X - rpcaX$S;
	od2 <- rpcaX$L - rpcaX$S;
	tmp[[m]] <- sort(apply(od1,1,function(x) sum(x^2))/apply(od2,1,function(x) sum(x^2)));
})
 
#4
Mutants <- sort(unique(data_per_animal$Strain)) 
m <- 59
df <- data_per_animal[data_per_animal$Strain %in% 'C57BL/6NJ', names(data_per_animal) %in% c('MouseID',Phenos.lin)];
tmp <- rospca::robpca(as.matrix(df[,-1]))

df.plot <- data.frame(sd = tmp$sd, od = tmp$od, Outlier = as.factor(1 - as.numeric(tmp$flag.all)))

ggplot(df.plot, aes(x = sd, y = od, color = Outlier)) + geom_point(size = 2) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1),as.character(df$MouseID),'')),size=6,box.padding=2) + 
geom_hline(yintercept = tmp$cutoff.od) + geom_vline(xintercept = tmp$cutoff.sd) + 
scale_color_manual(values = c('black','red')) + theme(legend.position = 'none') + labs(x = 'Score distance',
	y='Orthogonal distance')
df.plot <- list()
df.cutoff <- list()
Mutants <- setdiff(Mutants, 'C57BL/6NJ')
lapply(seq(length(Mutants)), function(m){
	df <- data_per_animal[data_per_animal$Strain %in% Mutants[m], names(data_per_animal) %in% c('MouseID',Phenos.lin)];
	tmp <- rospca::robpca(as.matrix(df[,-1]))

	df.plot[[m]] <<- data.frame(Strain = rep(paste0(Mutants[m]),nrow(df)), MouseID = df$MouseID, 
		sd = tmp$sd, od = tmp$od,Outlier = as.factor(1 - as.numeric(tmp$flag.all)))
	df.cutoff[[m]] <<- data.frame(Strain = Mutants[m], sdcut = rep(tmp$cutoff.sd, nrow(df)), odcut = rep(tmp$cutoff.od, nrow(df)))
	#ggplot(df.plot, aes(x = sd, y = od, color = Outlier)) + geom_point(size = 2) + 
	#ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1),as.character(df$MouseID),'')),size=6,box.padding=2) + 
	#geom_hline(yintercept = tmp$cutoff.od) + geom_vline(xintercept = tmp$cutoff.sd) + theme_bw(base_size=16) + 
	#scale_color_manual(values = c('black','red')) + theme(legend.position = 'none') + labs(x = 'Score distance',
	#	y='Orthogonal distance') + ggtitle(paste0(Mutants[m]))

})


df <- do.call(rbind,df.plot)
dfc <- do.call(rbind, df.cutoff)

ggplot(df, aes(x = sd, y = od, color = Outlier)) + geom_point(size = 1) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1),as.character(MouseID),'')),size=3,box.padding=1) + 
theme_bw(base_size=18) + geom_hline(data=dfc, aes(yintercept = odcut), linetype='dashed') + geom_vline(data=dfc, aes(xintercept = sdcut), linetype='dashed') + 
scale_color_manual(values = c('black','red')) + theme(legend.position = 'none') + labs(x = 'Score distance',
	y='Orthogonal distance') + facet_wrap(.~Strain)

ggsave('../Temp2/allstrains-out2.pdf', width=18, height=22)

df <- data_per_animal[data_per_animal$Strain %in% 'C57BL/6NJ', names(data_per_animal) %in% c('MouseID',Phenos.lin)];
tmp <- rospca::robpca(as.matrix(df[,-1]))

df$Outlier <- as.factor(1 - as.numeric(tmp$flag.all))

df.plot <- data.frame(sd = tmp$sd, od = tmp$od, Outlier = as.factor(1 - as.numeric(tmp$flag.all)))

ggplot(df.plot, aes(x = sd, y = od, color = Outlier)) + geom_point(size = 2) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1),as.character(df$MouseID),'')),size=6,box.padding=2) + theme_bw(base_size = 18) + 
geom_hline(yintercept = tmp$cutoff.od, linetype='dashed') + geom_vline(xintercept = tmp$cutoff.sd, linetype='dashed') + 
scale_color_manual(values = c('black','red')) + theme(legend.position = 'none') + labs(x = 'Score distance',
	y='Orthogonal distance') + ggtitle('C57BL/6NJ')
ggsave('../Temp2/ctrlstrain-out.pdf', width=8, height=8)

T <- (as.matrix(df[,-1]) - as.matrix(rep(1,nrow(df)))%*%t(as.matrix(tmp$center)))%*%tmp$loadings
od.out <- df[which((1 - as.numeric(tmp$flag.od))==1),'MouseID']
sd.out <- df[which((1 - as.numeric(tmp$flag.sd))==1),'MouseID']
tmp2 <- data.frame(MouseID = df$MouseID, PC1 = T[,1], PC2 = T[,2], 
	Outlier.all = as.factor(1 - as.numeric(tmp$flag.all)), Outlier.sd = as.factor(1 - as.numeric(tmp$flag.sd)),
	Outlier.x = as.factor(ifelse(df$MouseID %in% setdiff(sd.out,od.out),1,0)))
ggplot(tmp2, aes(x = PC1, y = PC2)) + geom_point(size = 3,aes(color = Outlier.all),alpha=0.7) + 
stat_ellipse(level=0.989, type='t', linetype=2) + 
ggrepel::geom_text_repel(aes(label=ifelse(Outlier.sd %in% c(1),as.character(df$MouseID),'')),size=6,box.padding=2) + theme_bw(base_size = 18) + 
scale_color_manual(values = c('black','red')) + theme(legend.position='none') + ggtitle('C57BL/6NJ')
ggsave('../Temp2/ctrlstrain-pc.pdf',width=8,height=8)

Mutants <- sort(unique(data_per_animal$Strain)) 
df.plot <- list()
df.cutoff <- list()
lapply(seq(length(Mutants)), function(m){
	df <- data_per_animal[data_per_animal$Strain %in% Mutants[m], names(data_per_animal) %in% c('MouseID',Phenos.lin)];
	tmp <- rospca::robpca(as.matrix(df[,-1]))

	df.plot[[m]] <<- data.frame(Strain = rep(paste0(Mutants[m]),nrow(df)), MouseID = df$MouseID, 
		sd = tmp$sd, od = tmp$od,Outlier = as.factor(1 - as.numeric(tmp$flag.all)))
	df.cutoff[[m]] <<- data.frame(Strain = Mutants[m], sdcut = rep(tmp$cutoff.sd, nrow(df)), odcut = rep(tmp$cutoff.od, nrow(df)))
	#ggplot(df.plot, aes(x = sd, y = od, color = Outlier)) + geom_point(size = 2) + 
	#ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(1),as.character(df$MouseID),'')),size=6,box.padding=2) + 
	#geom_hline(yintercept = tmp$cutoff.od) + geom_vline(xintercept = tmp$cutoff.sd) + theme_bw(base_size=16) + 
	#scale_color_manual(values = c('black','red')) + theme(legend.position = 'none') + labs(x = 'Score distance',
	#	y='Orthogonal distance') + ggtitle(paste0(Mutants[m]))

})


df <- do.call(rbind,df.plot)
data_per_animal$Outlier <- df$Outlier

ggplot(data_per_animal, aes(x = Strain, y = step_width)) + geom_boxplot() + geom_point(aes(color=Outlier), alpha = 0.7) + 
scale_color_manual(values = c('black','red')) + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position = 'none') + 
labs(x = 'Strain')



df.ctrl <- data_per_animal[data_per_animal$Strain == 'C57BL/6NJ',names(data_per_animal) %in% c(Phenos.lin)]
colnames(df.ctrl) <- Phenos.lin.Nomen
df.ctrl <- cbind(id = as.factor(1:nrow(df.ctrl)), df.ctrl)
df.melt <- reshape::melt(df.ctrl, id.vars = 'id')
df.melt <- cbind(Strain = rep('C57BL/6NJ', 1674), MouseID = rep(df[df$Strain %in% 'C57BL/6NJ','MouseID'], length(Phenos.lin)),
	df.melt, Outlier = rep(df[df$Strain %in% 'C57BL/6NJ','Outlier'], length(Phenos.lin)))
ggplot(df.melt, aes(x = value, y = Strain)) + geom_boxplot() + geom_point(aes(color=Outlier)) + theme_bw(base_size=20) + 
scale_color_manual(values = c('black','red')) + facet_wrap(.~variable, scales='free',ncol= 1, strip.position="right") + 
theme(legend.position='none',axis.text.y = element_blank(), strip.text.y = element_text(size = 9)) + 
labs(y = NULL, x = NULL) + ggtitle('C57BL/6NJ')
ggsave('../Temp2/ctrl-out-box.pdf', height=12,width=8)


df.ctrl <- data_per_animal[data_per_animal$Strain == 'C57BL/6NJ',names(data_per_animal) %in% c(Phenos.lin,'Outlier')]
t <- list(
  size = 14,
  color = 'black')
p <- plot_ly(df,x=~`stride_length`, y=~`step_width`, z=~`step_length1`, type='scatter3d', color=~`Outlier`,colors=c("black","red")) %>% 
layout(scene = list(
      xaxis = list(title = "Stride Length"),
      yaxis = list(title = "Step Width"),
      zaxis = list(title = "Step Length")
    ), font=t)

p <- plot_ly(df,x=~`nose_lateral_displacement`, y=~`tip_tail_lateral_displacement`, z=~`base_tail_lateral_displacement`, type='scatter3d', color=~`Outlier`,colors=c("black","red")) %>% 
layout(scene = list(
      xaxis = list(title = "Nose LD"),
      yaxis = list(title = "Tip Tail LD"),
      zaxis = list(title = "Base Tail LD")
    ), font=t)


plotly_IMAGE(p, width = 500, height = 500, format = "pdf", scale = 2,  out_file = "test.pdf")





