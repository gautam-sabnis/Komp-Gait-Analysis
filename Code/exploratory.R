#setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)

Phenos.lin <- c("angular_velocity","speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry","base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Angular Velocity","Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
"TS","Base Tail LD","Tip Tail LD","Nose LD")

data_per_stride <- read.delim(snakemake@input[[1]], stringsAsFactors = FALSE)
names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Genotype'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))
levels(data_per_stride$Strain)[1] <- "C57BL/6NJ"
levels(data_per_stride$Strain)[119] <- "Mrps22<tm1.1(KOMP)Vlcg> -/+"

#Remove Strains
toMatch <- c("Esrrb", "<em2J>/J COIN", "IMPC","Rik","RIK")
matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
Strains <- setdiff(unique(data_per_stride$Strain), matches)
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% 
vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= 
		summary(as.numeric(table(data_per_stride$MouseID)))[2]]),]
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin,'BodyLength')], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
data_per_animal <- cbind(Strain, TestAge, TestDate, Sex, data_per_animal)


#Remove outliers 
invisible(sapply(seq(length(unique(data_per_animal$Strain))), function(s) Map(function(p) {
	vals <- data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][[Phenos.lin[p]]];
	outliers <- boxplot.stats(vals)$out
	ids <- match(outliers, vals) 
	data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][paste0(Phenos.lin[p])][ids,] <<- NA

}, seq(length(Phenos.lin)))))

data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("Strain")], 
	FUN = function(x) c(mn = mean(x, na.rm=TRUE), md = median(x,na.rm=TRUE), std = sd(x,na.rm=TRUE)))

#<em1J> - em1J, (KOMP) - Mbp, Wtsi, Vlcg, (EUCOMM) - Hmgu 
data_per_animal$BG <- as.factor(ifelse(grepl("Mbp", data_per_animal$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_animal$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_animal$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_animal$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_animal$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))

komp_desc_plot <- function(df_animal=data_per_animal, df_strain = data_per_strain, phenotype, phenoname, align){
	df_animal$sorted_strain = make_sort_column(df_animal, 'Strain', phenotype ,FUN = mean, decreasing=FALSE)
	me <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenotype)][1]
	std <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenotype)][3]
	df <- as.data.frame(df_strain[,phenotype]) 
	df <- cbind(Strain = df_strain$Strain, df)

	if(align == 'h'){
		ggplot(df_animal, aes_string(x = "sorted_strain", y = phenotype)) + geom_point(alpha = 0.9, size = 3, aes(color=Sex)) + 
		geom_point(data = df, aes(x = Strain, y = df[,"mn"]), color = 'black', size = 4, alpha = 0.5) + 
		scale_color_manual(values=c("#E41A1C", "#377EB8")) + 
		geom_errorbar(data = df, aes(x = Strain, df[,"mn"], ymin = df[,"mn"]-df[,"std"],
		ymax = df[,"mn"]+df[,"std"]), color = 'grey') + geom_hline(yintercept = me - std, 
		color = "black", size = 1) +  geom_hline(yintercept = me, color = "black", size = 1) +  
		geom_hline(yintercept = me + std, color = "black", size = 1) + labs(x = 'Strain', y = paste0(phenoname)) + 
		theme_bw(base_size = 16) + ggtitle('') + theme_bw(base_size = 18) + theme(legend.position = c(0.05,0.85), 
			legend.text = element_text(size=13),legend.background = element_blank(),
			axis.text.x = element_text(angle = 85, vjust = 0.55, size = 12), legend.key.size = unit(0.1, "cm")) 
        
    } else {
		ggplot(df_animal, aes_string(y = "sorted_strain", x = phenotype)) + geom_point(alpha = 0.9, size = 3, aes(color=Sex)) + 
		geom_point(data = df, aes(y = Strain, x = df[,"mn"]), color = 'black', size = 4, alpha = 0.5) + 
		scale_color_manual(values=c("#E41A1C", "#377EB8")) + 
		geom_errorbarh(data = df, aes(y = Strain, df[,"mn"], xmin = df[,"mn"]-df[,"std"],
		xmax = df[,"mn"]+df[,"std"]), color = 'grey') + geom_vline(xintercept = me - std, 
		color = "black", size = 1,alpha = 0.6) +  geom_vline(xintercept = me, color = "black", size = 1,alpha = 0.5) +  
		geom_vline(xintercept = me + std, color = "black", size = 1,alpha = 0.6) + labs(y = 'Strain', x = paste0(phenoname)) +  
		theme_bw(base_size = 22) + ggtitle('') + theme(legend.position = c(0.80,0.05), 
			legend.text = element_text(size=14),legend.background = element_blank(),
			axis.text.y = element_text(size = 16,face="italic"), legend.key.size = unit(0.1, "cm"))
	}
}

invisible(lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal, 
	data_per_strain, Phenos.lin[x], Phenos.lin.Nomen[x], align='v'), inherits=TRUE)))

plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[1]], width=9, height=12)
plot_grid(p3,p4+labs(y=NULL),p5+labs(y=NULL),ncol=3)
ggsave(snakemake@output[[2]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[3]], width=9, height=12)

#dev.print(pdf,'../Temp3/sl-sw-strl-desc.pdf',width=15,height=20)

#Time based plots 
data_per_strain_BG <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("BG","TestDate")], 
	FUN = function(x) c(mn = mean(x, na.rm=TRUE), md = median(x,na.rm=TRUE), std = sd(x,na.rm=TRUE)))

time_plot <- function(df_strain,phenotype, phenoname, addPoints=FALSE){
	df <- as.data.frame(df_strain[,phenotype]) 
	df <- cbind(BG = df_strain$BG, TestDate = df_strain$TestDate, df)

	if (addPoints == FALSE){
		ggplot(df, aes(x = as.Date(TestDate), y = df[,"md"], color = BG)) + geom_point(size=3) +
      geom_line(size=1) + scale_x_date(date_breaks = "months" , date_labels = "%b-%y") + labs(x = 'TestDate', y=paste0(phenoname))
	} else {
		ggplot(data_per_animal, aes_string(x = as.Date(TestDate), y = phenotype)) + geom_jitter(alpha = 0.5,width = 0.50, aes(color=BG)) + 
		geom_point(data = df, aes(x = as.Date(TestDate), y = df[,"mn"]), size = 0.5) + 
		geom_line(data = df, aes(x = as.Date(TestDate), y = df[,"mn"], color=BG), size = 1) + 
		scale_x_date(date_breaks = "months" , date_labels = "%b-%y") + labs(x = 'TestDate', y=paste0(phenoname))}
}


invisible(lapply(seq((Phenos.lin)), function(x) assign(paste0("p",x), time_plot(data_per_strain_BG, Phenos.lin[x],
	Phenos.lin.Nomen[x],addPoints=TRUE), inherits=TRUE)))
legend <- get_legend(p1)
p <- plot_grid(p1+labs(x=NULL)+theme(legend.position='none'),p2+labs(x=NULL)+theme(legend.position='none'),p3+labs(x=NULL)+theme(legend.position='none'),
	p4+labs(x=NULL)+theme(legend.position='none'),p5+labs(x=NULL)+theme(legend.position='none'),
	p6+labs(x=NULL)+theme(legend.position='none'),p7+labs(x=NULL)+theme(legend.position='none'),
	p8+labs(x=NULL)+theme(legend.position='none'),p9+labs(x=NULL)+theme(legend.position='none'),
	p10+theme(legend.position='none'),nrow=10)
plot_grid(p,legend,rel_widths = c(3, .4))

ggsave(snakemake@output[[4]], width=9, height=12)


komp_time_plot <- function(df_animal=data_per_animal, df_strain=data_per_strain, phenotype, phenoname,colr){
	me <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenoname)][1]
	std <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenoname)][3]
	df <- as.data.frame(df_strain[,phenoname]) 
	df <- cbind(Strain = df_strain$Strain, df)

	df_animal$TestDate <- as.Date(df_animal$TestDate, format = "%Y-%m-%d")
	df_animal$Genotype <- ifelse(df_animal$Strain=='C57BL/6NJ', 'Control','Mutant')

	if (colr=='Sex'){
		ggplot(df_animal, aes_string(x = "TestDate", y = phenoname, color = 'Sex')) + geom_jitter(size = 2, alpha = 0.8) + 
		scale_color_manual(values=c("#E41A1C", "#377EB8")) +  geom_smooth(method='loess') + theme_bw(base_size=18) + 
		theme(axis.text.x = element_text(angle=45,hjust=1)) + scale_x_date(breaks=unique(df_animal$TestDate)) + 
		geom_hline(yintercept = me - std, color = "black", size = 1, alpha = 0.6) +  
		geom_hline(yintercept = me, color = "black", size = 1, alpha = 0.6) +  
		geom_hline(yintercept = me + std, color = "black", size = 1, alpha = 0.6) + labs(x = 'TestDate', y = paste0(phenotype))
	} else {
		ggplot(df_animal, aes_string(x = "TestDate", y = phenoname, color = 'Genotype')) + geom_jitter(size = 2, alpha = 0.8) + 
		scale_color_manual(values=c("#6a51a3", "#d94801")) +  geom_smooth(method='loess') + theme_bw(base_size=24) + 
		theme(axis.text.x = element_text(angle=45,hjust=1)) + scale_x_date(breaks=unique(df_animal$TestDate)) + 
		geom_hline(yintercept = me - std, color = "black", size = 1, alpha = 0.6) +  
		geom_hline(yintercept = me, color = "black", size = 1, alpha = 0.6) +  
		geom_hline(yintercept = me + std, color = "black", size = 1, alpha = 0.6) + labs(x = 'TestDate', y = paste0(phenotype))

	}
	
}
invisible(lapply(seq((Phenos.lin)), function(x) assign(paste0("p",x), komp_time_plot(data_per_animal, data_per_strain,
	Phenos.lin.Nomen[x],Phenos.lin[x],colr='Sex'), inherits=TRUE)))
legend <- get_legend(p1)
p <- plot_grid(p3+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank(),axis.ticks.x=element_blank()),
	p4+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank(),axis.ticks.x=element_blank()),
	p5+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank()), nrow=3)
plot_grid(p,legend,rel_widths = c(3, .4))
ggsave('../Temp5/time-sl-sw-stl.pdf', width=16.50,height=6)

p <- plot_grid(p1+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank(),axis.ticks.x=element_blank()),
	p2+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank(),axis.ticks.x=element_blank()),
	p6+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank()),
	p7+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank()),
	p8+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank()),
	p9+labs(x=NULL)+theme(legend.position='none',axis.text.x=element_blank()), nrow=6)

ggsave('../Temp5/time-supp.pdf', width=12.50,height=12)

#Density based plots
#By Strain
data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin,"BodyLength")], by = data_per_animal[c("Strain")], FUN = mean)
data_per_strain$BG <- as.factor(ifelse(grepl("Mbp", data_per_strain$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_strain$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_strain$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_strain$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_strain$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))

df <- data_per_strain[-1, Phenos.lin] #Remove C57BL/6NJ
names(df) <- Phenos.lin.Nomen
df <- cbind(id = 1:dim(df)[1], df)
df.melt <- reshape::melt(df, id.vars = 'id')
df.melt <- cbind(df.melt, BG = rep(data_per_strain$BG[-1], each = length(Phenos.lin)))
df.ctrl <- data_per_animal[data_per_animal$Strain=='C57BL/6NJ',Phenos.lin]
names(df.ctrl) <- Phenos.lin.Nomen
df.ctrl <- cbind(id = 1:dim(df.ctrl)[1], df.ctrl)
df.ctrl.melt <- reshape::melt(df.ctrl, id.vars = 'id')
df.ctrl.melt <- cbind(df.ctrl.melt, BG = rep('C57BL/6NJ', dim(df.ctrl)[2]-1))
df.melt <- rbind(df.melt,df.ctrl.melt)

ggplot(df.melt, aes(x = value, color = BG)) + geom_density(lwd=1) + 
facet_wrap(~ variable, scales = 'free') + 
scale_color_brewer(palette="Dark2") + labs(y = NULL, x = NULL) + theme_bw(base_size = 18)
ggsave(snakemake@output[[5]], width=12, height=12)


ggplot(df.melt, aes(y = value, fill = BG)) + geom_boxplot(alpha=0.5) + 
facet_wrap(~ variable, scales = 'free') + 
scale_color_brewer(palette="Dark2") + labs(y = NULL, x = NULL) + theme_bw(base_size = 18)
ggsave(snakemake@output[[6]], width=12, height=12)

#By animal
df <- data_per_animal[,Phenos.lin]
names(df) <- Phenos.lin.Nomen
df <- cbind(id = 1:dim(df)[1], df)
df.melt <- reshape::melt(df, id.vars = 'id')
df.melt <- cbind(df.melt, BG = rep(data_per_animal$BG, each = length(Phenos.lin)))

ggplot(df.melt, aes(x = value, color = BG)) + geom_density() + 
facet_wrap(~ variable, scales = 'free') + 
scale_color_brewer(palette="Dark2") + labs(y = NULL, x = NULL) + theme_bw(base_size = 18)
ggsave(snakemake@output[[7]], width=12, height=12)


ggplot(df.melt, aes(y = value, fill = BG)) + geom_boxplot(alpha=0.5) + 
facet_wrap(~ variable, scales = 'free') + 
scale_color_brewer(palette="Dark2") + labs(y = NULL, x = NULL) + theme_bw(base_size = 18)
ggsave(snakemake@output[[8]], width=12, height=12)

#Stride density plot
data_per_stride <- data_per_stride[data_per_stride$MouseID %in% data_per_animal$MouseID,]
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
df.tmp <- data.frame(Strides = as.numeric(table(data_per_stride$MouseID)))
Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) 
	data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
df.tmp <- cbind(df.tmp, Sex = Sex)
ggplot(df.tmp, aes(x = Strides, color = Sex)) + 
geom_density(lwd=1.1) + ggtitle('All speed bins with ang_vel = (-20,20)') + 
geom_vline(xintercept = as.numeric(summary(as.numeric(table(data_per_stride$MouseID)))[2]), linetype = 'dashed',
color='black',lwd=1) + labs(x = 'Number of Strides', y = 'Density') + theme_bw(base_size=14) + 
scale_color_manual(values=c("#E41A1C", "#377EB8"))
ggsave('../Temp7/central-angvel-bin-stride-density2.pdf', height=4,width=5)


#Pick atleast x animals 
tmp <- (table(as.numeric(table(data_per_animal$Strain))))
df.tmp <- data.frame(x = seq(length(tmp)), y = rev(cumsum(rev(tmp))))
df.tmp$highlight <- ifelse(df.tmp$x == 5, 'yes', 'no')
ggplot(df.tmp, aes(x,y, fill=highlight)) + geom_bar(stat='identity',color='black') + labs(x = 'Number of Animals', 
	y = 'Number of Mutant Lines') + scale_x_discrete(limits=rownames(tmp)) + 
geom_text(label = paste0(df.tmp$y), vjust = -.25, hjust = .5, size=3) + theme_bw(base_size=12) + 
scale_fill_manual(values = c( "yes"="red", "no"="grey50"), guide=FALSE)
ggsave('../Temp7/atleastx-hist2.pdf', height=4, width=5)

#Exploring batch effects 
df <- data_per_animal[data_per_animal$Strain == 'C57BL/6NJ',]
df$Strain <- droplevels(df$Strain)
ggplot(df, aes(x = TestDate, y=speed)) + geom_boxplot() + geom_point(aes(color=Sex), size=2,alpha=0.8) + 
scale_color_manual(values=c("#E41A1C", "#377EB8"))

tmp <- sapply(seq(length(Phenos.lin)), function(x) {fit <- lm(df[[Phenos.lin[x]]] ~ BodyLength + speed + Sex + TestDate, 
	data = df); (Anova(fit,type='II')$'Sum Sq'/sum(Anova(fit,type='II')$'Sum Sq')*100)[4]})



#####Legends#####
df.tmp <- data.frame(value = c(0,1,-1), Genotype = c('Mutant','Mutant (Outlier)','Control'))
p0 <- ggplot(df.tmp,aes(y=value,x=Genotype,color=Genotype)) + geom_point(alpha=0.4) + 
scale_color_manual(values=c("black","#6a51a3","#d94801")) + theme_bw(base_size=18) +  
theme(legend.position='top') + 
guides(color = guide_legend(override.aes = list(size = 5))) 
legend <- cowplot::get_legend(p0)
grid.newpage()
grid.draw(legend)
dev.print(pdf,'../Temp7/C-legend.pdf', width=6, height=1.5)
#theme(legend.title=element_blank(),legend.position='top')