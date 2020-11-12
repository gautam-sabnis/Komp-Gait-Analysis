Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry", "base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
	"TS","Base Tail LD","Tip Tail LD","Nose LD")
Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
"nose_lateral_displacement_phase")

#rm(list=ls())
#setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)


names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Genotype'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))
levels(data_per_stride$Strain)[1] <- "C57BL/6NJ"
levels(data_per_stride$Strain)[119] <- "Mrps22<tm1.1(KOMP)Vlcg> -/+"
levels(data_per_stride$Strain)[3] <- "Rik<em1J> -/-" 
#Remove Strains
toMatch <- c("Esrrb", "<em2J>/J COIN", "IMPC")
matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
Strains <- setdiff(unique(data_per_stride$Strain), matches)
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% 
vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= 
		summary(as.numeric(table(data_per_stride$MouseID)))[2]]),]
data_per_stride$Strain <- vapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]), FUN.VALUE="1")
data_per_stride$Strain <- vapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]), FUN.VALUE=character(1))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
BodyLength <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'BodyLength'][1])
Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
data_per_animal <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, data_per_animal)

#Filter Strains for which at least 8 animals were tested 
Strains8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 5]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strains8,] 
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains8,]
data_per_animal$Strain <- droplevels(data_per_animal$Strain)
data_per_stride$Strain <- droplevels(data_per_stride$Strain)

#Remove univariate outliers 
invisible(sapply(seq(length(unique(data_per_animal$Strain))), function(s) Map(function(p) {
	vals <- data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][[Phenos.lin[p]]];
	outliers <- boxplot.stats(vals)$out
	ids <- match(outliers, vals) 
	data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][paste0(Phenos.lin[p])][ids,] <<- NA

}, seq(length(Phenos.lin)))))

data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("Strain")], 
	FUN = function(x) c(mn = mean(x, na.rm=TRUE), md = median(x,na.rm=TRUE), std = sd(x,na.rm=TRUE)))

#Descriptive Plot
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
		theme_bw(base_size = 26) + ggtitle('') + theme(legend.position = c(0.60,0.05), 
			legend.text = element_text(size=14),legend.background = element_blank(),
			axis.text.y = element_text(size = 15,face="italic"), legend.key.size = unit(0.1, "cm"))
	}
}

invisible(lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal, 
	data_per_strain, Phenos.lin[x], Phenos.lin.Nomen[x], align='v'), inherits=TRUE)))
legend <- get_legend(p1)

p <- plot_grid(p3+theme(legend.position='none'),p4+labs(y=NULL)+theme(legend.position='none'),
	p5+labs(y=NULL)+theme(legend.position='none'),ncol=3)
plot_grid(p,legend,rel_widths = c(10, .1))

dev.print(pdf,'../Temp5/sl-sw-strl-desc.pdf',width=15,height=25)

p <- plot_grid(p1+theme(legend.position='none'),p2+labs(y=NULL)+theme(legend.position='none'),
	p6+labs(y=NULL)+theme(legend.position='none'),ncol=3)
plot(p)
dev.print(pdf,'../Temp5/sp-ldf-ts-desc.pdf',width=15,height=25)

p <- plot_grid(p7+theme(legend.position='none'),p8+labs(y=NULL)+theme(legend.position='none'),
	p9+labs(y=NULL)+theme(legend.position='none'),ncol=3)
plot(p)
dev.print(pdf,'../Temp5/bt-tt-nose-desc.pdf',width=15,height=25)