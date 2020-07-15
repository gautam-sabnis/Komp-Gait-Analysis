#setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#data_per_stride <- read.delim('../Data/kompdf', stringsAsFactors = FALSE)

Phenos.lin <- c("angular_velocity","speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry","base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Angular Velocity","Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
"TS","Base Tail LD","Tip Tail LD","Nose LD")

data_per_stride <- read.delim(snakemake@input[[1]], stringsAsFactors = FALSE)
names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Strain.Name'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'

#Remove Strains
toMatch <- c("B6.Cg-Esrrb<tm1(cre)Yba>/J", "<em2J>/J COIN")
matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
Strains <- setdiff(unique(data_per_stride$Strain), matches)
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
	'speed_25_ang_vel_neg20'),]
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
data_per_animal <- cbind(Strain, TestDate, data_per_animal)

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

komp_desc_plot <- function(df_animal, df_strain, phenotype, phenoname, align){
	df_animal$sorted_strain = make_sort_column(df_animal, 'Strain', phenotype ,FUN = mean)
	me <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenotype)][1]
	std <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenotype)][3]
	df <- as.data.frame(df_strain[,phenotype]) 
	df <- cbind(Strain = df_strain$Strain, df)

	if(align == 'h'){
		ggplot(df_animal, aes_string(x = "sorted_strain", y = phenotype)) + geom_point(alpha = 0.5) + 
		geom_point(data = df, aes(x = Strain, y = df[,"mn"]), color = 'red', size = 2.5) + 
		geom_errorbar(data = df, aes(x = Strain, df[,"mn"], ymin = df[,"mn"]-df[,"std"],
		ymax = df[,"mn"]+df[,"std"]), color = 'grey') + geom_hline(yintercept = me - std, 
		color = "red", size = 1,alpha = 0.3) +  geom_hline(yintercept = me, color = "red", size = 1,alpha = 0.5) +  
		geom_hline(yintercept = me + std, color = "red", size = 1,alpha = 0.3) + labs(x = 'Strain', y = paste0(phenoname)) + 
		theme_bw(base_size = 16) + ggtitle('') + theme(legend.position = 'none') + 
        theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6))
    } else {
		ggplot(df_animal, aes_string(y = "sorted_strain", x = phenotype)) + geom_point(alpha = 0.5) + 
		geom_point(data = df, aes(y = Strain, x = df[,"mn"]), color = 'red', size = 2.5) + 
		geom_errorbarh(data = df, aes(y = Strain, df[,"mn"], xmin = df[,"mn"]-df[,"std"],
		xmax = df[,"mn"]+df[,"std"]), color = 'grey') + geom_vline(xintercept = me - std, 
		color = "red", size = 1,alpha = 0.3) +  geom_vline(xintercept = me, color = "red", size = 1,alpha = 0.5) +  
		geom_vline(xintercept = me + std, color = "red", size = 1,alpha = 0.3) + labs(y = 'Strain', x = paste0(phenoname)) +  
		theme_bw(base_size = 16) + ggtitle('') + theme(legend.position = 'none') + 
           theme(axis.text.y = element_text(hjust = 1, size=5))
	}
}

invisible(lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal, 
	data_per_strain, Phenos.lin[x], Phenos.lin.Nomen[x], align='h'), inherits=TRUE)))

plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[1]], width=9, height=12)
plot_grid(p4,p5,p6,nrow=3)
ggsave(snakemake@output[[2]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[3]], width=9, height=12)

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
		ggplot(data_per_animal, aes_string(x = as.Date(TestDate), y = phenotype)) + geom_jitter(alpha = 0.5,width = 0.50) + 
		geom_point(data = df, aes(x = as.Date(TestDate), y = df[,"mn"], color=BG), size = 4) + 
		geom_line(data = df, aes(x = as.Date(TestDate), y = df[,"mn"], color=BG), size = 1) + 
		scale_x_date(date_breaks = "months" , date_labels = "%b-%y")}
}


invisible(lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), time_plot(data_per_strain_BG, Phenos.lin[x],
	Phenos.lin.Nomen[x]), inherits=TRUE)))
legend <- get_legend(p1)
p <- plot_grid(p1+labs(x=NULL)+theme(legend.position='none'),p2+labs(x=NULL)+theme(legend.position='none'),p3+labs(x=NULL)+theme(legend.position='none'),
	p4+labs(x=NULL)+theme(legend.position='none'),p5+labs(x=NULL)+theme(legend.position='none'),
	p6+labs(x=NULL)+theme(legend.position='none'),p7+labs(x=NULL)+theme(legend.position='none'),
	p8+labs(x=NULL)+theme(legend.position='none'),p9+labs(x=NULL)+theme(legend.position='none'),
	p10+theme(legend.position='none'),nrow=10)
plot_grid(p,legend,rel_widths = c(3, .4))

ggsave(snakemake@output[[4]], width=9, height=12)