#setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#data_per_stride <- read.delim('../Data/kompdf', stringsAsFactors = TRUE)

Phenos.lin <- c("angular_velocity","speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry","base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Angular Velocity","Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
"Temporal Symmetry","Base Tail LD","Tip Tail LD","Nose LD")

data_per_stride <- read.delim(snakemake@input[[1]], stringsAsFactors = FALSE)
names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Strain.Name'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
	'speed_25_ang_vel_neg20'),]
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
data_per_animal <- cbind(Strain, data_per_animal)

#Pick Strains with atleast 8 animals 
#Strain8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 8]
#data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strain8,]

#Remove outliers 
invisible(sapply(seq(length(unique(data_per_animal$Strain))), function(s) Map(function(p) {
	vals <- data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][[Phenos.lin[p]]];
	outliers <- boxplot.stats(vals)$out
	ids <- match(outliers, vals) 
	data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][paste0(Phenos.lin[p])][ids,] <<- NA

}, seq(length(Phenos.lin)))))

data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("Strain")], 
	FUN = function(x) c(mn = mean(x, na.rm=TRUE), md = median(x,na.rm=TRUE), std = sd(x,na.rm=TRUE)))

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

lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal, 
	data_per_strain, Phenos.lin[x], Phenos.lin.Nomen[x], align='h'), inherits=TRUE))

plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[1]], width=9, height=12)
plot_grid(p4,p5,p6,nrow=3)
ggsave(snakemake@output[[2]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[3]], width=9, height=12)


#Genetic backgrounds: <em1J>,(KOMP),(EUCOMM)
explore_strain_family <- function(family = c("em1J","KOMP","EUCOMM")){
	toMatch <- c("C57BL/6NJ",family)
	matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
	data_per_stride <- data_per_stride[data_per_stride$Strain %in% matches,]
	data_per_stride$Strain <- as.factor(data_per_stride$Strain)
	data_per_stride$Strain <- droplevels(data_per_stride$Strain)
	data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
	Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
	data_per_animal <- cbind(Strain, data_per_animal)
	data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("Strain")], 
	FUN = function(x) c(mn = mean(x, na.rm=TRUE), md = median(x,na.rm=TRUE), std = sd(x,na.rm=TRUE)))
	data_per_strain$Strain <- gsub('.*-(.*)<.*','\\1',data_per_strain$Strain)
	data_per_animal$Strain <- gsub('.*-(.*)<.*','\\1',data_per_animal$Strain)
	return(list(data_per_animal,data_per_strain)) 
}

BGs <- c("em1J","KOMP","EUCOMM")
invisible(lapply(seq(length(BGs)), function(x) assign(paste0("data_per_strain.",BGs[x]), explore_strain_family(BGs[x])[[2]],
	inherits=TRUE)))
invisible(lapply(seq(length(BGs)), function(x) assign(paste0("data_per_animal.",BGs[x]), explore_strain_family(BGs[x])[[1]],
	inherits=TRUE)))

#em1J
lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal.em1J, 
	data_per_strain.em1J, Phenos.lin[x], Phenos.lin.Nomen[x], align='h'), inherits=TRUE))
plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[4]], width=9, height=12)
plot_grid(p4,p5,p6,nrow=3)
ggsave(snakemake@output[[5]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[6]], width=9, height=12)

#KOMP
lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal.KOMP, 
	data_per_strain.KOMP, Phenos.lin[x], Phenos.lin.Nomen[x], align='h'), inherits=TRUE))
plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[7]], width=9, height=12)
plot_grid(p4,p5,p6,nrow=3)
ggsave(snakemake@output[[8]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[9]], width=9, height=12)

#EUCOMM
lapply(seq(length(Phenos.lin)), function(x) assign(paste0("p",x), komp_desc_plot(data_per_animal.EUCOMM, 
	data_per_strain.EUCOMM, Phenos.lin[x], Phenos.lin.Nomen[x], align='h'), inherits=TRUE))
plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[10]], width=9, height=12)
plot_grid(p4,p5,p6,nrow=3)
ggsave(snakemake@output[[11]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[12]], width=9, height=12)


