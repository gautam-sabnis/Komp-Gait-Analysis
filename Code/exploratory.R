
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
Strain8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 8]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strain8,]

#Remove outliers 
invisible(sapply(seq(length(unique(data_per_animal$Strain))), function(s) Map(function(p) {
	vals <- data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][[Phenos.lin[p]]];
	outliers <- boxplot.stats(vals)$out
	ids <- match(outliers, vals) 
	data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][paste0(Phenos.lin[p])][ids,] <<- NA

}, seq(length(Phenos.lin)))))

data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin)], by = data_per_animal[c("Strain")], 
	FUN = function(x) c(mn = mean(x, na.rm=TRUE), md = median(x,na.rm=TRUE), std = sd(x,na.rm=TRUE)))

komp_desc_plot <- function(phenotype, phenoname, align){
	data_per_animal$sorted_strain = make_sort_column(data_per_animal, 'Strain', phenotype ,FUN = mean)
	me <- data_per_strain[data_per_strain$Strain == 'C57BL/6NJ', names(data_per_strain) %in% c(phenotype)][1]
	std <- data_per_strain[data_per_strain$Strain == 'C57BL/6NJ', names(data_per_strain) %in% c(phenotype)][3]
	df <- as.data.frame(data_per_strain[,phenotype]) 
	df <- cbind(Strain = data_per_strain$Strain, df)

	if(align == 'h'){
		ggplot(data_per_animal, aes_string(x = "sorted_strain", y = phenotype)) + geom_point(alpha = 0.5) + 
		geom_point(data = df, aes(x = Strain, y = df[,"mn"]), color = 'red', size = 2.5) + 
		geom_errorbar(data = df, aes(x = Strain, df[,"mn"], ymin = df[,"mn"]-df[,"std"],
		ymax = df[,"mn"]+df[,"std"]), color = 'grey') + geom_hline(yintercept = me - std, 
		color = "red", size = 1,alpha = 0.3) +  geom_hline(yintercept = me, color = "red", size = 1,alpha = 0.5) +  
		geom_hline(yintercept = me + std, color = "red", size = 1,alpha = 0.3) + labs(x = 'Strain', y = paste0(phenoname)) + 
		theme_bw(base_size = 16) + ggtitle('') + theme(legend.position = 'none') + 
        theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6))
    } else {
		ggplot(data_per_animal, aes_string(y = "sorted_strain", x = phenotype)) + geom_point(alpha = 0.5) + 
		geom_point(data = data_per_strain, aes(y = Strain, x = df[,"mn"]), color = 'red', size = 2.5) + 
		geom_errorbarh(data = data_per_strain, aes(y = Strain, df[,"mn"], xmin = df[,"mn"]-df[,"std"],
		xmax = df[,"mn"]+df[,"std"]), color = 'grey') + geom_vline(xintercept = me - std, 
		color = "red", size = 1,alpha = 0.3) +  geom_vline(xintercept = me, color = "red", size = 1,alpha = 0.5) +  
		geom_vline(xintercept = me + std, color = "red", size = 1,alpha = 0.3) + labs(y = 'Strain', x = paste0(phenoname)) +  
		theme_bw(base_size = 16) + ggtitle('') + theme(legend.position = 'none') + 
           theme(axis.text.y = element_text(hjust = 1, size=5))
	}
}

p1 <- komp_desc_plot('step_length1', 'Step Length', align='h')
p2 <- komp_desc_plot('stride_length', 'Stride Length', align='h')
p3 <- komp_desc_plot('step_width', 'Step Width', align='h')
p4 <- komp_desc_plot('base_tail_lateral_displacement', 'Base Tail LD', align='h')
p5 <- komp_desc_plot('tip_tail_lateral_displacement', 'Tip Tail LD', align='h')
p6 <- komp_desc_plot('nose_lateral_displacement', 'Nose LD', align='h')
p7 <- komp_desc_plot('limb_duty_factor', 'Limb Duty', align='h')
p8 <- komp_desc_plot('temporal_symmetry', 'Temporal Symmetry', align='h')

plot_grid(p1,p2,p3,nrow=3)
ggsave(snakemake@output[[1]], width=9, height=12)
plot_grid(p4,p5,p6,nrow=3)
ggsave(snakemake@output[[2]], width=9, height=12)
plot_grid(p7,p8,nrow=2)
ggsave(snakemake@output[[3]], width=9, height=12)

