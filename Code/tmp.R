#####
data_per_animal$Genotype <- ifelse(data_per_animal$Strain=='C57BL/6NJ','Control','Mutant') 
df2 <- df[, Phenos.lin]
names(df2) <- Phenos.lin.Nomen
df2 <- cbind(id = 1:dim(df2)[1], df2)
df.melt <- reshape::melt(df2, id.vars = 'id')
df.melt <- cbind(df.melt, MouseID = rep(names(table(df$MouseID)), as.numeric(table(df$MouseID))*length(Phenos.lin)), 
	BodyLength = rep(data_per_animal[data_per_animal$MouseID %in% unique(df$MouseID), 'BodyLength'], 
		as.numeric(table(df$MouseID))*length(Phenos.lin)),
	Genotype = rep(data_per_animal[data_per_animal$MouseID %in% unique(df$MouseID), 'Genotype'], 
		as.numeric(table(df$MouseID))*length(Phenos.lin)))



tmp2 <- data_per_stride[data_per_stride$Strain %in% c('C57BL/6NJ',"C57BL/6NJ-Ccdc30<em1J>/J"),]



fit <- lmer(value ~ variable:(Genotype + BodyLength - 1) + (1|MouseID), data = df.melt,
	REML = FALSE,control = lmerControl(optimizer ='optimx', 
    				optCtrl=list(method='nlminb')))



df_animal <- aggregate(x = df[,names(df) %in% c(Phenos.lin,'BodyLength')], by = df[c('MouseID')], FUN = mean)
Strain <- sapply(seq(dim(df_animal)[1]), function(x) df[df$MouseID == df_animal$MouseID[x], 'Strain'][1])
TestAge <- sapply(seq(dim(df_animal)[1]), function(x) df[df$MouseID == df_animal$MouseID[x], 'TestAge'][1])
df_animal <- cbind(Strain,TestAge,df_animal)
df_animal$Genotype <- ifelse(df_animal$Strain=='C57BL/6NJ','Control','Mutant') 


fit <- lm(step_width ~ BodyLength + Genotype, data = df_animal)
fit2 <- lm(step_width ~ BodyLength + Genotype + TestAge, data = df_animal)
 S(lm(
 	~ BodyLength + Strain, data = df_animal))



TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
data_per_animal <- cbind(data_per_animal,TestAge)
data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin,"BodyLength")], by = data_per_animal[c("Strain")], FUN = mean)
TestAge <- sapply(seq(dim(data_per_strain)[1]), function(x) data_per_animal[data_per_animal$MouseID == data_per_strain$MouseID[x], 'TestAge'][1])
data_per_strain <- cbind(TestAge, data_per_strain)
data_per_strain$BG <- as.factor(ifelse(grepl("Mbp", data_per_strain$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_strain$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_strain$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_strain$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_strain$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))
mod <- lm(cbind(speed,limb_duty_factor,step_length1,step_width,stride_length,temporal_symmetry,
	base_tail_lateral_displacement,tip_tail_lateral_displacement,nose_lateral_displacement) ~ BG, data = data_per_strain)



sample.df <- function(df, n) df[sample(nrow(df), n), , drop = FALSE]


########
df <- data_per_animal[data_per_animal$BG %in% c('Mbp','C57BL/6NJ'),]
Ctrl <- 'C57BL/6NJ'
Mutant <- 'Cep126<tm1b(KOMP)Mbp> -/-'
df1 <- df[df$Strain %in% c(Ctrl,Mutant),]
df1$Strain <- droplevels(df1$Strain)

df2 <- rbind(df1[df1$Strain == Mutant,], sample.df(subset(df1, Strain == Ctrl), nrow(df1[df1$Strain == Mutant,])))
df2[,sapply(df2,is.numeric)] <- apply(df2[,sapply(df2,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
df2$Genotype <- ifelse(df2$Strain==Ctrl,'Control','Mutant')
df3 <- df2[,Phenos.lin]
rownames(df3) <- paste0(df2$MouseID,"(",df2$Genotype,")") 
kmu.list <- list(); tmp <- numeric(); nclusters <- 2
invisible(lapply(1:choose(62,nclusters), function(c) {kmu.list[[c]] <<- kmeans(df3, centers = df3[sample(nrow(df3),nclusters),], 
    nstart = 25);tmp[c] <<- as.numeric(kmu.list[[c]]['tot.withinss']);}))
mykmeans <- kmu.list[[which.min(tmp)]]
fviz_cluster(mykmeans, data = df3, palette = c("#999999", "#E69F00"),
    star.plot = TRUE, ggtheme = theme_minimal(), repel = TRUE, axes = c(1,2),
    geom = 'text',ellipse = TRUE, ellipse.type = 'convex', ellipse.level = 0.80) + 
theme_minimal(base_size = 18)






#postprocess 08/06/2020

hits <- sapply(seq(nrow(pvalGenSignif2)), function(x) sum(pvalGenSignif2[x,] != " "))
vignettes <- c(Mutants[which(hits==4)],Mutants[which(hits==3)])

tmp <- lapply(seq(length(Phenos.lin)), function(x) names(sort(abs(esizeGen2[,x]), decreasing=TRUE)[1:5])) 

vignette_plot <- function(phenotype, phenoname, esizeGen2 = tmp[[3]]){
	esize <- names(sort(abs(esizeGen2[,paste0(phenotype)]), decreasing=TRUE)[1:5])
	dflist <- list();
	invisible(lapply(seq(length(esize)), function(x) {CtrlIDs <- unique(subset(controlids.df,Strain == esize[x])$MouseID);
	df <<- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('Strain','Sex','TestDate',paste0(phenoname))];
	df['Genotype'] <<- ifelse(df$Strain == "C57BL/6NJ", 'Control','Mutant');
	df$Genotype <<- relevel(factor(df$Genotype), ref = "Control");
	df$Strain <<- droplevels(df$Strain);
	dflist[[x]] <<- df}))
	df <- do.call(rbind, dflist)
	df$Strain <- rep(esize, sapply(seq(length(esize)), function(x) dim(dflist[[x]])[1]))
	p <- ggplot(df, aes_string(y=phenoname, x='Genotype')) + geom_boxplot(outlier.shape=NA,width=0.3) + 
	geom_jitter(aes(color=Genotype), alpha=0.8, width=0.1) + scale_color_manual(values=c("#6a51a3", "#d94801")) + 
	theme_bw(base_size = 24) + theme(legend.position='none') + labs(y = paste0(phenotype)) + facet_grid(~Strain)
	return(p)
}







df <- data_per_animal[-which(data_per_animal$Strain == 'C57BL/6NJ'),]
ggplot(df,aes(Strain)) + geom_bar() + theme(axis.text.x = element_text(angle=90, hjust=1),axis.text.x.top = element_text(vjust = 0.5))