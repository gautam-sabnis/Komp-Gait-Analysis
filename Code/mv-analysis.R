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
names(data_per_stride)[names(data_per_stride) == 'OFA_Strain.Name'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
names(data_per_stride)[names(data_per_stride) == "OFA_Arena.ID"] <- 'Arena'
names(data_per_stride)[names(data_per_stride) == "OFA_Experimenter.ID"] <- 'Experimenter'

data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))

#Remove Strains
toMatch <- c("Esrrb", "<em2J>/J COIN", "IMPC")
matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
Strains <- setdiff(unique(data_per_stride$Strain), matches)
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
	'speed_25_ang_vel_neg20'),]
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin,'BodyLength')], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
Arena <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Arena'][1])
Experimenter <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Experimenter'][1])
data_per_animal <- cbind(Strain, TestDate, TestAge, data_per_animal)

data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.lin,"BodyLength")], by = data_per_animal[c("Strain")], FUN = mean)

#Remove multivariate outliers 
gBG <- c('C57BL/6NJ') #'em1J','Hmgu','Mbp','Vlcg','Wtsi'
invisible(lapply(seq(gBG), function(x) {
	df <- data_per_animal[data_per_animal$BG %in% c(gBG),]
	tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])
	ids <- df[tmp$outliers,'MouseID']
	data_per_animal <<- data_per_animal[-which((data_per_animal$BG == gBG) & (data_per_animal$MouseID %in% ids)),]
	data_per_stride <<- data_per_stride[-which((data_per_stride$BG == gBG) & (data_per_stride$MouseID %in% ids)),]
}))

#<em1J> - em1J, (KOMP) - Mbp, Wtsi, Vlcg, (EUCOMM) - Hmgu 
data_per_stride$BG <- as.factor(ifelse(grepl("Mbp", data_per_stride$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_stride$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_stride$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_stride$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_stride$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))

data_per_animal$BG <- as.factor(ifelse(grepl("Mbp", data_per_animal$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_animal$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_animal$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_animal$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("em2J", data_per_animal$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_animal$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ")))))))

data_per_strain$BG <- as.factor(ifelse(grepl("Mbp", data_per_strain$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_strain$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_strain$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_strain$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_strain$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))
tmp <- which(data_per_strain$BG == 'C57BL/6NJ')
data_per_strain[tmp[2],'BG'] <- "em1J"
data_per_strain[tmp[3],'BG'] <- "em1J"

lapply(Phenos.lin, function(x) ggplot(data_per_strain, aes(x = x, color=BG)) + geom_density())

data_per_strain <- data_per_strain[-which(data_per_strain$BG == 'C57BL/6NJ'),]
tmp <- data.frame(matrix(c(apply(data_per_animal[data_per_animal$Strain=='C57BL/6NJ' & data_per_animal$Sex=='Male',sapply(data_per_animal,is.numeric)],2,mean),
apply(data_per_animal[data_per_animal$Strain=='C57BL/6NJ' & data_per_animal$Sex=='Female',sapply(data_per_animal,is.numeric)],2,mean)),nrow=2))
names(tmp) <- names(data_per_strain)[2:11]
tmp <- cbind(Strain = c('C57BL/6NJ','C57BL/6NJ'),tmp)
tmp <- cbind(tmp,BG = c('C57BL/6NJ','C57BL/6NJ'))
data_per_strain <- rbind(data_per_strain, tmp)

#Test the effect of Arena, Experimenter
mvlm <- lm(cbind(speed,limb_duty_factor,step_length1,step_width,stride_length,temporal_symmetry,
	base_tail_lateral_displacement,tip_tail_lateral_displacement,nose_lateral_displacement) ~ Arena*Experimenter, data = data_per_animal)

#Test the effect of gBG
mvlm <- lm(cbind(speed,limb_duty_factor,step_length1,step_width,stride_length,temporal_symmetry,
	base_tail_lateral_displacement,tip_tail_lateral_displacement,nose_lateral_displacement) ~ BG, 
data = data_per_strain)

model1 <- MANOVA.RM::MANOVA.wide(cbind(stride_length,step_width,step_length1) ~ BG, data = data_per_strain)

formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", "BG + (1|Strain)"), simplify=TRUE)) 
fits <- sapply(formulas, function(x) lmer(formula=x, data = data_per_animal, REML = FALSE,
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))),simplify=FALSE)
sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['BG','Pr(>F)'])

tmp <- data.frame(Phenotype = Phenos.lin.Nomen, p = sapply(seq_along(fits), 
	function(x) anova(unname(fits[x])[[1]],type='II')['BG','Pr(>F)']))
p0 <- ggplot(tmp, aes(y = Phenotype, x = p)) + geom_point(size = 5) + geom_vline(xintercept = 0.05, linetype = 'dashed',
	color = 'red', lwd = 1.1) + theme_bw(base_size = 22) + labs(x = 'p-value') + 
scale_x_continuous(breaks = c(0.05,0.25,0.50,0.75),labels = c(0.05,0.25,0.50,0.75))

lapply(seq(Phenos.lin), function(x) {

	avg <- c(as.numeric(summary(fits[[x]])$coefficients[,1])[1], as.numeric(summary(fits[[x]])$coefficients[,1])[1] + 
	as.numeric(summary(fits[[x]])$coefficients[,1])[-1]);
	lwr <- as.numeric(c((summary(fits[[x]])$coefficients)[1,1] - (summary(fits[[x]])$coefficients)[1,2],
	avg[-1] - (summary(fits[[x]])$coefficients)[-1,2]));
	upr <- as.numeric(c((summary(fits[[x]])$coefficients)[1,1] + (summary(fits[[x]])$coefficients)[1,2],
	avg[-1] + (summary(fits[[x]])$coefficients)[-1,2]));
	BG <- c('C57BL/6NJ','em1J','Hmgu','Mbp','Vlcg','Wtsi');

	df.effect <- data.frame(BG = BG, avg = avg, lwr = lwr, upr = upr);
	assign(paste0("p",x), ggplot(df.effect, aes(x = BG, y = avg, ymin=lwr, ymax=upr)) + geom_point(size = 5) + geom_errorbar() + 
	labs(y = paste0(Phenos.lin.Nomen[x])) + theme_bw(base_size = 22) + 
	theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1)), inherits = TRUE)
	#ggsave(paste0('../Temp6/',Phenos.lin[x],'-BG-esize.pdf'), width = 4, height = 4)

})
p0 | (p1/p2) | (p3/p4) | (p5/p6) | (p7/p8) | p9

(p0/p1) | (p2/p3) | (p4/p5) | (p6/p7) | (p8/p9)
dev.print(pdf,'../Temp6/BG-pvalue-esize.pdf', width = 22, height = 9)

lapply(seq(Phenos.lin), function(x) {
	assign(paste0("p",x), redres::plot_resqq(fits[[x]]) + ggtitle(paste0(Phenos.lin.Nomen[x])) + 
		theme_bw(base_size = 22), inherits = TRUE)
	#ggsave(paste0('../Temp6/',Phenos.lin[x],'-BG-qq.pdf'), width = 4, height = 4);
})
#plot_grid(p3,p4,p5,p7,p8,p9,nrow = 2, ncol = 3)
#dev.print(pdf,'../Temp6/6qq.pdf', width = 10, height = 7)
(p1) | (p2/p3) | (p4/p5) | (p6/p7) | (p8/p9)
dev.print(pdf,'../Temp6/BG-lmm-assumption-qq.pdf', width = 22, height = 9)

lapply(seq(Phenos.lin), function(x) {
	bmat <- as.data.frame(lme4::ranef(model))

  	# converting each random effect vector into one line with nest
  	renest <- tidyr::nest(bmat, data = c("grp", "condval", "condsd"))

  	# generating list of ggplot objects for each random effect vector
  plots <- purrr::pmap(list(renest$data, renest$grpvar, renest$term),
                          function(a,b,c){
                            ggplot(data = a, aes_string(sample = "condval")) +
                              qqplotr::stat_qq_band(bandType = "pointwise",
                                                    distribution = "norm",
                                                    fill = "#FBB4AE", alpha = 0.4) +
                              qqplotr::stat_qq_line(distribution = "norm", colour = "#FBB4AE") +
                              qqplotr::stat_qq_point(distribution = "norm") +
                              xlab("Normal quantiles") + theme_bw(base_size = 22) +
                              ylab(paste(b,": ", c)) + ggtitle(paste0(Phenos.lin.Nomen[x]))
                          }
                          )
                          
  
	assign(paste0("p",x), plots[[1]], inherits = TRUE)
	#ggsave(paste0('../Temp6/',Phenos.lin[x],'-BG-ranef.pdf'), width = 4, height = 4);
})
plot_grid(p3,p4,p5,p7,p8,p9,nrow = 2, ncol = 3)
dev.print(pdf,'../Temp6/6ranef.pdf', width = 10, height = 7)

(p1) | (p2/p3) | (p4/p5) | (p6/p7) | (p8/p9)
dev.print(pdf,'../Temp6/BG-lmm-assumption-ranef.pdf', width = 22, height = 9)

heplot(mvlm, variables = c(4,5), xlab = 'Step Width', ylab = 'Stride Length',cex=0.75)
dev.print(pdf,'../Temp5/heplot-sl-sw.pdf',height=4.5,width=4.5)
heplot(mvlm, variables = c(4,3), xlab = 'Step Width', ylab = 'Step Length',cex=0.75)
dev.print(pdf,'../Temp5/heplot-sl1-sw.pdf',height=4.5,width=4.5)
heplot(mvlm, variables = c(1,2), xlab = 'Speed', ylab = 'Limb Duty Factor',cex=0.75)
dev.print(pdf,'../Temp5/heplot-sp-ldf.pdf',height=4.5,width=4.5)
heplot(mvlm, variables = c(6,7), xlab = 'Limb Duty Factor', ylab = 'Base Tail LD',cex=0.75)
dev.print(pdf,'../Temp5/heplot-bt-ldf.pdf',height=4.5,width=4.5)
heplot(mvlm, variables = c(8,9), xlab = 'Tip Tail LD', ylab = 'Nose LD',cex=0.75)
dev.print(pdf,'../Temp5/heplot-nose-tt.pdf',height=4.5,width=4.5)

lapply(seq(9), function(x) {tmp <-  c(summary(mvlm)[[x]]$coefficients[1,'Estimate'], 
	summary(mvlm)[[x]]$coefficients[1,'Estimate'] + summary(mvlm)[[x]]$coefficients[-1,'Estimate']); 
names(tmp)[1] <- "C57BL/6NJ"; tmp[order(as.numeric(tmp))]})

df <- data_per_animal[data_per_animal$BG == 'em1J',]
mvlm <- lm(cbind(speed,limb_duty_factor,step_length1,step_width,stride_length,temporal_symmetry,
	base_tail_lateral_displacement,tip_tail_lateral_displacement,nose_lateral_displacement) ~ Strain, data = df)
heplot(mvlm, variables = c(4,3), xlab = 'Step Width', ylab = 'Step Length')





formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", "BG + (1|Strain)"), simplify=TRUE)) 
fits <- sapply(formulas, function(x) lmer(formula=x, data = data_per_animal, REML = FALSE,control = lmerControl(optimizer ='optimx', 
    optCtrl=list(method='nlminb'))),simplify=FALSE2pvalGen <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['BG','Pr(>F)']) 

lapply(seq(9), function(x) {tmp <- c(summary(fits[[x]])$coefficients[1,'Estimate'], 
	summary(fits[[x]])$coefficients[1,'Estimate'] + summary(fits[[x]])$coefficients[-1,'Estimate']); 
names(tmp)[1] <- "C57BL/6NJ"; tmp[order(as.numeric(tmp))]})

formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", "Arena + Experimenter + (1|Strain)"), simplify=TRUE)) 
fits <- sapply(formulas, function(x) lmer(formula=x, data = data_per_animal, REML = FALSE,control = lmerControl(optimizer ='optimx', 
    optCtrl=list(method='nlminb'))),simplify=FALSE)
pvalGen.Arena <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Arena','Pr(>F)']) 
pvalGen.Experimenter <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Experimenter','Pr(>F)'])












#####
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
	'speed_25_ang_vel_neg20'),]
data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) > 1]),]
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
data_per_animal.mean <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin,'BodyLength')], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal.mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal.mean$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
data_per_animal.mean <- cbind(Strain, TestDate, data_per_animal.mean)
data_per_animal.var <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = var)
data_per_animal <- cbind(data_per_animal.mean,data_per_animal.var) 
data_per_animal <- data_per_animal[,-14]

#Mutants <- c('Cep126<tm1b(KOMP)Mbp> -/-', 'Rmnd5b<em1J> -/-',
#	'Arpc5l<tm1.1(KOMP)Vlcg> -/+','Mpv17<em1J> -/-','Mast3<tm1.1tat_ell(KOMP)Wtsi> -/-','Kcnd3<em1J> -/-',
#	'Rev3l<tm1b(EUCOMM)Hmgu> -/+','Prkcsh<tm1.1(KOMP)Vlcg> -/+')

#Mutants <- c(as.character(Mutants.out[1]),as.character(data_per_strain$Strain[sample(nrow(data_per_strain),1)]))

#Mutants <- c('Cep126<tm1b(KOMP)Mbp> -/-','Arpc5l<tm1.1(KOMP)Vlcg> -/+','Kcnd3<em1J> -/-',
#	'Mast3<tm1.1(KOMP)Wtsi> -/-','Mpv17<em1J> -/-','Mast3<tm1.1(KOMP)Wtsi> -/-')
Mutants.out <- df.out[df.out$Outlier==1,'Strain']
Mutants <- Mutants.out[c(1,12)]
df <- data.frame()
for (m in seq(Mutants)){
	cat("Mutant", paste0(Mutants[m]), "\n")
	CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
	df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
	df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
	df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 2)), ]
	df1$Strain <- droplevels(df1$Strain)
	df <- rbind(df,df1)
}

#df <- data_per_animal[data_per_animal$Strain %in% c(as.character(df.out[df.out$Outlier==1,'Strain']), 'C57BL/6NJ'),]
#df$Strain <- droplevels(df$Strain)

df <- df[,-6]
#tmp <- df[df$Strain == 'C57BL/6NJ',]
#tmp <- tmp[sample(nrow(tmp),17),]
#df <- rbind(df[df$Strain %in% Mutants, ], tmp)
df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
df_svd <- svd(df[,sapply(df,is.numeric)])
df_pca <- df_svd$u %*% diag(df_svd$d)
df_lda <- data.frame(Strain = df$Strain, df_pca)
fit_lda <- lda(Strain ~ ., data = df_lda)
lda_result <- data.frame(Strain = df$Strain, lda = predict(fit_lda)$x)
#C <- ggplot(lda_result, aes(x = lda.LD1, y = lda.LD2)) + stat_ellipse(geom = "polygon", alpha = 0.3, aes(color = Strain, fill = Strain)) + 
#geom_point(aes(color=Strain), alpha=0.8, size = 4) + labs(x='LD1', y='LD2') + 
#theme_classic(base_size=22) + theme(legend.position='top') 
PCLD_df <- as.data.frame(df_svd$v %*% fit_lda$scaling)
rownames(PCLD_df) <- Phenos.lin
PCLD1 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,1]))
CX <- ggplot(PCLD1, aes(x = reorder(Phenos, value), y = value)) + geom_bar(stat = 'identity', color = 'black') + 
theme(axis.text.x = element_text(angle=45,vjust=.5)) + labs(x = 'Phenotypes', y = 'Loadings')
PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
CY <- ggplot(PCLD2, aes(x = reorder(Phenos, value), y = value)) + geom_bar(stat = 'identity', color = 'black') + 
theme(axis.text.y = element_text(angle=45)) + labs(x = 'Phenotypes', y = 'Loadings') + 
coord_flip()
C <- ggord::ggord(fit_lda,df$Strain,veclsz=NA,labcol=NA) + theme_classic(base_size=22) + theme(legend.position = 'top')
#element_text(angle = 45, vjust = 0.5, hjust=1) #element_blank
blankPlot <- ggplot() + geom_blank(aes(1,1)) + theme_void()
CC <- gridExtra::grid.arrange(CY,C,blankPlot,CX, ncol=2,nrow=2,widths = c(1,4), heights=c(4,1))
dev.print(pdf,'../Temp3/lda-ex1.pdf',width=9,height=9)


ggsave('../Temp2/cluster.pdf',width=13,height=13)
ggpubr::ggscatter(lda_result,x='lda.LD1',y='lda.LD2',color='Strain') + 
stat_conf_ellipse(aes(color = Strain, fill = Strain), alpha = 0.1, geom = "polygon", bary=TRUE) + labs(x='LD1', y='LD2') 

kmu.list <- list(); tmp <- numeric(); nclusters <- 4
invisible(lapply(1:choose(20,nclusters/2), function(c) {kmu.list[[c]] <<- kmeans(df1, centers = df1[sample(nrow(df1),nclusters),], 
nstart = 25);tmp[c] <<- as.numeric(kmu.list[[c]]['tot.withinss']);}))
mykmeans <- kmu.list[[which.min(tmp)]]




#####
df <- data_per_strain[,names(data_per_strain) %in% c('Strain',Phenos.lin)]
rownames(df) <- df[,1]
df <- df[,-1]
names(df) <- Phenos.lin.Nomen
df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))

kmu.list <- list(); tmp <- numeric(); nclusters <- 2
invisible(lapply(1:choose(nrow(df),nclusters), function(c) {kmu.list[[c]] <<- kmeans(df, centers = df[sample(nrow(df),nclusters),], 
	nstart = 25);tmp[c] <<- as.numeric(kmu.list[[c]]['tot.withinss']);}))

fviz_cluster(kmu.list[[which.min(tmp)]], data = df, palette = c("#2E9FDF", "#00AFBB","#FC4E07"),
	star.plot = TRUE, ggtheme = theme_minimal(), repel = TRUE, axes = c(2,1), main = 'kmeans')

Mutants <- c("Tmco1<em1J> -/+", "Zbtb43<tm1b(KOMP)Mbp> -/-", "Arpc5l<tm1.1(KOMP)Vlcg> -/+")
df <- data_per_animal[data_per_animal$Strain %in% Mutants,]
df$Strain <- droplevels(df$Strain)