df <- data.frame(Strain = data_per_animal$Strain, MouseID = data_per_animal$MouseID, Sex = data_per_animal$Sex)
df <- cbind(df, apply(data_per_animal[,names(data_per_animal) %in% Phenos.lin],2,function(x) (x - mean(x))/sd(x)))
tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin])
df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, 
	Label = paste0(df[,'MouseID']," (", df[,'Strain'], ")"), MouseID = df[,'MouseID'], Strain = df[,'Strain'],
	Outlier = as.numeric(!tmp$wfinal01))
df.out[df.out$MouseID == 'J80962', 'Outlier'] <- -1
df.out[df.out$MouseID == 'J76119', 'Outlier'] <- -1
df.out$Outlier <- as.factor(df.out$Outlier)
out.df <- df.out[df.out$Outlier %in% c(1,-1), names(df.out) %in% c('Strain','MouseID')]
out.df <- out.df[out.df$Strain %in% names(table(data_per_animal$Strain)[table(data_per_animal$Strain) < 5]),]
out.df$Strain <- droplevels(out.df$Strain)
out.df$MouseID <- droplevels(out.df$MouseID)
total_animals <- sapply(seq(length(unique(out.df$Strain))), function(x) table(data_per_animal$Strain)[paste0(unique(out.df$Strain)[x])])
out_animals <- sapply(seq(length(unique(out.df$Strain))), function(x) table(out.df$Strain)[paste0(unique(out.df$Strain)[x])])
df <- data.frame(Proportion = out_animals/total_animals,
	Animals = total_animals)
df <- cbind(Strain = rownames(df), df)

Mutants <- c('Zbtb43-/-')
m <- 1
CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
dfa <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs, ]
dfa$Genotype <- ifelse(dfa$Strain == 'C57BL/6NJ', 'Control', 'Mutant')
dfa$MouseID <- droplevels(dfa$MouseID)
dfa <- dfa[dfa$TestDate %in% names(which(table(dfa$TestDate, dfa$Genotype)[,2] >= 1)), ]
dfa$TestDate <- droplevels(dfa$TestDate)
#CtrlIDs <- sample(unique(dfa[dfa$Genotype=='Control','MouseID']), length(unique(dfa[dfa$Genotype=='Mutant','MouseID'])))
CtrlIDs <- setdiff(unique(dfa[dfa$Genotype == 'Control', 'MouseID']), unique(out.df$MouseID))
df <- data_per_stride[data_per_stride$Strain == Mutants[m],c('MouseID','BodyLength','Sex',Phenos.lin)]
df <- rbind(df,data_per_stride[data_per_stride$MouseID %in% CtrlIDs, c('MouseID', 'BodyLength','Sex',Phenos.lin)])
df$MouseID <- droplevels(df$MouseID)
df$Strain <- ifelse(df$MouseID %in% CtrlIDs,'Control','Mutant')
df$Outlier <- as.factor(ifelse(df$MouseID %in% out.df$MouseID, 1,0)) 
#df$Outlier <- as.factor(ifelse(df$MouseID %in% 'J67783', 1,0)) 


formula <- 'BodyLength + Strain + (1|MouseID)'
formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", formula), simplify=TRUE)) 
fits <- sapply(formulas, function(x) lmer(formula=x, data = df, REML = FALSE,
    	control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))),
simplify=FALSE)
sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Strain','Pr(>F)']) 
sapply(seq_along(fits), function(x) summary(fits[[x]])$coefficients[3,'Estimate']) 


df.tmp <- df[df$Strain %in% 'Control',]
df.tmp <- rbind(df.tmp, df[df$MouseID %in% 'J60202',])
mod <- lmer(step_width ~ BodyLength + (1|MouseID), data = df.tmp)




lapply(seq(unique(out.df$Strain)), function(x) {

	if (df$Strain == unique(out.df$Strain)[x] & 





})









#########
df.tmp <- data_per_animal[,which(names(data_per_animal) %in% Phenos.lin)]
df.tmp <- scale(df.tmp, center=TRUE, scale=TRUE)
colnames(df.tmp) <- Phenos.lin.Nomen
tmp <- cor(df.tmp)
corrplot::corrplot(tmp, method = 'number',type='lower', tl.col='black',tl.cex=1.2, cl.sex=1.2)




#######
df <- data_per_strain
colnames(df) <- c('Strain',Phenos.lin.Nomen)
invisible(lapply(seq(length(Phenos.lin.Nomen)), function(x) assign(paste0("p",x),
	ggplot(df,aes_string(y=paste0(Phenos.lin.Nomen[x]))) + 
	geom_boxplot() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) , inherits=TRUE)))

p1|p2|p3|p4|p5|p6|p7|p8|p9

#########
Mutants <- c('Arpc5l-/+','Fam120b')

#########
Strains1 <- as.character(Strains1)
Strains2 <- as.character(Strains2)
Strains3 <- as.character(Strains3)
Strains1v <- as.character(Strains1v)
Strains2v <- as.character(Strains2v)
Strains3v <- as.character(Strains3v)
MvoutStrains <- as.character(MvoutStrains)

AllStrains <- c(Strains1,Strains2,Strains3,Strains1v,Strains2v,Strains3v,MvoutStrains)
AllStrains <- sort(unique(AllStrains))
dfStrains <- data.frame(Strains = AllStrains)

dfStrains$M1 <- as.factor(ifelse(dfStrains$Strains %in% Strains1, 1, 0)) 
dfStrains$M2 <- as.factor(ifelse(dfStrains$Strains %in% Strains2, 1, 0))
dfStrains$M3 <- as.factor(ifelse(dfStrains$Strains %in% Strains3, 1, 0))
dfStrains$V1 <- as.factor(ifelse(dfStrains$Strains %in% Strains1v, 1, 0))
dfStrains$V2 <- as.factor(ifelse(dfStrains$Strains %in% Strains2v, 1, 0))
dfStrains$V3 <- as.factor(ifelse(dfStrains$Strains %in% Strains3v, 1, 0))
dfStrains$Mvout <- as.factor(ifelse(dfStrains$Strains %in% MvoutStrains, 1, 0))