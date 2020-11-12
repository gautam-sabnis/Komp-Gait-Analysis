Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry", "base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
	"TS","Base Tail LD","Tip Tail LD","Nose LD")
Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
"nose_lateral_displacement_phase")

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
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
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

komp_lmer <- function(CtrlStrain="C57BL/6NJ", model){	
	if (model == 'M2' || model == 'M3'){
		Phenos.lin <- setdiff(Phenos.lin,"speed");
		Phenos.lin.Nomen <- setdiff(Phenos.lin.Nomen,"Speed");
	} else {
		Phenos.lin <- Phenos.lin;
		Phenos.lin.Nomen <- Phenos.lin.Nomen;
	}

	Mutants <- setdiff(unique(data_per_animal$Strain),"C57BL/6NJ")
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	pvalSex <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeSex <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	#pvalint <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	#esizeint <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		if (model == 'M2' || model == 'M3'){
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength','speed',Phenos.lin)]	
		} else {
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength',Phenos.lin)]
		
		}
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
		#if (sum(table(df$TestDate,df$Genotype)[,1] & table(df$TestDate,df$Genotype)[,2] >= 1) >= 1){
    	#	df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
    	#} else {
    	#	df <- df
    	#}
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    	if (model=='M3')
    		{FReffects.lmer <- 'BodyLength + speed + Sex + Genotype + (1|TestDate)';
    		 FReffects <- 'BodyLength + speed + Sex + Genotype'}
    	else if (model == 'M2')
    		{FReffects.lmer <- 'speed + Sex + Genotype + (1|TestDate)';
    		 FReffects <- 'speed + Sex + Genotype'}
    	else {FReffects.lmer <- 'BodyLength + Sex + Genotype + (1|TestDate)';
    		 FReffects <- 'BodyLength + Sex + Genotype'}
    	#if (Mutants[s] == "Bzw1-/-") {FReffects <- 'BodyLength + Genotype'}
    	if (length(levels(df$TestDate)) >= 2){
    		formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects.lmer), simplify=TRUE)) 
    		fits <- sapply(formulas, function(x) lmer(formula=x, data = df, REML = FALSE,
    			control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))),
    		simplify=FALSE)
    		pvalGen[s,] <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Genotype','Pr(>F)']) 
			esizeGen[s,] <- sapply(fits, function(mod) S(mod)$fixed.effects['GenotypeMutant','Estimate']) 
			if (min(table(df$Genotype,df$Sex)['Mutant',]) > 1){
				pvalSex[s,] <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Sex','Pr(>F)']) 
				esizeSex[s,] <- sapply(seq_along(fits), function(x) S(fits[[x]])$fixed.effects['SexFemale','Estimate']) 
				#pvalint[s,] <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Sex:Genotype','Pr(>F)']) 
				#esizeint[s,] <- sapply(seq_along(fits), function(x) S(fits[[x]])$fixed.effects['SexFemale:GenotypeMutant','Estimate']) 
			} else {
				pvalSex[s,] <- rep(1,length(Phenos.lin)) 
				esizeSex[s,] <- rep(0,length(Phenos.lin)) 
				#pvalint[s,] <- rep(0,length(Phenos.lin))
				#esizeint[s,] <- rep(0,length(Phenos.lin))
			}
			
    	} else {
    		formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    		fits <- lapply(formulas, lm, data = df)
    		pvalGen[s,] <- sapply(seq_along(fits), function(x) Anova(fits[[x]])['Genotype','Pr(>F)']) 
			esizeGen[s,] <- sapply(seq_along(fits), function(x) unname(fits[[x]]$coefficients['GenotypeMutant'])) 
			pvalSex[s,] <- sapply(seq_along(fits), function(x) Anova(fits[[x]])['Sex','Pr(>F)']) 
			esizeSex[s,] <- sapply(seq_along(fits), function(x) unname(fits[[x]]$coefficients['SexFemale'])) 

    	}
	}

	pvalGen <- pvalGen[complete.cases(pvalGen),]
	esizeGen <- esizeGen[complete.cases(esizeGen),]
	pvalSex <- pvalSex[complete.cases(pvalSex),]
	esizeSex <- esizeSex[complete.cases(esizeSex),]
	#Genotype
	Mutants <- setdiff(unique(data_per_animal$Strain),"C57BL/6NJ")
	esizeGen <- as.matrix(esizeGen)
	rownames(pvalGen) <- Mutants
	colnames(pvalGen) <- Phenos.lin.Nomen
	rownames(esizeGen) <- Mutants
	colnames(esizeGen) <- Phenos.lin.Nomen
	pvalGenadj <- sapply(Phenos.lin.Nomen, function(r) p.adjust(pvalGen[,r],method = 'fdr'))
	pvalGenSignif <- apply(pvalGenadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodGen <- apply(pvalGenadj,2, function(x)-log(x))
    rownames(pvalGenSignif) <- Mutants
    
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.lin.Nomen
	lodGen <- lodGen[order(row.names(lodGen)),]
	esizeGen <- esizeGen[order(row.names(esizeGen)),]
	pvalGenSignif <- pvalGenSignif[order(row.names(pvalGenSignif)),]
	lodGen2 <- lodGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	esizeGen2 <- esizeGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	pvalGenSignif2 <- pvalGenSignif[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	Mutants <- rownames(lodGen2)
	row <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][1]
	col <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][2]
	esizeGen2[row,col] <- tail(sort(abs(esizeGen2)),1)[1]

	summ.df.Gen <- data.frame(Phenos.lin.Nomen, 
	pink = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99 & abs(esizeGen2[,x]) >= 0.2))), 
	black = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99))))

	pinkStrainsGen <- names(which(apply((lodGen2 > 2.99 | lodGen2 > 2.99) & (abs(esizeGen2) >= 0.2 | abs(esizeGen2) >= 0.2),1,any)))
	blackStrainsGen <- names(apply(lodGen2 > 2.99,1,any))

	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.p <- Heatmap((lodGen2), row_names_gp = gpar(fontsize = 16, fontface='italic',
		col=ifelse(apply(((lodGen2 > 2.99 | lodGen2 > 2.99) & (abs(esizeGen2) >= 0.2 | abs(esizeGen2 >= 0.2))),1, function(x) any(x == TRUE)), '#dd3497', 'black')),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title=paste0(model,': Variance Phenotypes'))

	col_fun <- circlize::colorRamp2(c(min(esizeGen2),0,max(esizeGen2)),c("#542788","#f7f7f7","#d73027"))
	ht.Gen.e <- Heatmap((esizeGen2), row_names_gp = gpar(fontsize = 16, fontface="italic"),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(round(min(esizeGen2),2),0,round(max(esizeGen2),2)),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen2[i, j]) * 0.8 * (max(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e

	#Sex
	#Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Bzw1-/-"))
	Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ"))
	esizeSex <- as.matrix(esizeSex)
	rownames(pvalSex) <- Mutants
	colnames(pvalSex) <- Phenos.lin.Nomen
	rownames(esizeSex) <- Mutants
	colnames(esizeSex) <- Phenos.lin.Nomen
	pvalSexadj <- sapply(Phenos.lin.Nomen, function(r) p.adjust(pvalSex[,r],method = 'fdr'))
	pvalSexSignif <- apply(pvalSexadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodSex <- apply(pvalSexadj,2, function(x)-log(x))
    rownames(pvalSexSignif) <- Mutants
    
    rownames(lodSex) <- Mutants
	colnames(lodSex) <- Phenos.lin.Nomen
	lodSex <- lodSex[order(row.names(lodSex)),]
	esizeSex <- esizeSex[order(row.names(esizeSex)),]
	pvalSexSignif <- pvalSexSignif[order(row.names(pvalSexSignif)),]
	lodSex2 <- lodSex[sapply(seq(nrow(lodSex)), function(x) any(lodSex[x,] >= 2.98)),]
	esizeSex2 <- esizeSex[sapply(seq(nrow(lodSex)), function(x) any(lodSex[x,] >= 2.98)),]
	pvalSexSignif2 <- pvalSexSignif[sapply(seq(nrow(lodSex)), function(x) any(lodSex[x,] >= 2.98)),]
	Mutants <- rownames(lodSex2)
	row <- which(abs(esizeSex2)==max(abs(esizeSex2)), arr.ind=TRUE)[1:2][1]
	col <- which(abs(esizeSex2)==max(abs(esizeSex2)), arr.ind=TRUE)[1:2][2]
	esizeSex2[row,col] <- tail(sort(abs(esizeSex2)),1)[1]

    summ.df.Sex <- data.frame(Phenos.lin.Nomen, 
	pink = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodSex2[,x] > 2.99 & abs(esizeSex2[,x]) >= 0.1))), 
	black = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodSex2[,x] > 2.99))))

	pinkStrainsSex <- names(which(apply((lodSex2 > 2.99 | lodSex2 > 2.99) & (abs(esizeSex2) >= 0.1 | abs(esizeSex2) >= 0.1),1,any)))
	blackStrainsSex <- names(apply(lodGen2 > 2.99,1,any))
	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Sex.p <- Heatmap((lodSex2), row_names_gp = gpar(fontsize = 16, fontface='italic',
		col=ifelse(apply(((lodSex2 > 2.99 | lodSex2 > 2.99) & (abs(esizeSex2) >= 0.2 | abs(esizeSex2 >= 0.2))),1, function(x) any(x == TRUE)), '#dd3497', 'black')),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalSexSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title=paste0(model,': Variance Phenotypes'))

	col_fun <- circlize::colorRamp2(c(min(esizeSex2),0,max(esizeSex2)),c("#542788","#f7f7f7","#d73027"))
	ht.Sex.e <- Heatmap((esizeSex2), row_names_gp = gpar(fontsize = 16, fontface="italic"),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(round(min(esizeSex2),2),0,round(max(esizeSex2),2)),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeSex2[i, j]) * (max(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeSex2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Sex <- ht.Sex.p + ht.Sex.e


	return(list(ht.Gen,summ.df.Gen,pinkStrainsGen,blackStrainsGen,ht.Sex,summ.df.Sex,pinkStrainsSex,
		blackStrainsSex))

}

tmp1v <- komp_lmer(model = 'M1')
tmp2v <- komp_lmer(model = 'M2')
tmp3v <- komp_lmer(model = 'M3')


Strains1 <- tmp1v[[4]]
Strains2 <- tmp2v[[4]]
Strains3 <- tmp3v[[4]]

Strains1 <- gsub("*./.*","",Strains1)
Strains2 <- gsub("*./.*","",Strains2)
Strains3 <- gsub("*./.*","",Strains3)



tmp7 <- komp_lmer(model='M3')
tmp14 <- komp_lmer(model='M3')
tmp21 <- komp_lmer(model='M3')
tmp28 <- komp_lmer(model='M3')
#Volcano plot 
dfp.melt <- reshape::melt(tmp[[2]], id.vars = 'Strain')
names(dfp.melt) <- c('Strain','Phenotype','LOD')
dfe.melt <- reshape::melt(tmp[[3]], id.vars = 'Strain')
names(dfe.melt) <- c('Strain', 'Phenotype', 'Esize')
df <- dfp.melt
df <- cbind(df,esize=dfe.melt$Esize)
ggplot(df, aes(x = esize, y=LOD, color=Phenotype)) + geom_point(shape = 1,size=2,stroke=1) + 
ggrepel::geom_text_repel(aes(label=ifelse((LOD >= 2.99 & abs(esize) >= 0.2),paste0(Strain), "")),box.padding = 0.5) + 
geom_vline(xintercept = c(-0.2,0.2), linetype='dashed') + geom_hline(yintercept = 2.99, linetype='dashed') + theme_bw(base_size=14) + 
theme(legend.position = 'top') + labs(x = 'Effect Size', y = 'LOD Score')
ggsave('../Temp5/volcano-var-M3.pdf', width=8, height=8)

dev.print(pdf,'../Temp5/mean-lmer-M3-5.pdf', height = 11.75, width = 6.4)
dev.print(pdf,'../Temp5/mean-lmer-M1-5.pdf', height = 11.4, width = 6.5)
dev.print(pdf,'../Temp5/mean-lmer-M2-5.pdf', height = 11.43, width = 5.72)
dev.print(pdf,'../Temp4/variance-lmer-M2.pdf', height = 13, width = 8.5)
dev.print(pdf,'../Temp4/variance-lmer-M1.pdf', height = 14, width = 8.5)
dev.print(pdf,'../Temp5/variance-lmer-M3-5.pdf', height = 11.75, width = 6.4)


vignette_plot <- function(phenotype, phenoname, esizeGen2 = tmp[[3]], model){
	if (model == 'M2' || model == 'M3'){
		Phenos.lin <- setdiff(Phenos.lin,"speed");
		Phenos.lin.Nomen <- setdiff(Phenos.lin.Nomen,"Speed");
	} else {
		Phenos.lin <- Phenos.lin;
		Phenos.lin.Nomen <- Phenos.lin.Nomen;
	}
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
	p <- ggplot(df, aes_string(y=phenoname, x='Genotype')) + geom_boxplot(outlier.shape=NA,width=0.7) + 
	geom_jitter(aes(color=Genotype), width=0.1, size = 5, shape = 1, stroke=1.5) + scale_color_manual(values=c("#6a51a3", "#d94801")) + 
	theme_bw(base_size = 32) + theme(legend.position='none') + labs(y = paste0(phenotype)) + facet_grid(~Strain)
	return(p)
}

lapply(seq(length(Phenos.lin)), function(x) {vignette_plot(Phenos.lin.Nomen[x],Phenos.lin[x],model='M3'); 
ggsave(paste0('../Temp5/Vignettes/vignette-var-M3-', Phenos.lin.Nomen[x],'.pdf'), width=16, height=4)})

check_assumptions <- function(CtrlStrain = 'C57BL/6NJ', model){
	if (model == 'M2' || model == 'M3'){
		Phenos.lin <- setdiff(Phenos.lin,"speed");
		Phenos.lin.Nomen <- setdiff(Phenos.lin.Nomen,"Speed");
	} else {
		Phenos.lin <- Phenos.lin;
		Phenos.lin.Nomen <- Phenos.lin.Nomen;
	}

	Mutants <- setdiff(unique(data_per_animal$Strain),"C57BL/6NJ")
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	df.diag <- list()
	df.re <- list()
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		if (model == 'M2' || model == 'M3'){
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength','speed',Phenos.lin)]	
		} else {
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength',Phenos.lin)]
		
		}
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
    	#if (sum(table(df$TestDate,df$Genotype)[,1] & table(df$TestDate,df$Genotype)[,2] >= 1) >= 1){
    	#	df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
    	#} else {
    	#	df <- df
    	#}
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)

    	if (model=='M3')
    		{FReffects.lmer <- 'BodyLength + speed + Genotype + Sex + (1|TestDate)';
    		 FReffects <- 'BodyLength + speed + Genotype + Sex'} else if (model == 'M2')
    		{FReffects.lmer <- 'speed + Genotype + Sex + (1|TestDate)';
    		 FReffects <- 'speed + Genotype + Sex'} else {FReffects.lmer <- 'BodyLength + Genotype + Sex + (1|TestDate)';
    		 FReffects <- 'BodyLength + Genotype + Sex'}
    	#cat("Strain: , TestDate:", paste(Mutants[s],length(levels(df$TestDate))), "\n")
    	list.diagnostics <- list()
    	list.re <- list()
    	if (length(levels(df$TestDate)) >= 2){
    		formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects.lmer), simplify=TRUE)) 
    		fits <- sapply(formulas, function(x) lmer(formula=x, data = df, REML = FALSE,
    			control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))),
    		simplify=FALSE)
    		list.diagnostics <- lapply(seq(length(Phenos.lin)), function(x) data.frame(raw_cond = residuals(fits[[x]]), f = fitted(fits[[x]]), 
    		Phenotype = rep(Phenos.lin.Nomen[x], each=length(residuals(fits[[x]]))), Strain = rep(Mutants[s], length(residuals(fits[[x]]))),
    		BodyLength = df[(!is.na(df[Phenos.lin[x]]) & !is.na(df$speed)),'BodyLength'], Speed = df[(!is.na(df[Phenos.lin[x]]) & !is.na(df$speed)),'speed']))
    		list.re <- lapply(seq(length(Phenos.lin)), function(x) data.frame(TestDate = rownames(ranef(fits[[1]])$TestDate), re = ranef(fits[[x]])$TestDate[[1]],
    			Phenotype = rep(Phenos.lin.Nomen[x], each=length(ranef(fits[[x]])$TestDate[[1]])), 
    			Strain = rep(Mutants[s],length(ranef(fits[[x]])$TestDate[[1]]))))
    		df.diag[[s]] <- do.call(rbind, list.diagnostics)
    		df.re[[s]] <- do.call(rbind, list.re)
    		
    	} else {
    		formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    		fits <- lapply(formulas, lm, data = df)
    	}
	}
	return(list(df.diag,df.re)) 
}

tmp <- check_assumptions(model='M3')

#Normality of residuals
#Phenos.lin<- setdiff(Phenos.lin, 'speed') uncomment for M3,M2
resid.pvals <- sapply(tmp[[1]], function(l) sapply(seq(length(Phenos.lin)), function(x) 
	shapiro.test(l[l$Phenotype==Phenos.lin.Nomen[x],'raw_cond'])$p.value))
rownames(resid.pvals) <- Phenos.lin.Nomen
colnames(resid.pvals) <- Mutants
(which(resid.pvals <= 0.05, arr.ind=TRUE))
tmp2 <- t(resid.pvals)
tmp2 <- tmp2[order(rownames(tmp2)),]
df1 <- tmp2[1:44,]
df2 <- tmp2[45:88,]
df3 <- tmp2[89:134,]

#Normality of random effects
ranef.pvals <- sapply(tmp[[2]], function(l) sapply(seq(length(Phenos.lin)), function(x) 
	ifelse(var(l[l$Phenotype==Phenos.lin.Nomen[x],'re']) == 0, NA,
		shapiro.test(l[l$Phenotype==Phenos.lin.Nomen[x],'re'])$p.value)))
rownames(ranef.pvals) <- Phenos.lin.Nomen
colnames(ranef.pvals) <- Mutants
which(ranef.pvals <= 0.05, arr.ind=TRUE)

tmp2 <- do.call(rbind,tmp[[1]][1:length(Mutants)])

#Checking linearity assumption of covariate Speed in the model
sapply(seq(length(Phenos.lin)), function(x) {p <- ggplot(tmp2[tmp2$Phenotype %in% c(paste0(Phenos.lin.Nomen[x])),], aes(x = Speed, y = raw_cond)) + geom_point() + 
geom_smooth(method='loess',color='red') + facet_wrap(~Strain , scales='free') + 
ggtitle(paste0(Phenos.lin.Nomen[x])) + labs(x = 'Speed', y = 'Conditional Raw Residuals');
plot(p)})

#Checking constant variance assumption in the model
sapply(seq(length(Phenos.lin)), function(x) {p <- ggplot(tmp2[tmp2$Phenotype %in% c(paste0(Phenos.lin.Nomen[x])),], aes(x = f, y = raw_cond)) + geom_point() + 
geom_smooth(method='loess',color='red') + facet_wrap(~Strain , scales='free') + 
ggtitle(paste0(Phenos.lin.Nomen[x])) + labs(x = 'Fitted values', y = 'Conditional Raw Residuals');
plot(p)})


#Post-process plots
tmp.df <- data.frame(Phenos.lin.Nomen, 
	pink = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99 & abs(esizeGen2[,x]) >= 0.2))), 
	black = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99))))