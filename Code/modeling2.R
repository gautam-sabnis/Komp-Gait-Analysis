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
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
	'speed_25_ang_vel_neg20'),]
data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) > 100]),]
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = var)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
BodyLength <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'BodyLength'][1])
Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
data_per_animal <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, data_per_animal)

#Filter Strains for which at least 8 animals were tested 
Strains8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 8]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strains8,] 
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains8,]
data_per_animal$Strain <- droplevels(data_per_animal$Strain)
data_per_stride$Strain <- droplevels(data_per_stride$Strain)

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
						ifelse(grepl("Hmgu", data_per_animal$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))

#Remove multivariate outliers from Control Strain
df <- data_per_animal[data_per_animal$Strain %in% "C57BL/6NJ",]
tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])
ids <- df[tmp$outliers,'MouseID']
data_per_animal <<- data_per_animal[-which((data_per_animal$Strain == "C57BL/6NJ") & (data_per_animal$MouseID %in% ids)),]
data_per_stride <<- data_per_stride[-which((data_per_stride$Strain == "C57BL/6NJ") & (data_per_stride$MouseID %in% ids)),]


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

komp_lmer <- function(CtrlStrain="C57BL/6NJ"){	
	#Mutants <- if(gBG=='all') setdiff(unique(data_per_stride$Strain),CtrlStrain) else
	#	setdiff(unique(data_per_stride[data_per_stride$BG==gBG,'Strain']),CtrlStrain)
	Mutants <- setdiff(unique(data_per_animal$Strain),"C57BL/6NJ")
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength',Phenos.lin)]
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
		#df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    	Phenos.lin2 <- setdiff(Phenos.lin,"speed")
    	if (length(levels(df$TestDate)) >= 2){
    		FReffects <- 'speed + BodyLength + Genotype + Sex + (1|TestDate)'
    		formulas <- unname(sapply(Phenos.lin2 ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    		fits <- sapply(formulas, function(x) lmer(formula=x, data = df, REML = FALSE,
    			control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))),
    		simplify=FALSE)
    		pvalGen[s,] <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Genotype','Pr(>F)']) 
			esizeGen[s,] <- sapply(fits, function(mod) S(mod)$fixed.effects['GenotypeMutant','Estimate']) 
    	} else {
    		FReffects <- 'speed + BodyLength + Genotype + Sex'
    		formulas <- unname(sapply(Phenos.lin2 ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    		fits <- lapply(formulas, lm, data = df)
    		pvalGen[s,] <- sapply(seq_along(fits), function(x) Anova(fits[[x]])['Genotype','Pr(>F)']) 
			esizeGen[s,] <- sapply(seq_along(fits), function(x) unname(fits[[x]]$coefficients['GenotypeMutant'])) 
    	}
	}
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

	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.p <- Heatmap((lodGen2), row_names_gp = gpar(fontsize = 10, fontface='italic'),row_names_side = "left", column_names_gp = gpar(fontsize = 12),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title='Mean Phenotypes')

	col_fun <- circlize::colorRamp2(c(min(esizeGen2),0,max(esizeGen2)),c("#542788","#f7f7f7","#d73027"))
	ht.Gen.e <- Heatmap((esizeGen2), row_names_gp = gpar(fontsize = 8, fontface="italic"),row_names_side = "left", column_names_gp = gpar(fontsize = 12),
    	heatmap_legend_param = list(at = c(round(min(esizeGen2)) - 1,0,round(max(esizeGen2))+0.5),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen2[i, j]) * 2 * (min(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e
	return(list(ht.Gen,lodGen2, esizeGen2))

}

tmp <- komp_lmer()

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
	p <- ggplot(df, aes_string(y=phenoname, x='Genotype')) + geom_boxplot(outlier.shape=NA,width=0.7) + 
	geom_jitter(aes(color=Genotype), width=0.1, size = 5, shape = 1, stroke=1.5) + scale_color_manual(values=c("#6a51a3", "#d94801")) + 
	theme_bw(base_size = 26) + theme(legend.position='none') + labs(y = paste0(phenotype)) + facet_grid(~Strain)
	return(p)
}

lapply(seq(length(Phenos.lin)), function(x) {vignette_plot(Phenos.lin.Nomen[x],Phenos.lin[x]); 
ggsave(paste0('../Temp2/vignette-mean-', Phenos.lin.Nomen[x],'.pdf'), width=16, height=4)})


komp_mv <- function(CtrlStrain="C57BL/6NJ", gBG='all'){	
	Mutants <- if(gBG=='all') setdiff(unique(data_per_stride$Strain),CtrlStrain) else
		setdiff(unique(data_per_stride[data_per_stride$BG==gBG,'Strain']),CtrlStrain)
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength',Phenos.lin)]
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
		df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    	fits <- lm(cbind(speed,limb_duty_factor,step_length1,step_width,stride_length,
				temporal_symmetry,base_tail_lateral_displacement,tip_tail_lateral_displacement,
				nose_lateral_displacement) ~ BodyLength + Genotype + Sex, data = df)
    	pvalGen[s,] <- sapply(seq_along(Phenos.lin), function(x) summary(fits)[[x]]$coefficients['GenotypeMutant','Pr(>|t|)']) 
		esizeGen[s,] <- sapply(seq_along(Phenos.lin), function(x) summary(fits)[[x]]$coefficients['GenotypeMutant','Estimate'])
	}

	esizeGen <- as.matrix(esizeGen)
	rownames(pvalGen) <- Mutants
	colnames(pvalGen) <- Phenos.lin.Nomen
	rownames(esizeGen) <- Mutants
	colnames(esizeGen) <- Phenos.lin.Nomen
	pvalGenadj <- pvalGen
	#pvalGenadj <- sapply(Phenos.lin.Nomen, function(r) p.adjust(pvalGen[,r],method = 'fdr'))
	pvalGenSignif <- apply(pvalGenadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodGen <- apply(pvalGenadj,2, function(x)-log(x))
    
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.lin.Nomen
	lodGen <- lodGen[order(row.names(lodGen)),]
	esizeGen <- esizeGen[order(row.names(esizeGen)),]
	pvalGenSignif <- pvalGenSignif[order(row.names(pvalGenSignif)),]
	lodGen2 <- lodGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	esizeGen2 <- esizeGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	pvalGenSignif2 <- pvalGenSignif[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]

	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.p <- Heatmap((lodGen2), row_names_gp = gpar(fontsize = 12),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif2[i,j]),x, y, gp = gpar(fontsize = 8))})

	col_fun <- circlize::colorRamp2(c(min(esizeGen2),0,max(esizeGen2)),c("#542788","#f7f7f7","#d73027"))
	ht.Gen.e <- Heatmap((esizeGen2), row_names_gp = gpar(fontsize = 8),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(round(min(esizeGen2)) - 1,0,round(max(esizeGen2))+0.5),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = 2*(esizeGen2[i, j]) * (min(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen2[i,j]))))}, rect_gp = gpar(type = "none"))
	ht.Gen <- ht.Gen.p + ht.Gen.e
	return(ht.Gen)
}
