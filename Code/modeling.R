Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry","base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
Phenos.lin.Nomen <- c("Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
"Temporal Symmetry","Base Tail LD","Tip Tail LD","Nose LD")
Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
"nose_lateral_displacement_phase")

data_per_stride <- read.delim(snakemake@input[[1]], stringsAsFactors = FALSE)
names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
names(data_per_stride)[names(data_per_stride) == 'OFA_Strain.Name'] <- 'Strain'
names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))

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

#Filter Strains for which at least 8 animals were tested 
Strains8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 8]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strains8,] 
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains8,]

#Remove outliers  
invisible(sapply(seq(length(unique(data_per_animal$Strain))), function(s) Map(function(p) {
	vals <- data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],c('MouseID',Phenos.lin[p])];
	outliers <- boxplot.stats(vals[,2])$out
	ids <- match(outliers, vals[,2]) 
	data_per_animal[data_per_animal$Strain == unique(data_per_animal$Strain)[s],][paste0(Phenos.lin[p])][ids,] <<- NA
	data_per_stride[(data_per_stride$Strain == unique(data_per_animal$Strain)[s]) & 
	(data_per_stride$MouseID %in% vals[ids,1]), Phenos.lin[p]] <<- NA
}, seq(length(Phenos.lin)))))

#<em1J> - em1J, (KOMP) - Mbp, Wtsi, Vlcg, (EUCOMM) - Hmgu 
data_per_stride$BG <- as.factor(ifelse(grepl("Mbp", data_per_stride$Strain, fixed = TRUE), "Mbp",
						ifelse(grepl("Wtsi", data_per_stride$Strain, fixed = TRUE), "Wtsi",
						ifelse(grepl("Vlcg", data_per_stride$Strain, fixed = TRUE),"Vlcg",
						ifelse(grepl("em1J", data_per_stride$Strain, fixed = TRUE),"em1J",
						ifelse(grepl("Hmgu", data_per_stride$Strain, fixed = TRUE), "Hmgu","C57BL/6NJ"))))))

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

controlids.df <- komp_select_controls(CtrlStrain="C57BL/6NJ-Ccdc30<em1J>/J",tw=31)

komp_lmer <- function(CtrlStrain="C57BL/6NJ", freffects='BodyLength + Genotype + TestAge + (1|MouseID)', 
	gBG='all'){	
	Mutants <- if(gBG=='all') setdiff(unique(data_per_stride$Strain),CtrlStrain) else
		setdiff(unique(data_per_stride[data_per_stride$BG==gBG,'Strain']),CtrlStrain)
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		ifelse(length(CtrlIDs) < 10, print("Number of Controls < 10"), {
			df <- data_per_stride[data_per_stride$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','BodyLength',Phenos.lin)]
			df$Strain <- droplevels(df$Strain)
			possibleError <- tryCatch({
				df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    			df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
    			df$Sex <- relevel(factor(df$Sex), ref = "Male")
    			df$TestAge <- as.factor(df$TestAge)
    			df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
    			FReffects <- freffects
    			ormulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    			fits <- sapply(formulas, function(x) lmer(formula=x, data = df, REML = FALSE,control = lmerControl(optimizer ='optimx', 
    				optCtrl=list(method='nlminb'))),simplify=FALSE)
			}, error = function(e) e
			)
			if (!inherits(possibleError, "error")){
    			pvalGen[s,] <- sapply(seq_along(fits), function(x) anova(unname(fits[x])[[1]],type='II')['Genotype','Pr(>F)']) 
				esizeGen[s,] <- sapply(fits, function(mod) S(mod)$fixed.effects['GenotypeMutant','Estimate']) 
			} else{
				pvalGen[s,] <- 0
				esizeGen[s,] <- 0 
			}
    		}
		)
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
    
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.lin.Nomen
	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#FFFFFF","#FFEDA9",'#FFE42B','#FF8303'))
	ht.Gen.p <- Heatmap((lodGen), row_names_gp = gpar(fontsize = 12),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif[i,j]),x, y, gp = gpar(fontsize = 8))})

	col_fun <- colorRamp2(c(min(esizeGen),0,max(esizeGen)),c("#542788","#f7f7f7","#d73027"))
	ht.Gen.e <- Heatmap((esizeGen), row_names_gp = gpar(fontsize = 8),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(round(min(esizeGen)) - 1,0,round(max(esizeGen))+0.5),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen[i, j]) * (min(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e
	return(ht.Gen)

}
