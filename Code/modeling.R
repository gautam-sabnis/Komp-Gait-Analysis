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
#data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% c('speed_20_ang_vel_neg20',
#	'speed_25_ang_vel_neg20'),]
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% 
vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]

data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= median(as.numeric(table(data_per_stride$MouseID)))]),]
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin,'BodyLength')], by = data_per_stride[c("MouseID")], FUN = mean)
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
data_per_animal <- cbind(Strain, TestAge, TestDate, Sex, data_per_animal)

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

#Around 30 C57BL/6NJs are removed as mvoutliers. 

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

#This function uses a linear mixed model approach at the stride level
komp_lmer <- function(CtrlStrain="C57BL/6NJ", freffects='BodyLength + Genotype + Sex + (1|TestDate)', 
	gBG='all'){	
	Mutants <- if(gBG=='all') setdiff(unique(data_per_stride$Strain),CtrlStrain) else
		setdiff(unique(data_per_stride[data_per_stride$BG==gBG,'Strain']),CtrlStrain)
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		ifelse(length(CtrlIDs) < 10, print("Number of Controls < 10"), {
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength',Phenos.lin)]
			df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    		df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
			df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 2)), ]
			df$Strain <- droplevels(df$Strain)
			df$TestDate <- droplevels(df$TestDate)
			possibleError <- tryCatch({
    			df$Sex <- relevel(factor(df$Sex), ref = "Male")
    			df$TestAge <- as.factor(df$TestAge)
    			df$TestDate <- as.factor(df$TestDate)
    			df$TestDate <- droplevels(df$TestDate)
    			df$MouseID <- droplevels(df$MouseID)
    			df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    			#df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
    			}, error = function(e) e
			)
			if (!inherits(possibleError, "error")){
    			FReffects <- freffects
    			formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    			fits <- sapply(formulas, function(x) lmer(formula=x, data = df, REML = FALSE,control = lmerControl(optimizer ='optimx', 
    				optCtrl=list(method='nlminb'))),simplify=FALSE)
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
	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.p <- Heatmap((lodGen), row_names_gp = gpar(fontsize = 12),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif[i,j]),x, y, gp = gpar(fontsize = 8))})

	col_fun <- circlize::colorRamp2(c(min(esizeGen),0,max(esizeGen)),c("#542788","#f7f7f7","#d73027"))
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

#This function uses a multivariate linear model approach at the animal level. 
komp_mv <- function(CtrlStrain="C57BL/6NJ", freffects='BodyLength + Genotype + Sex', 
	gBG='all'){	
	Mutants <- if(gBG=='all') setdiff(unique(data_per_animal$Strain),CtrlStrain) else
		setdiff(unique(data_per_animal[data_per_animal$BG==gBG,'Strain']),CtrlStrain)
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		ifelse(length(CtrlIDs) < 10, print("Number of Controls < 10"), {
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','BodyLength',Phenos.lin)]
			df$Strain <- droplevels(df$Strain)
			df$Strain <- relevel(df$Strain, ref = CtrlStrain)
			possibleError <- tryCatch({
				df$Sex <- relevel(factor(df$Sex), ref = "Male")
    			df$TestAge <- as.factor(df$TestAge)
    			df$MouseID <- as.factor(df$MouseID)
    			df$MouseID <- droplevels(df$MouseID)
    			df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
    			}, error = function(e) e
			)
			if (!inherits(possibleError, "error")){
    			fits <- lm(cbind(speed,limb_duty_factor,step_length1,step_width,stride_length,
						base_tail_lateral_displacement,tip_tail_lateral_displacement,nose_lateral_displacement) ~ BodyLength + 
    					Strain, data = df)
    			pvalGen[s,] <- sapply(seq_along(Phenos.lin), function(x) summary(fits)[[x]]$coefficients[paste0('Strain',Mutants[s]),'Pr(>|t|)']) 
				esizeGen[s,] <- sapply(seq_along(Phenos.lin), function(x) summary(fits)[[x]]$coefficients[paste0('Strain',Mutants[s]),'Estimate'])
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
	pvalGenadj <- pvalGen
	#pvalGenadj <- sapply(Phenos.lin.Nomen, function(r) p.adjust(pvalGen[,r],method = 'fdr'))
	pvalGenSignif <- apply(pvalGenadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodGen <- apply(pvalGenadj,2, function(x)-log(x))
    
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.lin.Nomen
	col_fun <- circlize::colorRamp2(c(2.98,2.99,4,5),c("#FFFFFF","#FFEDA9",'#FFE42B','#FF8303'))
	ht.Gen.p <- Heatmap((lodGen), row_names_gp = gpar(fontsize = 12),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif[i,j]),x, y, gp = gpar(fontsize = 8))})

	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.e <- Heatmap((esizeGen), row_names_gp = gpar(fontsize = 8),row_names_side = "left", column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    	heatmap_legend_param = list(at = c(round(min(esizeGen)) - 1,0,round(max(esizeGen))+0.5),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen[i, j]) * 0.5*(min(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e
	return(ht.Gen)

}

controlids.df <- komp_select_controls(CtrlStrain="C57BL/6NJ", tw=0) #tw: time window

komp_cluster_byline <- function(CtrlStrain="C57BL/6NJ", gBG='all', method='lda'){	
	Mutants <- if(gBG=='all') setdiff(unique(data_per_animal$Strain),CtrlStrain) else
		setdiff(unique(data_per_animal[data_per_animal$BG==gBG,'Strain']),CtrlStrain)
	if (method=='lda'){
		for (s in 1:length(Mutants)){
			cat('Analyzing Strain', paste0(Mutants[s]), "\n")
			CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
			cat('Number of controls:', length(CtrlIDs), "\n")
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','BodyLength',Phenos.lin)]
			df$Strain <- droplevels(df$Strain)
			df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
			df_svd <- svd(df[,sapply(df,is.numeric)])
			df_pca <- df_svd$u %*% diag(df_svd$d)
			df_lda <- data.frame(Strain = df$Strain, df_pca)
			fit_lda <- lda(Strain ~ ., data = df_lda)
			lda_result <- data.frame(Strain = df$Strain, lda = predict(fit_lda)$x)
			if (ncol(lda_result) > 2){
				print(ggplot(lda_result, aes(lda.LD1,lda.LD2, color=Strain)) + geom_point(alpha=0.5,size=3) + labs(x='LD1',y='LD2') + 
				theme(legend.position='none') + ggtitle(paste0(Mutants[s])))
			} else {
				print(ggplot(lda_result, aes(y = LD1, x = seq(nrow(lda_result)), color = Strain)) + geom_point(size = 3) + 
				labs(x = 'Index') + ggtitle(paste0(Mutants[s])))
			}			
		}} else if(method=='gmm'){
			for (s in 1:length(Mutants)){
				cat('Analyzing Strain', paste0(Mutants[s]), "\n")
				CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
				cat('Number of controls:', length(CtrlIDs), "\n")
				df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','BodyLength',Phenos.lin)]
				df$Strain <- droplevels(df$Strain)
				df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
				df_svd <- svd(df[,sapply(df,is.numeric)])
				df_pca <- df_svd$u %*% diag(df_svd$d)
				mod <- Mclust(df_pca, G = 2)
				mod$classification <- df$Strain
				print(plot(mod, what = 'classification', dimens = c(1,2), xlab = 'PC1', ylab = 'PC2'))
		}} else if(method=='kmeans'){
			for (s in 1:length(Mutants)){
				cat('Analyzing Strain', paste0(Mutants[s]), "\n")
				CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
				cat('Number of controls:', length(CtrlIDs), "\n")
				df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','BodyLength',Phenos.lin)]
				df$Strain <- droplevels(df$Strain)
				df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)], 2, function(x) (x - mean(x))/sd(x))
				df$Genotype <- ifelse(df$Strain=='C57BL/6NJ','Control','Mutant')
				df1 <- df[,Phenos.lin]
				rownames(df1) <- paste0(df$MouseID,"(",df$Genotype,")") 
				kmu.list <- list(); tmp <- numeric(); nclusters <- 3
				invisible(lapply(1:choose(62,nclusters), function(c) {kmu.list[[c]] <<- kmeans(df1, centers = df1[sample(nrow(df1),nclusters),], 
    				nstart = 25);tmp[c] <<- as.numeric(kmu.list[[c]]['tot.withinss']);}))
				mykmeans <- kmu.list[[which.min(tmp)]]
				print(factoextra::fviz_cluster(mykmeans, data = df1, palette = c("#999999", "#E69F00","#000000"),
    			star.plot = TRUE, ggtheme = theme_minimal(), repel = TRUE, axes = c(1,2),
    			geom = 'text',ellipse = TRUE, ellipse.type = 'convex', ellipse.level = 0.80,
    			main = paste0(Mutants[s])) + theme_minimal(base_size = 18))
				}}
}



