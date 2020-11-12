Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
"nose_lateral_displacement_phase")
Phenos.circ.Nomen <- c('Base Tail Phase', 'Tip Tail Phase', 'Nose Phase')

#setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#data_per_stride <- read.delim('../Data/kompdf-corr', stringsAsFactors = FALSE)


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
data_per_stride[names(data_per_stride) %in% Phenos.circ] <- lapply(data_per_stride[names(data_per_stride) %in% Phenos.circ],
    function(x) x*(2*pi))
data_per_stride[names(data_per_stride) %in% Phenos.circ] <- lapply(data_per_stride[names(data_per_stride) %in% Phenos.circ],
    function(x) circular(x, type = 'angles', units = 'radians'))

#Focus on certain speed bins 
data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% 
vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= 
		summary(as.numeric(table(data_per_stride$MouseID)))[2]]),]
data_per_stride$Strain <- vapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]), FUN.VALUE="1")
data_per_stride$Strain <- vapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]), FUN.VALUE=character(1))
data_per_stride$Strain <- as.factor(data_per_stride$Strain)
data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.circ)], by = data_per_stride[c("MouseID")], FUN = mean)
BodyLength <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'BodyLength'][1])
speed <- sapply(seq(dim(data_per_animal)[1]), function(x) mean(data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'speed']))
Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
data_per_animal <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, speed, data_per_animal)

#Filter Strains for which at least 8 animals were tested 
Strains8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 5]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strains8,] 
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains8,]
data_per_animal$Strain <- droplevels(data_per_animal$Strain)
data_per_stride$Strain <- droplevels(data_per_stride$Strain)

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

komp_circular <- function(CtrlStrain="C57BL/6NJ", model){	
	if (model == 'M1'){
		#Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Enpp5-/-","Elavl1-/+",
		#"Bcl11b-/+","Rmnd5b-/-","Top3a-/+","Mrps22-/+","Zfp579-/-","Slc50a1-/-","Pdgfrl-/-","Cdc14a-/+"))
		Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ", "Spin1-/+","Enpp5-/-","Rxfp3-/-",
			"Mrps22-/+","Zfp579-/-","Ofcc1-/-","Slc50a1-/- ","Slc50a1-/-","Rrs1-/+","Tmco1-/+","Xpnpep3-/+",
			"Rai14-/+")) #Variance
	} else if (model == 'M2'){
		#Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Rad51ap1-/-","Gpr15-/-",
		#"Mettl10-/-","Rxfp3-/-","Pde4c-/-","Zfp422-/-","Hpse2-/+","Mrps22-/+","Cdc14a-/+",
		#"Xpnpep3-/+","Atp5o-/+"))
		Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Spin1-/+","Zfp422-/-","Mrps22-/+",
			"Tomm40-/+","Rai14-/+")) #Variance
	} else {
		Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Tmem222-/+","Alg11-/+","Rad51ap1-/-",
			"Kifc2-/-","Ing2-/+","Myom1-/-","Cldn13-/-","Steap3-/-","Coq8b-/+","Fbxo30-/+","Tada1-/+",
			"Rab3gap2-/-","Ppip5k2-/+","Abi3-/-","Xpnpep3-/+","Dctn6-/+","Adamts14-/-","Sdf2l1-/-",
			"Tomm22-/+","E130311K13RIK-/-","1700001O22Rik-/-"))
		#Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Tbc1d5-/-","Spin1-/+","Gm572-/+",
		#	"Mast3-/-","Ofcc1-/-","Ust-/+","Tmco1-/+","Xpnpep3-/+","Ylpm1-/+","Atp5o-/+","Tomm40-/+",
		#	"Rai14-/+","Pskh1-/+","Rapgef1-/+","Pip5k1c-/+","Rec8-/-","Sesn1-/-","Ccdc28a-/-",
		#	"Slc6a20b-/-","Ints10-/+","Kcnd3-/-","Steap3-/-","Mpv17-/-","Zfp329-/-","Fam19a1-/-",
		#	"E130311K13RIK-/-","Igfl3-/-","Tada1-/+","Keg1-/-","Ppip5k2-/+","Tex2-/+","Best2-/+",
		#	"Dctn6-/+","Dalrd3-/-","Dbn1-/-","Bzw1-/-","Med10-/+","Prpf4b-/+","Tomm22-/+",
		#	"Sdf2l1-/-","Adamts14-/-","Mcmbp-/+","St6galnac3-/-","Prkcsh-/+","Rusc1-/-","Mrpl58-/+",
		#	"BC052040-/+","Cdc14a-/+","Abi3-/-","Cmtr2-/+","Nr2f1-/+","Rrs1-/+","1700001O22Rik-/-",
		#	"Clec3b-/-","Pdgfrl-/-","Slc50a1-/-","Zfp579-/-","Tox4-/+","Arhgap29-/+",
		#	"Tlk1-/-")) #Variance
		#Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ","Tmem222-/+","Alg11-/+"))
	}
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.circ)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.circ)))
	for (s in 1:length(Mutants)){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength','speed',Phenos.circ)]	
		df[names(df) %in% Phenos.circ] <- lapply(df[names(df) %in% Phenos.circ], function(x) 
			circular(x, type = 'directions', units = 'radians'))
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
		df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    	df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)],2, function(x) (x-mean(x))/sd(x))
    	if (model=='M3')
    	{X <- cbind(df$BodyLength, df$Genotype, df$speed, df$Sex)} else if (model == 'M2')
    	{X <- cbind(df$speed, df$Genotype, df$Sex)} else 
    	{X <- cbind(df$BodyLength, df$Genotype, df$Sex)}
    	lapply(Phenos.circ, function(c) print(watson.test(as.numeric(unlist(df[paste0(c)])),alpha = 0.05, dist = "vonmises")))
    	fits <- lapply(Phenos.circ, function(p) lm.circular(y=df[[paste0(p)]], x = X, init=rep(0,ncol(X)), type='c-l'))
    	pvalGen[s,] <- sapply(seq_along(fits), function(x) fits[[x]]$p.values[2]) 
		esizeGen[s,] <- sapply(seq_along(fits), function(x) fits[[x]]$coefficients) 
    }
	esizeGen <- as.matrix(esizeGen)
	rownames(pvalGen) <- Mutants
	colnames(pvalGen) <- Phenos.circ.Nomen
	rownames(esizeGen) <- Mutants
	colnames(esizeGen) <- Phenos.circ.Nomen
	pvalGenadj <- sapply(Phenos.circ.Nomen, function(r) p.adjust(pvalGen[,r],method = 'fdr'))
	pvalGenSignif <- apply(pvalGenadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodGen <- apply(pvalGenadj,2, function(x)-log(x))
    rownames(pvalGenSignif) <- Mutants
    
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.circ.Nomen
	lodGen <- lodGen[order(row.names(lodGen)),]
	esizeGen <- esizeGen[order(row.names(esizeGen)),]
	pvalGenSignif <- pvalGenSignif[order(row.names(pvalGenSignif)),]
	lodGen2 <- lodGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	esizeGen2 <- esizeGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	pvalGenSignif2 <- pvalGenSignif[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	Mutants <- rownames(lodGen2)
	row <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][1]
	col <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][2]
	esizeGen2[row,col] <- tail(sort(abs(esizeGen2)),2)[1]

	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.p <- Heatmap((lodGen2), row_names_gp = gpar(fontsize = 16, fontface='italic'),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title=paste0(model,': Variance Phenotypes'))

	col_fun <- circlize::colorRamp2(c(min(esizeGen2),0,max(esizeGen2)),c("#542788","#f7f7f7","#d73027"))
	ht.Gen.e <- Heatmap((esizeGen2), row_names_gp = gpar(fontsize = 16, fontface="italic"),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(round(min(esizeGen2)) - 1,0,round(max(esizeGen2))+0.5),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen2[i, j]) * 1/2 * (min(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e
	return(list(ht.Gen,lodGen2, esizeGen2))

}

tmp <- komp_circular(model='M3')

dev.print(pdf,'../Temp4/circular-mean-M1.pdf',width= 5.6, height=7.5)
dev.print(pdf,'../Temp4/circular-mean-M2.pdf',width= 5 , height=14)
dev.print(pdf,'../Temp5/circular-mean-M3-5.pdf',width= 4.4, height=11.4)

dev.print(pdf,'../Temp4/circular-variance-M1.pdf',width= 4.6, height=8.6)
dev.print(pdf,'../Temp4/circular-variance-M2.pdf',width= 4.6, height=8.6)
dev.print(pdf,'../Temp5/circular-variance-M3-5.pdf',width= 4.5, height=11.4)


dfp.melt <- reshape::melt(tmp[[2]], id.vars = 'Strain')
names(dfp.melt) <- c('Strain','Phenotype','LOD')
dfe.melt <- reshape::melt(tmp[[3]], id.vars = 'Strain')
names(dfe.melt) <- c('Strain', 'Phenotype', 'Esize')
df <- dfp.melt
df <- cbind(df,esize=dfe.melt$Esize)
ggplot(df, aes(x = esize, y=LOD, color=Phenotype)) + geom_point(shape = 1,size=2,stroke=1) + 
ggrepel::geom_text_repel(aes(label=ifelse((LOD >= 2.99 & abs(esize) >= 0.2),paste0(Strain), "")),box.padding = 0.5) + 
geom_vline(xintercept = c(-1,1), linetype='dashed') + geom_hline(yintercept = 2.99, linetype='dashed') + theme_bw(base_size=14) + 
theme(legend.position = 'top') + labs(x = 'Effect Size', y = 'LOD Score')
ggsave('../Temp5/volcano-var-circular-M3.pdf', width=8, height=8)

sort(abs(tmp[[3]][,3])[names(which(tmp[[2]][,3] >= 2.50))])

#Further analysis 
data_per_strain <- aggregate(x = data_per_animal[,names(data_per_animal) %in% c(Phenos.circ)], by = data_per_animal[c('Strain')], FUN = mean)
cdist <- apply(data_per_strain[,sapply(data_per_strain,is.numeric)], 2, function(x) ifelse(abs(x - median.circular(x)) <= pi,
abs(x - median.circular(x)), 2*pi - x + median.circular(x)))


Mutants <- c('Rab3gap2-/-')
	df <- data.frame()
	for (m in seq(Mutants)){
		#cat("Mutant", paste0(Mutants[m]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
		df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
		df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
		df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 2)), ]
		df1$Strain <- droplevels(df1$Strain)
		df <- rbind(df,df1)
	}
df <- df[,-which(names(df) %in% c('BodyLength','Speed'))] #Remove BodyLength
tmp <- skmeans(as.matrix(df[,names(df) %in% c(Phenos.circ)]), k = 2)