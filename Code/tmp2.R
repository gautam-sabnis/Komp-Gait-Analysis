#Filter Strains for which at least 8 animals were tested 
Strains8 <- names(table(data_per_animal$Strain))[table(data_per_animal$Strain) >= 5]
data_per_animal <- data_per_animal[data_per_animal$Strain %in% Strains8,] 
data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains8,]
data_per_animal$Strain <- droplevels(data_per_animal$Strain)
data_per_stride$Strain <- droplevels(data_per_stride$Strain)

df <- table(data_per_animal$Strain, data_per_animal$Sex)
df <- cbind(Strain = rownames(df), df[,1:2], B = sapply(sort(unique(data_per_animal$Strain)), function(x) length(unique(data_per_animal[data_per_animal$Strain == x, 'TestDate'])))) 
rownames(df) <- NULL
df <- as.data.frame(df)
df1 <- df[1:45,]
df1 <- df1[-14,]
df2 <- df[46:89,]
df3 <- df[90:133,]
dff <- cbind(df1,df2,df3)


tmp <- unique(df.out[df.out$Outlier == 1, 'Strain'])
tmp2 <- setdiff(unique(data_per_animal$Strain), Strains8)[which(setdiff(unique(data_per_animal$Strain), Strains8) %in% tmp)]


length(intersect(Strains1,mvoutStrains))/length(mvoutStrains)

df <- data.frame(Analysis = c('M1','M2','M3','Mvoutlier'), Proportion = c(length(intersect(Strains1,mvoutStrains))/length(mvoutStrains),
	length(intersect(Strains2,mvoutStrains))/length(mvoutStrains),length(intersect(Strains3,mvoutStrains))/length(mvoutStrains),
	length(intersect(mvoutStrains,mvoutStrains))/length(mvoutStrains)), Total = c(length(Strains1), length(Strains2), length(Strains3),
	length(mvoutStrains))) 
ggplot(df, aes(x=Analysis,y=Proportion)) + geom_bar(stat='identity') + theme_bw(base_size=22) +  
theme( 
	axis.text.x.top = element_text(vjust=0.5)) + labs(y = 'Proportion') + 
geom_text(label = paste0(df$Total),vjust = -.01, hjust = 0.5,size=7)
ggsave('../Temp6/pinkStrains.pdf', width=9, height=9)


##

mvoutStrains <- gsub("*./.*","",mvoutStrains)

AllStrains <- union(union(union(Strains1,Strains2),Strains3),mvoutStrains)

df <- data.frame(Strains = AllStrains, M1 = rep(0,length(AllStrains)), M2 = rep(0,length(AllStrains)), M3 = rep(0,length(AllStrains)), 
	Mvout = rep(0,length(AllStrains)))

df$M1 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains1, 1, 0))
df$M2 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains2, 1, 0))
df$M3 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains3, 1, 0))
df$Mvout <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% mvoutStrains, 1, 0))

UpSetR::upset(df, text.scale = 4, point.size = 4, line.size = 2)



phenList <- read.csv('../Data/phenotypeHitsPerParameterAndProcedure.csv', header = TRUE)
allstats <- read.csv('../Data/statistical-results-ALL.csv', header = TRUE)
parameter_id <- phenList[phenList$Parameter.Name %in% c('Gait','Gait (inc. ataxia)'),'Parameter.Id']
df0 <- allstats[allstats$parameter_stable_id %in% parameter_id, c('parameter_name','marker_symbol','p_value')]

cutoff <- 0.05
df <- df0[df0$p_value < cutoff,]
#df$pipeline_name <- droplevels(df$pipeline_name)

tmp <- data.frame(Pipeline = names(table(df$pipeline_name)), Proportion = as.numeric(table(df$pipeline_name))/as.numeric(table(df0$pipeline_name)), 
	Total = as.numeric(table(df0$pipeline_name)))
tmp <- tmp[complete.cases(tmp), ]
tmp <- tmp[-c(6,9),]
tmp$Pipeline <- c('BCM','CCP','German Mouse Clinic','Harwell','Harwell Interval','ICS',
	'IMPC','JAX','MGP','TCP','UCD')

tmp2 <- data.frame(Pipeline = rep(tmp$Pipeline,each=2), 
	Percent = as.numeric(sapply(seq(tmp$Pipeline), function(x) c(tmp[x,2], 100 - tmp[x,2]))), 
	Condition = rep(c('Gait', 'No Gait'), length = 22))


ggplot(tmp2, aes(fill = Condition, y = Percent, x = Pipeline)) + geom_bar(position = 'fill', stat = 'identity') + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
	axis.text.x.top = element_text(vjust=0.5))
tmp$Pipeline <- with(tmp, reorder(Pipeline,-Proportion))

p2 <- ggplot(tmp, aes(x = Pipeline, y = Proportion)) + geom_bar(position = 'stack',stat = 'identity') + theme_bw(base_size = 22) + 
geom_text(label = paste0(tmp$Total),vjust = -.1, hjust = 0.5,size=7.5) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none',
	axis.text.x.top = element_text(vjust=0.5)) + ggtitle(paste0('p-value < ',cutoff))

data_per_animal$Strain <- gsub("*./.*","",data_per_animal$Strain)
common_lines <- unique(data_per_animal$Strain)[unique(data_per_animal$Strain) %in% df$marker_symbol]
sapply(seq(common_lines), function(x) nrow(data_per_animal[data_per_animal$Strain == common_lines[x],]))

################################################################
setwd("/Users/sabnig/Documents/Projects/Komp/Temp")
#experiment <- read.csv('../Data/experiment.csv',header=TRUE)
df0 <- read.csv('../Data/statistical-result.csv', header=TRUE)
cutoff <- 0.05
df <- df0[df0$p_value < cutoff,names(df0) %in% c('phenotyping_center','marker_symbol','p_value')]
df <- df[complete.cases(df),]
rownames(df) <- 1:nrow(df) 

tmp <- data.frame(Center = names(table(df$phenotyping_center)), Proportion = as.numeric(table(df$phenotyping_center))/as.numeric(table(df0$phenotyping_center)), 
	Total = as.numeric(table(df0$phenotyping_center)))
tmp$Center <- with(tmp,reorder(Center,-Proportion))

ggplot(tmp, aes(x = Center, y = Proportion)) + geom_bar(stat = 'identity') + 
geom_text(label = paste0(tmp$Total),vjust = -.1, hjust = 0.5,size=7.5) + theme_bw(base_size = 22) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none',
	axis.text.x.top = element_text(vjust=0.5)) + ggtitle(paste0('p-value < ',cutoff))

mutants_gait <- sort(df$marker_symbol[complete.cases(df$marker_symbol)]) 