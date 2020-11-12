#Remove multivariate outliers 
gBG <- c('C57BL/6NJ') #'em1J','Hmgu','Mbp','Vlcg','Wtsi'
invisible(lapply(seq(gBG), function(x) {
	df <- data_per_animal[data_per_animal$BG %in% c(gBG),]
	tmp <- mvoutlier::dd.plot(df[,names(df) %in% Phenos.lin])
	ids <- df[tmp$outliers,'MouseID']
	data_per_animal <<- data_per_animal[-which((data_per_animal$BG == gBG) & (data_per_animal$MouseID %in% ids)),]
	data_per_stride <<- data_per_stride[-which((data_per_stride$BG == gBG) & (data_per_stride$MouseID %in% ids)),]
}))

tmp <- komp_lmer()

df_animal <- data_per_animal[data_per_animal$Strain %in% c('C57BL/6NJ',tmp[[2]]),]
df_animal$sorted_strain <- make_sort_column(df_animal, 'Strain', 'step_length1' ,FUN = median,decreasing=FALSE)
ggplot(df_animal, aes_string(x = 'sorted_strain', y = 'step_length1')) + geom_boxplot() + 
geom_point(aes(color=Sex)) + scale_color_manual(values=c("#E41A1C", "#377EB8")) + theme_bw(base_size=11) + 
theme(axis.text.x=element_text(angle=90, vjust=.5,size=11)) + labs(x = 'Strain', y = 'Step Length') 
ggsave('sw-mean.pdf',width=15,height=5)

