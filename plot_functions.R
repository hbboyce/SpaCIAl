library(cetcolor)
library(ggplot2)
library(extrafont)
library(ggpolypath)
library(dplyr)
library(bayesplot)
#font_import()
loadfonts()

ps <- 2
fs <- 24
mp <- 5
scale_min <- -3.25
scale_max <- 3.25
plot_width <- 9.66
plot_height <- 6
#figure.dir <- "/scratch/users/hboyce/CellDIVE/figures"
figure.dir <- "figures"
file_type <- "png"


PlotEpsilon <- function(feature, epsilon) {
	g <- data.frame(x=seq_along(feature), y=feature) %>%
		ggplot(aes(x=x, y=y)) +
		geom_line() +
		geom_hline(yintercept=epsilon, linetype="dashed", color="#85144B", size=2) +
		geom_text(label=sprintf("eps: %.3f", epsilon), x=0, y=epsilon + 5, hjust=0, vjust=0, size=fs*0.5, color="#85144B") +
		xlab("Index") +
		ylab("6-nn distances") +
		theme_bw() +
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.margin=margin(0,0,0,0,"cm"),
			  #title=element_blank(),
			  axis.text=element_text(size=fs),
			  legend.key.size=unit(1,"cm"),
			  legend.text=element_text(size=fs),
			  legend.text.align=1)
	
	return(g)
}

PlotHeatmap <- function(X, Y, area, marker, expression, sample_id, figure.dir=figure.dir, plot_object=FALSE) {
	
	# extract names
	sample_id <- unique(sample_id)[1]
	marker <- unique(marker)[1]
	
	# create directory and file name
	file.dir <- sprintf("%s/heatmaps/%s", figure.dir, marker)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/%s_subject_%s.%s", file.dir, marker, sample_id, file_type)
	
	# scale the expression for visualization purposes
	#expression <- scale(expression)
	expression[expression < scale_min] <- scale_min
	expression[expression > scale_max] <- scale_max
	
	# isoluminant color palette
	perceptually.uniform <- colorRampPalette(cet_pal(256, name="rainbow"))
	# plotting object
	g <- tibble(value=as.numeric(expression), X=X, Y=Y, area=area) %>%
		#tibble(value=as.numeric(expression), X=X, Y=Y) %>%
		ggplot(aes(X,Y)) +
		#geom_point(aes(color=value), size=ps) +
		geom_point(aes(color=value, size=area), alpha=0.75) +
		scale_color_gradientn(colors=perceptually.uniform(256), limits=c(scale_min,scale_max), name=sprintf("Normalized\nexpression")) +
		labs(x=expression("Cell Center X ["*mu*m*"]"), y=expression("Cell Center Y ["*mu*m*"]"),
			 title=sprintf("Sample ID: %s\nProtein: %s", sample_id, marker)) +
		scale_y_reverse() +
		coord_fixed() +
		scale_size(guide="none", range=c(1,3)) +
		theme_bw() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.title=element_blank(),
			  plot.margin=margin(0,0,0,0,"cm"),
			  #title=element_blank(),
			  axis.text=element_text(size=fs),
			  legend.key.size=unit(1,"cm"),
			  legend.text=element_text(size=fs),
			  legend.text.align=1,
			  #axis.title=element_blank(),
			  panel.grid=element_blank())
	
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)	
		ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)	
	}
	
	# save plot to pdf
	#CairoPDF(fn, width=8, height=6, bg="transparent", family="Arial")
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	# pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	# print(g)
	# dev.off()
}

### Plot the single cell cutoffs for hypoxia
PlotSingleCell <- function(X, Y, area, sample_id, thresholded, figure.dir=figure.dir, plot_object=FALSE) {
	# extract names
	sample_id <- unique(sample_id)[1]
	thresholded <- factor(thresholded, levels=c(TRUE, FALSE), labels=c("Hypoxic", "Normoxic"))
	
	# save file name
	file.dir <- sprintf("%s/regions/singlecell", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/singlecell_subject_%s.%s", file.dir, sample_id, file_type)
	
	# specifify the color of points in each group
	cell.colors <- c(Hypoxic="#85144B", Normoxic="#cccccc")
	
	# plot
	g <- tibble(value=thresholded, X=X, Y=Y, area=area) %>%
		ggplot(aes(X,Y)) +
		#geom_point(aes(color=value), size=ps) +
		geom_point(aes(color=value, size=area), alpha=0.75) +
		scale_color_manual(values=cell.colors, name="Region") +
		labs(x=expression("Cell Center X ["*mu*m*"]"), y=expression("Cell Center Y ["*mu*m*"]"),
			 title=sprintf("Sample ID: %s", sample_id)) +
		scale_y_reverse() +
		coord_fixed() +
		scale_size(guide="none", range=c(1,3)) +
		theme_bw() +
		guides(color=guide_legend(override.aes=list(size=ps*mp))) +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.title=element_blank(),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  #axis.title=element_blank(),
			  panel.grid=element_blank())
	
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)
		ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	}
	
	# write plot to pdf
	#pdf(fn, width=8, height=6, bg="transparent")
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
	#return(g)
}


### plot DBSCAN clusters
PlotClusters <- function(X, Y, area, sample_id, cluster, cell.colors, figure.dir=figure.dir, plot_object=FALSE) {
	# extract names
	sample_id <- unique(sample_id)[1]
	cluster[cluster == 0] <- "Normoxic"
	cluster <- factor(cluster)
	
	# save file name
	file.dir <- sprintf("%s/regions/clusters", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/clusters_subject_%s.%s", file.dir, sample_id, file_type)
	
	# plot
	g <- tibble(value=cluster, X=X, Y=Y, area=area) %>%
		ggplot(aes(X,Y)) +
		#geom_point(aes(color=value), size=ps) +
		geom_point(aes(color=value, size=area), alpha=0.75) +
		scale_color_manual(values=cell.colors, name="Cluster") +
		labs(x=expression("Cell Center X ["*mu*m*"]"), y=expression("Cell Center Y ["*mu*m*"]"),
			 title=sprintf("Sample ID: %s", sample_id)) +
		scale_y_reverse() +
		coord_fixed() +
		scale_size(guide="none", range=c(1,3)) +
		theme_bw() +
		guides(color=guide_legend(override.aes=list(size=ps*mp))) +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.title=element_blank(),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  #axis.title=element_blank(),
			  panel.grid=element_blank())
	
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)
		ggsave(fn, g, width=plot_width, height=plot_height, device="png")
	}
	# write plot to pdf
	#pdf(fn, width=8, height=6, bg="transparent")
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	# pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	# print(g)
	# dev.off()
	#return(g)
}


### plot regions of hypoxia
PlotRegions <- function(X, Y, area, sample_id, region, shape, figure.dir=figure.dir, plot_object=FALSE) {
	# extract names
	sample_id <- unique(sample_id)[1]
	# if (!(sample_id %in% shapes$sample_id)) {
	# 	return(0)
	# }
	shape <- shape[[1]]
	region <- factor(region, levels=c(TRUE, FALSE), labels=c("Hypoxic", "Normoxic"))
	
	# save file name
	file.dir <- sprintf("%s/regions/regions", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/regions_subject_%s.%s", file.dir, sample_id, file_type)
	
	# specifify the color of points in each group
	cell.colors <- c(Hypoxic="#85144B", Normoxic="#cccccc")
	
	# generate data for plot
	point.data <- tibble(value=region, X=X, Y=Y)
	sid <- sample_id
	#shape.data <- fortify((shapes %>% filter(sample_id == sid) %>% pull(shape))[[1]])
	shape.data <- fortify(shape)
	
	# plot
	g <- ggplot() +
		#geom_point(data=point.data, aes(x=X, y=Y, color=value), size=ps) +
		geom_point(data=point.data, aes(x=X, y=Y, color=value, size=area), alpha=0.75) +
		scale_color_manual(values=cell.colors, name="Region") +
		ggpolypath::geom_polypath(data=shape.data, aes(y=lat, x=long, group=group), fill="#85144B", color="black", size=ps, alpha=0.40) +
		labs(x=expression("Cell Center X ["*mu*m*"]"), y=expression("Cell Center Y ["*mu*m*"]"),
			 title=sprintf("Sample ID: %s", sample_id)) +
		scale_y_reverse() +
		#scale_x_reverse() +
		scale_size(guide="none", range=c(1,3)) +
		coord_fixed() +
		theme_bw() +
		guides(color=guide_legend(override.aes=list(size=ps*mp))) +
		#coord_flip() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.title=element_blank(),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  #axis.title=element_blank(),
			  panel.grid=element_blank())
	
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)
		ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	}
	
	# write plot to pdf
	#pdf(fn, width=8, height=6, bg="transparent")
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	# pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	# print(g)
	# dev.off()
	#return(g)
}


PlotTrace <- function(posterior, parameter) {
	file.dir <- sprintf("%s/bayes/trace", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/trace_%s.%s", file.dir, tolower(parameter), file_type)
	
	p <- strsplit(parameter, "_")[[1]]
	if (p[2] == "A") {
		hist <- "Adenocarcinoma"
	} else if (p[2] == "S") {
		hist <- "Squamous Cell Carcinoma"
	} else {
		hist <-  "Other"
	}
	
	color_scheme_set("mix-blue-red")
	g <- mcmc_trace(posterior, pars=c(parameter))
	g <- g + labs(x="Iteration", y=expression(beta[j]),
			 title=(sprintf("Trace for %s in %s", p[1], hist))) +
		guides(color=guide_legend(override.aes=list(size=5))) +
		ylim(-2,2) +
		theme_bw() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  panel.grid=element_blank())
	
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
}


PlotPosteriorArea <- function(posterior, histology) {
	file.dir <- sprintf("%s/bayes/area", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/area_%s.%s", file.dir, tolower(strsplit(histology, " ")[[1]][1]), file_type)
	
	pars <- dimnames(posterior)[[3]]
	dimnames(posterior)[[3]] <- purrr::reduce(lapply(strsplit(pars, "_"), "[[", 1), c)
	print(pars)
	
	means <- c()
	for (i in 1:dim(posterior)[3]) {
		means <- c(means, mean(posterior[,,i]))
	}
	or <- order(means, decreasing=TRUE)
	posterior <- posterior[,,or]
	
	g <- mcmc_areas(posterior, prob=0.8) +
		labs(x="log2 Fold Change", title=sprintf("Posterior Densities for %s", histology)) +
		xlim(-2, 2) +
		theme_bw() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  panel.grid.minor.x=element_blank(),
			  panel.grid.major.x=element_blank())
	
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
}

PlotRhat <- function(fit, histology) {
	
	file.dir <- sprintf("%s/bayes/rhat", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/rhat.%s", file.dir, file_type)
	
	color_scheme_set("brightblue")
	
	rhats <- rhat(fit)
	rhats <- rhats[grepl("region", names(rhats)) & !grepl("Sigma", names(rhats))]
	
	region_pars <- names(rhats)
	#histology <- substr(lapply(strsplit(region_pars, ":"), "[[", 4), 1, 1)
	protein <- purrr::reduce(lapply(strsplit(region_pars, ":"), "[[", 2), c)
	pars <- paste(protein, histology, sep="_")
	
	names(rhats) <- pars
	
	g <- mcmc_rhat(rhats, size=2) +
		yaxis_text(hjust=1) +
		theme_bw() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  panel.grid.minor.x=element_blank(),
			  panel.grid.major.x=element_blank())
	
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
	
}

PlotNeff <- function(fit, histology) {
	
	file.dir <- sprintf("%s/bayes/neff", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/neff.%s", file.dir, file_type)
	
	color_scheme_set("brightblue")
	
	neff <- neff_ratio(fit)
	neff <- neff[grepl("region", names(neff)) & !grepl("Sigma", names(neff))]
	
	region_pars <- names(neff)
	#histology <- substr(lapply(strsplit(region_pars, ":"), "[[", 4), 1, 1)
	protein <- purrr::reduce(lapply(strsplit(region_pars, ":"), "[[", 2), c)
	pars <- paste(protein, histology, sep="_")
	
	names(neff) <- pars
	
	g <- mcmc_neff(neff, size=2) +
		yaxis_text(hjust=1) +
		theme_bw() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  panel.grid.minor.x=element_blank(),
			  panel.grid.major.x=element_blank())
	
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
}

PlotAcf <- function(posterior, parameter) {
	file.dir <- sprintf("%s/bayes/acf", figure.dir)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/acf_%s.%s", file.dir, tolower(parameter), file_type)
	
	p <- strsplit(parameter, "_")[[1]]
	i <- length(p)
	if (p[i] == "A") {
		hist <- "Adenocarcinoma"
	} else if (p[i] == "S") {
		hist <- "Squamous Cell Carcinoma"
	} else {
		hist <-  "Other"
	}
	
	color_scheme_set("brightblue")
	g <- mcmc_acf(posterior, pars=c(parameter))
	g <- g + labs(title=(sprintf("Autocorrelation for %s in %s", p[1], hist))) +
		guides(color=guide_legend(override.aes=list(size=5))) +
		theme_bw() +
		#theme(text=element_text(size=fs, family="Arial"),
		theme(text=element_text(size=fs),
			  title=element_text(size=fs),
			  plot.margin=margin(0,0,0,0,"cm"),
			  axis.text=element_text(size=fs),
			  legend.text=element_text(size=fs),
			  legend.text.align=0,
			  strip.text.x=element_blank(),
			  panel.grid=element_blank())
	
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent", family="Arial")
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
}

PlotDifferentialCorrelation <- function(df, subtype="", data_source="GE_Lung") {
	
	file.dir <- sprintf("%s/%s/differential", figure.dir, data_source)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/diffcorr_%s.%s", file.dir, subtype, file_type)
	
	rownames(df) <- NULL
	cDiff <- df$cdiff
	Norm <- df$norm
	
	# color palette for plot
	color.palette <- colorRampPalette(c("#0072B2", "#56B4E9", "white", "#E69F00", "#E34234"))(n = 500)
	breaks.hm <- seq(-1, 1, 0.01)
	
	color_box <- colorRampPalette(cet_pal(5, name="d2"))(n=500)
	
	g <- #ggplot(df, aes(marker1, ordered(marker2, levels=rev(levels(marker2))))) +
		ggplot(df, aes(marker1, forcats::fct_rev(factor(marker2)))) +
		#geom_tile(color="black", fill="white") +
		geom_tile(aes(fill=norm), color="black", alpha=.9) +
		geom_point(aes(size=sds, color=cdiff, alpha=overlap)) +#, stroke=0) +
		geom_point(aes(size=sds, alpha=overlap), shape=1, color="black") +
		#geom_point(data=df[is.na(df$sds),], color="white", size=9, shape=15) +
		scale_alpha(range=c(0,1)) + 
		scale_color_gradientn(colors=color.palette, breaks=c(-7,0,7), limits=c(-7,7), na.value="white") +
		scale_fill_gradientn(colors=color_box, breaks=c(-1,0,1), limits=c(-1,1), na.value="white") +
		#scale_y_discrete(limits=rev(df$marker1)) +
		scale_size(range=c(0,10)) +
		guides(alpha=FALSE, size=FALSE) +
		coord_fixed(ratio=1) +
		theme(text=element_text(size=fs),#, family="Arial"),
			  legend.position="none",
			  axis.text.x=element_text(angle=90, size=fs, hjust=1, vjust=0.5),
			  axis.text.y=element_text(size=fs),
			  plot.title=element_text(size=fs),
			  axis.title=element_blank(),
			  legend.key.size=unit(2,"cm"),
			  legend.title=element_blank(),
			  legend.text=element_text(size=fs),
			  panel.grid.major=element_blank(),
			  panel.grid.minor=element_blank(),
			  panel.background=element_rect(fill="transparent", color=NA),
			  plot.background=element_rect(fill="transparent", color=NA))
	
	# create data.frame for the color spectrum
	#hist.data <- seq(-7,7,0.1)
	#df <- data.frame(x=hist.data, y=1)
	
	# create data.frame for the density plot overlay
	#df1 <- data.frame(x=cDiff)
	
	r <- 7
	
	df_circ <- expand.grid(x=seq(-r,r,0.1), y=seq(-r,r,0.1)) %>%
		filter(x**2 + y**2 <= r**2)
	dens <- density(cDiff)
	df1 <- data.frame(x=dens$x, y=dens$y*20) %>%
		filter(abs(x) < 7)
	p <- df_circ %>%
		ggplot() +
		geom_tile(aes(x=x, y=y, fill=x)) + 
		scale_fill_gradientn(colors=color.palette, breaks=c(-7,0,7), limits=c(-7,7), na.value="white") + 
		geom_ribbon(data=df1, aes(x=x, ymax=y-3), ymin=-3, alpha=0.5, color="black", fill="black") + 
		theme_bw() +
		coord_fixed() + 
		labs(title="Differential \nCorrelation (dz)", x="dz") + 
		theme(text=element_text(size=fs-4),
			axis.text.x=element_text(angle=45, size=fs-4, hjust=1),# vjust=0.5),
			axis.text.y=element_blank(),
			axis.line=element_blank(),
			axis.ticks.y=element_blank(),
			axis.ticks.length=unit(0.25,"cm"),
			panel.background=element_blank(),
			axis.title=element_blank(),
			plot.background=element_blank(),
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			panel.border=element_blank(),
			legend.position="none")
	
	# Legend plot
	# p <- ggplot() +
	# 	geom_area(data=df, aes(x=x, y=y)) +
	# 	#geom_bar(data=df, aes(x=x, y=y, fill=x), stat="identity", show.legend=FALSE, position=position_dodge2(padding=0)) +
	# 	scale_y_continuous(expand=c(0,0)) +
	# 	scale_fill_gradientn(colors=color.palette, breaks=c(-7,0,7), limits=c(-7,7), na.value="white") +
	# 	geom_density(data=df1, aes(x=x, y=..scaled..), fill="black", alpha=0.4) + 
	# 	xlim(c(-7.1,7.1)) +
	# 	theme_bw() +
	# 	coord_fixed(ratio=2) +
	# 	labs(title="Differential Correlation") +
	# 	theme(#text=element_text(size=fs, family="Arial"),
	# 		text=element_text(size=fs),
	# 		  axis.text.x=element_text(angle=0, size=fs),#, hjust=1, vjust=0.5),
	# 		  axis.text.y=element_blank(),
	# 		  axis.line=element_blank(),
	# 		  axis.ticks.y=element_blank(),
	# 		  axis.ticks.length=unit(0.25,"cm"),
	# 		  panel.background=element_blank(),
	# 		  axis.title=element_blank(),
	# 		  plot.background=element_blank(),
	# 		  panel.grid.major=element_blank(),
	# 		  panel.grid.minor=element_blank(),
	# 		  panel.border=element_blank())
	# 
	
	#df_box <- expand.grid(x=seq(-1,1,0.01), y=seq(-1,1,0.01))
	df_box <- expand.grid(x=seq(-1,1,0.01), y=seq(-0.57,1.43,0.01))
	df2 <- data.frame(x=Norm)
	color.palette <- colorRampPalette(c("#0072B2", "#56B4E9", "white", "#E69F00", "#E34234"))(n = 500)
	q <- df_box %>%
		ggplot() +
		geom_tile(aes(x=x, y=y, fill=x)) +
		scale_fill_gradientn(colors=color_box, breaks=c(-1,0,1), limits=c(-1,1), na.value="white") + 
		geom_density(data=df2, aes(x=x, y=..scaled..), fill="black", alpha=0.5) +
		scale_x_continuous(n.breaks=3) +
		theme_bw() +
		coord_fixed() +
		labs(title="Normoxic \ncorrelation", x="Correlation") +
		theme(text=element_text(size=fs-4),
			axis.text.x=element_text(angle=45, size=fs-4, hjust=1),#, vjust=0.25),
			axis.text.y=element_blank(),
			axis.line=element_blank(),
			axis.ticks.y=element_blank(),
			axis.ticks.length=unit(0.25,"cm"),
			panel.background=element_blank(),
			axis.title=element_blank(),
			plot.background=element_blank(),
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			panel.border=element_blank(),
			legend.position="none")
	
	
	# hist.data <- seq(-1,1,0.01)
	# df <- data.frame(x=hist.data, y=1)
	# df2 <- data.frame(x=Norm)
	# q <- ggplot() +
	# 	geom_bar(data=df, aes(x=x, y=y, fill=x), stat="identity", show.legend=FALSE, position=position_dodge2(padding=0)) +
	# 	scale_y_continuous(expand=c(0,0)) +
	# 	scale_fill_gradientn(colors=color.palette, breaks=c(-1,0,1), limits=c(-1,1), na.value="white") +
	# 	geom_density(data=df2, aes(x=x, y=..scaled..), fill="black", alpha=0.4) + 
	# 	xlim(c(-1.015,1.015)) +
	# 	theme_bw() +
	# 	coord_fixed(ratio=2/7) +
	# 	labs(title="Normoxic Correlation") +
	# 	theme(#text=element_text(size=fs, family="Arial"),
	# 		text=element_text(size=fs),
	# 		axis.text.x=element_text(angle=0, size=fs),#, hjust=1, vjust=0.5),
	# 		axis.text.y=element_blank(),
	# 		axis.line=element_blank(),
	# 		axis.ticks.y=element_blank(),
	# 		axis.ticks.length=unit(0.25,"cm"),
	# 		panel.background=element_blank(),
	# 		axis.title=element_blank(),
	# 		plot.background=element_blank(),
	# 		panel.grid.major=element_blank(),
	# 		panel.grid.minor=element_blank(),
	# 		panel.border=element_blank())
	
	
	# place legend on correlation plot
	plt <- g +
		annotation_custom(grob=ggplotGrob(p), xmin=11, xmax=19, ymin=20, ymax=29) +
		annotation_custom(grob=ggplotGrob(q), xmin=20, xmax=28, ymin=20, ymax=29)
		#annotation_custom(grob=ggplotGrob(p), xmin=11, xmax=17, ymin=22, ymax=28) +
		#annotation_custom(grob=ggplotGrob(q), xmin=19, xmax=25, ymin=22, ymax=28)
		
	plot_width<-12
	plot_height<-12
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#ggsave(filename=fn, plot=plt, device="pdf", width=plot_width, height=plot_height)#, bg="transparent")
	ggsave(fn, plt, width=plot_width, height=plot_height, device=file_type)
	#print(g)
	#dev.off()
}

PlotDiagnostics <- function(res_fit, analysis="corr", histology, marker1, marker2=NULL, data_source="GE_Lung") {
	file.dir <- sprintf("%s/%s/mcmc_diagnostics_%s/%s", figure.dir, data_source, analysis, histology)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	
	
	# plot neff_ratios
	g <- mcmc_neff(res_fit$neff) + yaxis_text(TRUE)
	
	if (is.null(marker2)) {
		fn <- sprintf("%s/%s_%s.%s", file.dir, "neff", marker1, file_type)
	} else {
		fn <- sprintf("%s/%s_%s_%s.%s", file.dir, "neff", marker1, marker2, file_type)
	}
	
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
	
	# plot rhat
	g <- mcmc_rhat(res_fit$rhat) + yaxis_text(TRUE)
	if (is.null(marker2)) {
		fn <- sprintf("%s/%s_%s.%s", file.dir, "neff", marker1, file_type)
	} else {
		fn <- sprintf("%s/%s_%s_%s.%s", file.dir, "neff", marker1, marker2, file_type)
	}
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
	
	
	# extract chains
	chains <- res_fit$chains
	colnames(chains) <- 1:4
	
	# plot acf
	auto <- apply(chains, 2, function(chain) as.double(acf(chain, lag.max=20, plot=FALSE)$acf)) %>%
		as_tibble() %>%
		gather("Chain", "Autocorrelation")
	
	
	g <- auto %>%
		group_by(Chain) %>%
		mutate(Lag=seq_along(Autocorrelation)) %>%
		ungroup() %>%
		ggplot() + 
		geom_path(aes(Lag, Autocorrelation), color="#0072B2") +
		facet_grid(Chain ~ .) +#, switch="y") +
		geom_hline(yintercept=0, alpha=0.5, color="#0072B2") +
		geom_segment(aes(x=Lag, xend=Lag, yend=Autocorrelation), y=0, size=0.25, alpha=0.25, color="#0072B2") +
		ylim(-1.1, 1.1) +
		bayesplot_theme_get() +
		theme(strip.text.y = element_text(angle=0)) 
	
	if (is.null(marker2)) {
		fn <- sprintf("%s/%s_%s.%s", file.dir, "neff", marker1, file_type)
	} else {
		fn <- sprintf("%s/%s_%s_%s.%s", file.dir, "neff", marker1, marker2, file_type)
	}
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
	
	# plot trace
	g <- chains %>%
		as_tibble() %>%
		gather("Chain", "mu") %>%
		group_by(Chain) %>%
		mutate(Iteration = seq_along(mu)) %>%
		ungroup() %>%
		ggplot(aes(x=Iteration, y=mu, color=Chain)) +
		geom_path(alpha=0.5) +
		geom_hline(yintercept=0, alpha=0.5) +
		scale_color_brewer(palette="Paired") +
		bayesplot_theme_get()
	
	if (is.null(marker2)) {
		fn <- sprintf("%s/%s_%s.%s", file.dir, "neff", marker1, file_type)
	} else {
		fn <- sprintf("%s/%s_%s_%s.%s", file.dir, "neff", marker1, marker2, file_type)
	}
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
	
}

PlotDifferentialExpression <- function(df, histology, data_source="GE_Lung") {
	
	file.dir <- sprintf("%s/%s/differential", figure.dir, data_source)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/diffexpr_%s.%s", file.dir, histology, file_type)
	
	df <- df %>%
		group_by(protein) %>%
		mutate(overlap=ifelse(prod(range(x)) <= 0, TRUE, FALSE)) %>%
		ungroup()
	
	g <- df %>%
		mutate(protein=forcats::fct_reorder(protein, desc(median))) %>%
		ggplot() +
		#geom_ribbon(aes(x=x, ymax=y), ymin=0, alpha=0.5, color="#0072B2", fill="#0072B2") +
		geom_ribbon(aes(x=x, ymax=y, color=overlap, fill=overlap), ymin=0, alpha=0.5) +
		scale_color_manual(values=c("#0072B2", "#D5D8DC")) +
		scale_fill_manual(values=c("#0072B2", "#D5D8DC")) +
		geom_segment(aes(x=median, xend=median, yend=yend), y=0, size=1) +
		geom_vline(xintercept=0, size=1.5) +
		labs(fill="", x="Log 2 Fold Change") +
		guides(color=FALSE, alpha=FALSE, fill=FALSE) +
		xlim(-1.25, 1.25) +
		geom_hline(yintercept=0, alpha=0.25) +
		geom_vline(xintercept=c(-1,-0.5,0,0.5,1), alpha=0.25) +
		ylim(0,1.5) +
		facet_grid( protein ~ ., switch="y") + 
		theme(panel.background=element_blank(),
			  panel.spacing=unit(0, "lines"),
			  legend.background=element_blank(),
			  legend.text=element_text(size=fs),
			  legend.key.size=unit(1,"cm"),
			  axis.title.y=element_blank(),
			  axis.ticks.y=element_blank(),
			  axis.text.y=element_blank(),
			  axis.text.x=element_text(size=fs),
			  axis.title.x=element_text(size=fs),
			  strip.background=element_blank(),
			  strip.text.y.left=element_text(angle=0, size=fs, hjust=1, vjust=0.5))
	
	plot_width<-12
	plot_height<-12
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	#pdf(fn, width=plot_width, height=plot_height, bg="transparent")
	#print(g)
	#dev.off()
}

PlotRegionCount <- function(n_shapes, plot_object=FALSE, data_source="GE_Lung") {
	
	file.dir <- sprintf("%s/%s/regions_desc", figure.dir, data_source)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/sample_count.%s", file.dir, file_type)
	
	bins <- max(n_shapes$n_shapes)
	g <- n_shapes %>%
		ggplot(aes(n_shapes)) +
		geom_histogram(bins=bins+1, color="black", alpha=0.5) +
		scale_x_continuous("Number of hypoxic regions", labels=as.character(seq(0,bins,4)), breaks=seq(0,bins,4)) +
		ylab("Sample Count") +
		theme_bw() +
		theme(panel.grid=element_blank(),
			  axis.text=element_text(size=fs*1.5),
			  axis.title=element_text(size=fs*1.5))
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height)
		ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	}
}

PlotCircularity <- function(dta, plot_object=FALSE, data_source="GE_Lung") {
	
	file.dir <- sprintf("%s/%s/regions_desc", figure.dir, data_source)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/circularity.%s", file.dir, file_type)
	
	ra <- as.integer(range(dta$perimeter))
	x <- ra[1]:ra[2]
	a <- (x**2)/(4*pi)
	circ <- log(a/x)
	circ <- data.frame(x=x, circ=circ)
	p <- 2*pi*500
	hline <- data.frame(p=as.integer(p):ra[2], hline=log((pi*500**2)/p))
	vline <- data.frame(p=p, vline=seq(min(dta$circularity), log(((p**2)/(4*pi))/p), 0.01))
	
	g <- dta %>%
		ggplot(aes(perimeter, circularity)) +
		geom_point(pch=21, fill="black", alpha=0.5, size=3) +
		geom_line(data=circ, aes(x, circ), size=2) +
		#geom_line(data=hline, aes(p, hline), size=2) +
		#geom_line(data=vline, aes(p, vline), size=2) +
		labs(x="Region Perimeter (\u03bcm)", y="Circularity") +
		theme_bw() +
		theme(panel.grid=element_blank(),
			  axis.text=element_text(size=fs*1.5),
			  axis.title=element_text(size=fs*1.5))
	
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)
		ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	}
}

PlotProportion <- function(dta, plot_object=FALSE, data_source="GE_Lung") {
	
	file.dir <- sprintf("%s/%s/regions_desc", figure.dir, data_source)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/proportion.%s", file.dir, file_type)
	
	g <- ggplot(dta, aes(proportion)) +
		geom_histogram(bins=50, color="black", alpha=0.5) +
		labs(x="Hypoxic proportion", y="Count") +
		theme_bw() +
		theme(panel.grid=element_blank(),
			  axis.text=element_text(size=fs*1.5),
			  axis.title=element_text(size=fs*1.5))
	
	if (plot_object) {
		return(g)
	} else {
		#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)
		ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
	}
	
}

PlotKaplanMeier <- function(g, data_source="GE_Lung", type="") {
	file.dir <- sprintf("%s/%s/survival", figure.dir, data_source)
	dir.create(file.dir, recursive=TRUE, showWarnings=FALSE)
	fn <- sprintf("%s/km_%s.%s", file.dir, type, file_type)
	
	#ggsave(fn, g, width=plot_width, height=plot_height, device=cairo_pdf)
	ggsave(fn, g, width=plot_width, height=plot_height, device=file_type)
}
