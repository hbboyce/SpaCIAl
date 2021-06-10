# require(plyr)
# require(raster)
# require(rgeos)
# 
# library(Rcpp)
# library(dplyr)
# library(dbscan)
# library(geometry) # delaunayn
# library(prevR) # points in spatial polygons
# library(sp)



#' Get the knee of a set of feature values
#' 
#' \code{GetKnee} returns the value for the knee of a set of passed in values.
#' 
#' This function sorts the passed in feature in increasing order and fits a locally polynomial regression (improves curvature
#' estimation) and calculates the maximum curvature of the fitted lines using finite differences.
#' 
#' @param feature numeric vector. 
#' @return knee. the value where curvature is maximized for the local polynomial regression fit of the feature.
GetKnee <- function(feature, plot_object=FALSE) {
	feature <- sort(feature)
	x <- seq_along(feature)
	
	# smooth data to get accurate curvature calculation
	lo <- stats::loess(feature ~ x, span=5)
	d <- stats::predict(lo, x)
	
	# calculate the curvature
	first.order <- unlist(sapply(c(2:length(d)-1), function(i) (d[i+1] - d[i-1])/2))
	second.order <- unlist(sapply(c(2:length(d)-1), function(i) (d[i+1] - 2*d[i] + d[i-1])/4))
	curvature <- second.order/((1 + first.order^2)^1.5)
	
	# knee is where curvature is maximised
	knee <- d[which.max(curvature)]
	
	g <- PlotEpsilon(feature, knee)
	
	if (plot_object) {
		return(g)
	} else {
		return(knee)	
	}
	
}


#' Determine density parameter epsilon
#' 
#' \code{GetEpsilon} returns a scalar value (epsilon) that is a proxy for the density of a two-dimensional point set
#' 
#' This function groups points in the point set into tiles based on their position, selects the median density tile 
#' based on number of cells per tile, and calculates a quantity epsilon based on the (minPts-1)-nearest neighbors.
#' 
#' Specifically, the nearest neighbor distances are sorted and the curvature is calculated using first- and second-order
#' finite differences. Epsilon is chosen to be the distance where curvature is maximized. This results in a parameter whose
#' value does not drastically affect the number of points included in the density estimation, but is not completely
#' insensitive to the number of points affected.
#' 
#' @param X numeric vector. first dimension coordinate
#' @param Y numeric vector. second dimension coordinate
#' @param RANGE numeric. the maximum range of the X and Y coordinates. Assumes the same range for both X and Y and that the
#' minumum range is 0 
#' @param TILE numeric. how many regions to tile the point space into across one dimension. In two-dimensions the number of
#' tiles will be TILE**2
#' @param minPts numeric. the minimum number of points within a radius for a point to be considered a dense region
#' @param plot_object logical. Whether to return the k-nn/epsilon plot instead of the epsilon value.
#' @return Epsilon. A numeric value that represents the radius that corresponds to a region that, on average, encapsulates
#' at least minPts number of points in the point set.
#' @export
GetEpsilon <- function(X, Y, RANGE=NULL, TILE=5, minPts=7, plot_object=FALSE) {
	if (is.null(RANGE)) {
		RANGE <- diff(range(c(X,Y)))
	}

	distances <- c()

	# find which tile has the median cells
	for (x in seq(RANGE/TILE, RANGE, RANGE/TILE)) {
		for (y in seq(RANGE/TILE, RANGE, RANGE/TILE)) {
			tile <- x - RANGE/TILE < X &
				X <= x &
				y - RANGE/TILE < Y &
				Y <= y
			if (sum(tile) >= 25) {
				tile_x <- X[tile]
				tile_y <- Y[tile]
				distances <- c(distances, as.vector(dbscan::kNNdist(data.frame(X=tile_x, Y=tile_y), k=minPts-1)))
			}

		}
	}
	if (length(distances) < 25) {
		return(NA)
	}
	epsilon <- GetKnee(distances, plot_object)

	return(epsilon)

}


#' Normalize features across groups
#' 
#' \code{FeatureNormalization} returns a matrix that contains the passed in feature values scaled across groups
#' 
#' This function splits the features matrix into groups based on the group parameter. Each group is scaled to have 
#' 0 mean and 1 standard deviation. The matrix is then ungrouped and returned as a tibble.
#' 
#' @param features numeric matrix. Each column corresponds to a different feature and the row correspond to the points that
#' the features take on
#' @param group factor. grouping to scale by
#' @return A tibble with the passed in featuers scaled by group.
#' @export
FeatureNormalization <- function(features, group) {
	features <- as.data.frame(features)
	columns <- colnames(features)
	normalized <- features %>%
		dplyr::mutate(group=group) %>%
		dplyr::group_by(group) %>%
		dplyr::mutate_at(columns, scale) %>%
		dplyr::ungroup() %>%
		dplyr::select(-group)
	return(normalized)
}

#' Threshold features
#' 
#' \code{ThresholdFeature} returns a logical vector the same length as (feature)
#' 
#' This function simply checks values in the (feature) vector to see if they are greater
#' than or equal to a (threshold) value.
#' 
#' @param feature numeric vector. 
#' @param threshold numeric. Threshold to check (feature) vector against.
#' @return A logical vector the same length as (feature) 
ThresholdFeature <- function(feature, threshold) {
	return(feature >= threshold)
}

#' Perform density-based clustering
#' 
#' \code{DensityClustering} Determine which thresholded cells belong to a dense clustering.
#' 
#' Performs DBSCAN clustering on the passed in coordinate values. Core points and reachable points
#' are considered the same. Small clusters (fewer than 25 cells) are set as noise.
#' 
#' @param X numeric vector. first dimension coordinate
#' @param Y numeric vector. second dimension coordinate
#' @param positive logical vector. Same length as (X) and (Y). Element is TRUE if that cell meets previous thresholding step.
#' @param eps numeric. Epsilon argument for DBSCAN algorithm.
#' @param k.values integer. minPts argument for DBSCAN algorithm.
#' @return data frame with columns X, Y, and dense cluster number. 0 is assigned for cells not belonging to a dense cluster.
#' @export
DensityClustering <- function(X, Y, positive, eps, k.values=7) {
	
	# initialize coordinates
	coordinates <- cbind.data.frame(X, Y)
	coordinates <- dplyr::mutate(coordinates, cluster=0)
	
	# if no cells meet the threshold return the initialized zero values
	if (sum(positive) == 0) {
		return(coordinates)
	}
	
	# if epsilon could not be calculated return initialized zero values
	if (is.na(eps)) {
		return(coordinates)
	}
	
	# do density-based clustering
	distance.matrix <- stats::dist(coordinates[positive,])
	set.seed(1)
	res <- dbscan::dbscan(x=distance.matrix, eps=eps, minPts=k.values)
	
	# DBSCAN found no dense clusters
	if (all(res$cluster == 0)) {
		return(coordinates)
	}
	
	# we want all DBSCAN clusters to have at least k members
	# ties break this so set all of these clusters to noise (0)
	cluster.counts <- plyr::count(res$cluster)
	small.clusters <- which(cluster.counts$freq < 25)
	small.cluster.numbers <- cluster.counts$x[small.clusters]
	res$cluster[res$cluster %in% small.cluster.numbers] <- 0
	
	# shift cluster numbers for lost small clusters
	cluster.numbers <- sort(unique(res$cluster))
	new.clusters <- seq_along(cluster.numbers) - 1
	
	# update the cluster with the DBSCAN cluster number
	coordinates$cluster[positive] <- new.clusters[match(res$cluster, cluster.numbers)]
	return(coordinates)
}


#' Calculate the circumradius of a triangle
#' 
#' \code{Circumradius} Calculate the circumradius of a pair of 3-tuples
#' 
#' @param points numeric vector. Length 6. First 3 elements correspond to X coordinates of the 3-tuple. Second 3 to Y coordinates.
#' @return numeric. The circumradius of the triangle.
Circumradius <- function(points) {
	x1 <- points[1]
	x2 <- points[2]
	x3 <- points[3]
	y1 <- points[4]
	y2 <- points[5]
	y3 <- points[6]
	
	a <- det(matrix(c(x1,x2,x3,
					  y1,y2,y3,
					  1, 1, 1), nrow=3))
	bx <- -det(matrix(c(x1^2 + y1^2, x2^2 + y2^2, x3^2 + y3^2,
						y1, y2, y3,
						1, 1, 1), nrow=3))
	by <- det(matrix(c(x1^2 + y1^2, x2^2 + y2^2, x3^2 + y3^2,
					   x1, x2, x3,
					   1, 1, 1), nrow=3))
	c <- -det(matrix(c(x1^2 + y1^2, x2^2 + y2^2, x3^2 + y3^2,
					   x1, x2, x3,
					   y1, y2, y3), nrow=3))
	
	rad <- sqrt(bx^2 + by^2 - (4*a*c))/(2*abs(a))
	return(rad)
}


#' Perform alpha shapes segmentation
#' 
#' \code{AlphaShapes} Get a SpatialPolygon of the alpha shape of dense regions.
#' 
#' Remove triangles from the set of triangles given by the Delaunay triangulation which have circumradii larger than (alpha).
#' 
#' @param X numeric vector. first dimension coordinates of cells for which an alpha shape should be drawn.
#' @param Y numeric vector. second dimension coordinates of cells for which an alpha shape should be drawn.
#' @param alpha numeric. Threshold for circumradius of Delaunay triangulation triangles. Recommended to use (epsilon) from \code{GetEpsilon}.
#' @return SpatialPolygon. Aggregate of polygons identified by alpha shape borders.
#' @export
AlphaShapes <- function(X, Y, alpha) {
	# initialize coordinates
	coordinates <- cbind(X, Y)
	
	# delaunay triangulation
	tri <- geometry::delaunayn(coordinates)
	
	# which delaunay simplices have circumraidus smaller than alpha
	radii <- apply(tri, 1, function(x) Circumradius(coordinates[x,]))
	alpha.index <- which(radii <= alpha)
	
	# get alpha complex and create polygons
	alpha.complex <- tri[alpha.index,]
	alpha.complex <- cbind(alpha.complex, alpha.complex[,1])
	polygons <- apply(alpha.complex, 1, function(x) sp::Polygon(coordinates[x,]))
	polygons <- lapply(seq_along(polygons), function(i) sp::Polygons(list(polygons[[i]]), i))
	alpha.shape <- raster::aggregate(sp::SpatialPolygons(polygons))
	
	return(alpha.shape)
}


#' Determine if a cell falls within a SpatialPolygon.
#' 
#' \code{CellsinRegion} Check which cells in a sample fall within the area of a SpatialPolygon aggregate for that sample.
#' 
#' @param X numeric vector. first dimension coordinate
#' @param Y numeric vector. second dimension coordinate
#' @param shapes SpatialPolygon. aggregate of SpatialPolygons to check.
#' @param SID integer. Sample id that is being checked.
#' @return logical vector. Same length as (X) and (Y). TRUE if cell falls within SpatialPolygon. FALSE otherwise.
#' @export
CellsInRegion <- function(X, Y, shapes, SID) {
	# extract alpha shape and check if it exists
	alpha.shape <- shapes %>%
		dplyr::filter(sample_id == SID) %>%
		dplyr::pull(shape)
	if (length(alpha.shape) == 0) {
		return(rep(FALSE, length(X)))
	}
	# logical vector if a coordinate is the in the shape or not
	prevR::point.in.SpatialPolygons(X, Y, alpha.shape[[1]])
}


#' Get the rough area of the entire sample.
#' 
#' \code{SampleArea} Use the convex hull to get an estimate of the area of the entire sample.
#' 
#' @param X numeric vector. first dimension coordinate
#' @param Y numeric vector. second dimension coordinate
#' @return numeric. Area of the sample.
#' @export
SampleArea <- function(X, Y) {
	ch <- grDevices::chull(X, Y)
	ch <- c(ch, ch[1])
	coord <- cbind(X,Y)
	ch_poly <- sp::Polygon(coord[ch,], hole=FALSE)
	return(ch_poly@area)
}


#' Extract morphometrics of shapes.
#' 
#' \code{GetShapeParams} For a single SpatialPolygon extract its area and its perimeter.
#' 
#' @param shape SpatialPolygon. A single shape
#' @return data frame with a column for area and a column for perimeter.
#' @export
GetShapeParams <- function(shape) {
	shape <- shape[[1]]
	area <- rgeos::gArea(shape, byid=TRUE)
	perimeter <- rgeos::gLength(shape, byid=TRUE)
	return(data.frame(area=area, perimeter=perimeter))
}

DensityBasedSegmentation <- function(features, coordinates, marker, cell_type=NULL, threshold=NULL, epsilon=NULL, PLOT=FALSE) {
	
	# set default quantile threshold on expression to 0.75
	if (is.null(threshold)) {
		threshold <- 0.75
	}
	
	# normalize expression
	if (is.null(cell_type)) {	# normalization when all cell types are the same
		expr <- features %>%
			select(sample_id, cell_id, one_of(marker))
		norm_feat <- expr %>%
			select(one_of(marker)) %>%
			mutate_all(function(x) x >= stats::quantile(x, probs=threshold)) %>%
			mutate(positive = rowSums(.) > 0)
		expr <- expr %>%
			mutate(positive = norm_feat$positive)
	} else {					# normalization when cell types are given
		expr <- features %>%
			select(sample_id, cell_id, one_of(marker)) %>%
			left_join(cell_type, by=c("sample_id", "cell_id"))
		
		norm_feat <- FeatureNormalization(select(expr, one_of(marker)), pull(expr, cell_type))
		norm_feat <- norm_feat %>%
			mutate_all(function(x) x >= stats::quantile(x, probs=threshold)) %>%
			mutate(positive = rowSums(.) > 0)
		expr <- expr %>%
			mutate(positive = norm_feat$positive)
	}
	
	if (is.null(epsilon)) {
		epsilon <- coordinates %>%
			group_by(sample_id) %>%
			summarise(epsilon = GetEpsilon(X, Y)) 
	}
	
	
	## density-based clustering
	clusters <- coordinates %>%
		bind_cols(select(expr, positive)) %>%
		left_join(epsilon, by=c("sample_id")) %>%
		group_by(sample_id) %>%
		do(DensityClustering(.$X, .$Y, .$positive, unique(.$epsilon))) 
	
	
	## alpha shapes
	shapes <- clusters %>%
		left_join(epsilon, by=c("sample_id")) %>%
		filter(cluster > 0) %>%
		group_by(sample_id, cluster) %>%
		do(shape=AlphaShapes(.$X, .$Y, unique(.$epsilon))) %>%
		ungroup() %>%
		group_by(sample_id) %>%
		do(shape=if (nrow(.) > 1)
			do.call(raster::bind, .$shape)
			else .$shape[[1]]) %>%
		ungroup()
	
	## find cells that are in each region
	region <- coordinates %>%
		group_by(sample_id) %>%
		mutate(region=CellsInRegion(X, Y, shapes, unique(sample_id))) %>%
		ungroup()
	
	
	## get area and perimeter of each region
	shape_params <- shapes %>%
		group_by(sample_id) %>%
		do(GetShapeParams(.$shape)) %>%
		ungroup()
	
	## get sample area for each sample and combined features
	shape_params <- coordinates %>%
		group_by(sample_id) %>%
		summarise(sample_area=SampleArea(X,Y)) %>%
		left_join(shape_params, by=c("sample_id")) %>%
		mutate(proportion=area/sample_area,
			   circularity=log(area/perimeter)) %>%
		left_join(shapes, by=c("sample_id"))
	
	##
	sample_morphometrics <- shape_params %>%
		group_by(sample_id) %>%
		summarise(area=sum(area),
				  perimeter=sum(perimeter),
				  proportion=sum(proportion),
				  circularity=median(circularity),
				  n_shapes = ifelse(is.na(area), 0, n()),
				  .groups="drop") %>%
		tidyr::replace_na(list(area=0, perimeter=0, proportion=0, circularity=0))

	segmentation <- coordinates %>%
		select(sample_id, cell_id) %>%
		mutate(positive=expr$positive,
			   cluster=clusters$cluster,
			   region=region$region)
	
	return(list(segmentation=segmentation, sample_morphometrics=sample_morphometrics, shape_params=shape_params))
}

