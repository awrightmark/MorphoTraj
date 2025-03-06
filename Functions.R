# Calculates pairwise distances between all points in multi-dimensional space using a matrix (dimensions as columns, points as rows) and a factor labeling groups as inputs; only calculates distances between points not part of the same group
pairwiseDistances <- function(matrix, group){
  dist.mtrx <- as.matrix(dist(matrix, method = "euclidean"))
  
  inter_group_distances <- matrix(NA, nrow = length(group), ncol = length(group))
  rownames(inter_group_distances) <- colnames(inter_group_distances) <- rownames(matrix)
  
  for (i in 1:length(group)) {
    for (j in 1:length(group)) {
      if (group[i] != group[j]) {
        inter_group_distances[i, j] <- dist.mtrx[i, j]
      }
    }
  }
  
  heatmap(inter_group_distances[group != group[length(group)],group != group[1]], Rowv = NA, Colv = NA, scale = "none",
          col = rev(heat.colors(256)), main = "Euclidean Distance Heatmap",
          na.rm = TRUE)
  
  inter_group_distances[lower.tri(inter_group_distances, diag = TRUE)] <- NA  # Remove duplicates & self-distances
  dist_flat <- as.data.frame(as.table(inter_group_distances))  # Convert to long format
  dist_flat <- na.omit(dist_flat)  # Remove NA values
  # dist_flat$Distance <- round(dist_flat$Distance, 2) # Round extra decimal places
  colnames(dist_flat) <- c("Point1", "Point2", "Distance")
  rownames(dist_flat) <- 1:nrow(dist_flat)
  
  rank_order_dist <- dist_flat[order(dist_flat$Distance), ]
  rownames(rank_order_dist) <- 1:nrow(rank_order_dist)
  return(rank_order_dist)
}

# Function for drawing a morphogenetic trajectory based on sum of squares minimized distances to the line per stage and corresponding arrow position
draw.traj <- function(labeledMatrix, species, line.col){
  # Under the Hood
  pc.local <- prcomp(labeledMatrix[rownames(labeledMatrix) == species,])
  pc.local.vector <- pc.local$rotation[whichPCs,1]
  pc.local.slope <- pc.local$rotation[whichPCs[2],1]/pc.local$rotation[whichPCs[1],1]
  pc.local.intercept <- pc.local$center[whichPCs[2]] - pc.local.slope*pc.local$center[whichPCs[1]]
  series.minx <- min(labeledMatrix[rownames(labeledMatrix) == species,whichPCs[1]])
  series.maxx <- max(labeledMatrix[rownames(labeledMatrix) == species,whichPCs[1]])
  series.miny <- min(labeledMatrix[rownames(labeledMatrix) == species,whichPCs[2]])
  series.maxy <- max(labeledMatrix[rownames(labeledMatrix) == species,whichPCs[2]])
  y_at_xmin <- pc.local.slope * series.minx + pc.local.intercept
  y_at_xmax <- pc.local.slope * series.maxx + pc.local.intercept
  x_at_ymin <- (series.miny - pc.local.intercept) / pc.local.slope
  x_at_ymax <- (series.maxy - pc.local.intercept) / pc.local.slope
  valid_points <- data.frame(x = c(series.minx, series.maxx, x_at_ymin, x_at_ymax), y = c(y_at_xmin, y_at_xmax, series.miny, series.maxy))
  valid_points <- valid_points[valid_points$x >= series.minx & valid_points$x <= series.maxx & valid_points$y >= series.miny & valid_points$y <= series.maxy, ]
  
  # Draw arrow in direction that minimizes sum of squares distance between each point in the series and that respect position along the arrow
  numPts <- nrow(labeledMatrix[rownames(labeledMatrix) == species,])
  seq_x <- seq(valid_points$x[1], valid_points$x[2], length.out = numPts) # change 6 here to not be specific to this case
  seq_y <- pc.local.slope * seq_x + pc.local.intercept
  forward_dist_sq <- sum((labeledMatrix[rownames(labeledMatrix) == species,whichPCs[1]] - seq_x)^2 + (labeledMatrix[rownames(labeledMatrix) == species,whichPCs[2]] - seq_y)^2)
  reverse_dist_sq <- sum((labeledMatrix[rownames(labeledMatrix) == species,whichPCs[1]] - rev(seq_x))^2 + (labeledMatrix[rownames(labeledMatrix) == species,whichPCs[2]] - rev(seq_y))^2)
  best_direction <- ifelse(forward_dist_sq < reverse_dist_sq, "Forward", "Reverse")
  Mean_SSq_perPt <- min(forward_dist_sq, reverse_dist_sq)/nrow(labeledMatrix[rownames(labeledMatrix) == species,])
  
  # Output or Product
  if (best_direction == "Forward") {
    arrows(seq_x[1], seq_y[1], seq_x[numPts], seq_y[numPts], col = line.col, lwd = 2, length = 0.15, angle = 25) # 6 here too
  } else {
    arrows(seq_x[numPts], seq_y[numPts], seq_x[1], seq_y[1], col = line.col, lwd = 2, length = 0.15, angle = 25)
  }
  return(Mean_SSq_perPt)
}
