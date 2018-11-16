multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


print.plots <- function( name, plot, plots.on.page, seg.id, type)
{
   total.pages <- ceiling(length(plot) / plots.on.page)

   first <- 1
   last <- plots.on.page
   for (i in 1:total.pages)
   {
     prefilename <- paste(name, ".",   sep="")
     filename <- paste(prefilename, type, sep="")
     file1    <- paste(filename, ".", sep="")
     file2    <- paste(file1, seg.id, sep="")
     file3    <- paste(file2, "."   , sep="")
     file4    <- paste(file3, i     , sep="")
     file5    <- paste(file4, ".png", sep="")
     
     allplots <- plot[first:last]

     png(file5, width = 1200, height = 800)
     multiplot( plotlist=allplots, cols=3)
     dev.off()
     
     first <- last + 1
     last <- last + plots.on.page
   }
}

create_title <- function( name, i, pos, j, total_segments, number, cna.mean, distance, purity, area_diff) {

  n_distance <- suppressWarnings(as.numeric(distance))
  n_purity <- suppressWarnings(as.numeric(purity))
  n_area <- suppressWarnings(as.numeric(area_diff))

  r_distance <- round(n_distance,  digits=2)
  r_purity   <- round(n_purity,    digits=2)
  r_area     <- round(n_area, digits=2)

  titleA <- paste("file=", name, sep="")
  titleB <- paste(titleA, "\nseg.id=", sep=" ")
  title1 <- paste(titleB, i, sep="")
  title2 <- paste(title1, "Pos=", sep=" ")
  title3 <- paste(title2, pos, sep="")
  title4 <- paste(title3, "N=", sep= " ")
  title5 <- paste(title4, number, sep= "")
  title6 <- paste(title5, "Conserting=", sep= " ")
  title7 <- paste(title6, cna.mean, sep= "")
  title8 <- paste(title7, "\nblock", sep= " ")
  title9 <- paste(title8, j, sep= " ")
  title10 <- paste(title9, "of", sep= " ")
  title11 <- paste(title10, total_segments, sep= " ")
  title12 <- paste(title11, "Distance=", sep= " ")
  title13 <- paste(title12, r_distance, sep= "")
  title14 <- paste(title13, "Purity=", sep= " ")
  title15 <- paste(title14, r_purity, sep= "")
  title16 <- paste(title15, "Area_diff=", sep= " ")
  title17 <- paste(title16, r_area, sep= "")
  return(title17)
}

create_density_plot <- function( title, heatmap_vector, left_point, right_point)
{
  dp <- data.frame(heatmap_vector)
  rnames <- 1:100
  snames <- rnames / 100
  dp$row <- snames

  max_height <- max(heatmap_vector)

  p <- ggplot(dp, aes(x=row, y=heatmap_vector)) + 
       geom_point(shape=1) + 
       geom_smooth( method=loess, span = 0.2) + 
       geom_vline(colour = "black", xintercept = left_point) +
       annotate("text", x=left_point - 0.05, y=max_height, label=left_point) +  
       geom_vline(colour = "black", xintercept = right_point) + 
       annotate("text", x=right_point + 0.05, y=max_height, label=right_point) +  
       theme_bw() +
       scale_x_continuous(name="Tumor B-Allele Frequency") + 
       ylab("Count") +
       ggtitle(title)
  
  return(p)
}

create_heatmap_plot <- function( title, heatmap_matrix, left_point, right_point)
{
  hm <- data.frame(heatmap_matrix)
  rnames <- 1:100
  snames <- rnames / 100
  colnames(hm) <- snames
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  df <- cbind(snames, hm)
  df.m <- melt(df, id.vars = "snames")
  p <-
     ggplot( df.m, aes(x = variable, y = snames, fill = value)) +
     geom_tile() +
     geom_vline(colour = "white", xintercept = left_point * 100) +
     annotate("text", x=left_point * 100, y=1.03, label=left_point) +  
     geom_vline(colour = "white", xintercept = right_point * 100) + 
     annotate("text", x=right_point * 100, y=1.03, label=right_point) +  
     scale_fill_gradientn(colours = myPalette(100)) +
     theme_bw() +
     scale_x_discrete(
        name="Tumor BAF",
        breaks=seq(0, 1, 0.25)) +
     ylab("Germline BAF") +
     ggtitle(title)
  return(p)
}

create_purity_plot <- function( purity.set)
{
  # Create consistant color scale between graphs
  library(RColorBrewer)
  my.colors <- brewer.pal(4, "Set1")
  names(my.colors) <- c("Gain","Loss","CnLoh","BalAmp")
  col.scale <- scale_colour_manual(name = "Type", values = my.colors)

  # Create a plot of the final purity values 
  vis_plot <-
    ggplot(purity.set, aes(x = position, y = as.numeric(purity), colour = Type)) +
    geom_point(shape=16, size=3) +
    ylab("Tumor Purity") +
    theme(panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank(),
       axis.line = element_line(colour="black"),
       axis.text=element_text(size=14),
       axis.title=element_text(size=14),
       legend.text = element_text(size=12)) +
    scale_x_continuous(breaks=NULL) +
    scale_y_continuous(limits = c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1.0)) +
    ggtitle( sample.name) +
    theme(axis.title.x=element_blank())

    data_plot <- vis_plot + col.scale
  return(data_plot)
}

print_purity_plot <- function( purity_plot, sample_name, working_directory)
{
  # print final purity estimates
  var.name <- paste("purity_", sample_name, sep="")
  image.name <- paste(var.name, ".png", sep="")
  tmp.loc <- paste(working_directory, "/", sep="")
  image.loc <- paste(tmp.loc, image.name, sep="")
  png(image.loc, width = 600, height = 600)
  plot(purity_plot)
  dev.off()

  output.name <- paste(var.name, ".csv", sep="")
  output.loc  <- paste(tmp.loc, output.name, sep="")

  # Write Data behind purity estimate plot
  write.csv(purity.set, file=output.loc)
}


split_segment <- function(baf.on.segment, baf.window)
{
  # Compute the size of the last segment
  rem.size <- length(baf.on.segment) %% baf.window
  
  # Compute the total number of segments in the block os size baf.window
  segment.list.size <- ceiling(length(baf.on.segment) / baf.window)

  # If the last segment is smaller than half the size of baf window, combinine it 
  # With the previous segment by decreasing the segment list size by one
  if( rem.size < (baf.window / 2))
  {
    segment.list.size <- segment.list.size - 1 
  }
  
  # Pre allocate a vector to hold the ids of each segment
  baf.segments.list <- vector('list', segment.list.size)

  # Initialize the left and right boundaries
  left  <- 0
  right <- left + baf.window 

  # Handle the case of only one segment
  if ( segment.list.size == 1 )
  {
    right <- length(baf.on.segment)
  }

  # Iterate over the final number of segments
  for( i in 1:segment.list.size)
  {
    # Each Segment is given by the interval between the left and right boundary
    baf.segments.list[[i]] <- baf.on.segment[left:right]
    
    # Update the left boundary
    left  <- right + 1
    
    # On the iteration prior to the last iteration the right boundary is either
    #     1) The final end point: length(baf.on.segment)
    #     2) The window size: baf.window
    if( i == segment.list.size - 1)
    {
      right <- length(baf.on.segment)
    }
    else
    {
      right <- right + baf.window
    }
    
  }
  # Returns a vector of lists, where each list contains bafwindow size of BAF Ids
  return(baf.segments.list)
}

max.clusters <- function( total_count, max_start, step)
{
  max <- max_start
  iter_max <- floor(total_count / step)
  if( iter_max > 0)
  {
    for ( i in seq(from=1,to=total_count, by=step))
    {
      if(i > 1)
      {
        print( i)
        max = max + 1
      }
    }
  }
  
# Maximum of 20 clusters
  if (max > 20)
  {
    max <- 20
  }
  return(max)
}

find.mode <- function( heatmap_vector)
{
  # Fit loess line through points
  dp <- data.frame(heatmap_vector)
  rnames <- 1:100
  snames <- rnames / 100
  dp$row <- snames

  # fit a line through the points using loess
  fit = predict(loess(heatmap_vector~row, dp, span=0.2, degree=2), dp$row)
  
  # Put the fit line into a data frame
  fit.df <- data.frame(fit)
  fit.df$row <- snames
  
  # Convert any negative number to zero on the fit line
  fit.df$fit[fit.df$fit < 0] <- 0

  # Split the set into 2 halves
  left_half <- fit.df$fit[1:51]
  right_half <- fit.df$fit[49:100]
   
  # Find the maximum point on left side and right side
  max_left <- max(left_half)
  max_right <- max(right_half)

  # Find the index of the max points on the left side and right side
  left_index <- match( max_left, left_half)
  right_index <- match( max_right, right_half)
  right_index <- right_index + 48

  # Compute the area of the left and right halves
  lh.x <- seq(0.01, 0.51, length=51)
  rh.x <- seq(0.49, 1.00, length=52)

  left_auc  <- sum(diff(lh.x) * (head(left_half,  -1)+tail(left_half,  -1)))/2
  right_auc <- sum(diff(rh.x) * (head(right_half, -1)+tail(right_half, -1)))/2

  # Compute the difference between halves
  area_diff <- abs(left_auc - right_auc)

  lp <- left_index / 100
  ml <- max_left
  rp <- right_index / 100
  mr <- max_right
  auc.dif <- area_diff

  result = data.frame(lp, rp, auc.dif, ml, mr)
  return(result)
}





