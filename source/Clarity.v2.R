library("gplots")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("grid")

require(mclust)
  
###############################################################################################
#                                                                                             #
#  Commandline Example to run purity estimate:                                                #
#  R CMD BATCH --no-save --no-restore '--args file.seg="conserting_file"  file.baf="baf_file" #
#    sample.name="SAMPLE_NAME"' Clarity.v2.R SJHGG001.out                                     #
#                                                                                             #
###############################################################################################

# read in the arguments listed at the command line
args=(commandArgs(TRUE))

# args is a list of character vectors
# First check to see if the argments are passed
# Then cycle thorugh each element of the list and evaluate the expressions

if(length(args)==0){
  print("No arguments supplied. Try again")
  quit(save = "no", status = 1, runLast = TRUE)
} else {
    for(i in 1:length(args)){
     eval(parse(text=args[[i]]))
    }
}

#base_dir <- "/home/dputnam/vcf2cna"
#working_directory <- "/home/dputnam/purity_benchmark/analyzed_files/wes/SJNBL011450_D1_G1"
#sample.name <- "SJNBL011450_D1_G1"
#file.name <- "SJNBL011450_D1_G1.high_20.out"
#analysis.type <- "fast" 

if ( exists("working_directory")) {
  setwd(working_directory)
} else {
  print("The working directory must be set, Program will abort.")
  q()
}
getwd()

file.baf=paste(sample.name, ".germ.baf", sep="")
file.seg=paste("Result/", file.name, "_CONSERTING_Mapability_100.txt", sep="")

# Create Analysis type to implement analysis on Full Conserting Runs ( "full")  vs VCF2CNA or EXOME ( "fast")

function_file = (paste(base_dir, "/source/Clarity.fun.R", sep=""))

source(function_file)
baf.window = 0


if (analysis.type == "fast")
{
  # Analysis type "fast" is for High20 files or VCF files
  baf.window = 300
} else if ( analysis.type == "full") {
  # Analysis type "full" is for counts derived from BAM files
  baf.window = 1000
} else {
  print("Incorrect analysis type provided.  Please use either fast or full")
  quit(save = "no", status = 1, runLast = TRUE)
}
  
print_id = FALSE
print.heatmap.plot = FALSE
print.density.plot = FALSE
seperation_threshold = 0.10
print.purity.plot = TRUE 
  
  
  ############################################################################
  #                                                                          #
  #  Read filenames into variables that will be operated on to compute       #
  #     tumor purity:                                                        #
  #        sample.baf         B-Allele Frequencey Data                       #
  #        sample.seg         CONCERTING Output                              #
  #        id.seg:            subset of sample.seg                           #
  #                                                                          #
  ############################################################################
  
  sample.baf         = read.table(file.baf,         header=T)
  sample.seg         = read.table(file.seg,         header=T)

  autosomes <- subset(sample.seg, chrom <= 22)
  weighted.gmean <- weighted.mean(autosomes$GMean, autosomes$num.mark)
  
  id.seg = 0  

  if( analysis.type == "fast")
  {
    id.seg = which(
    sample.seg$chrom < 23 & 
      (sample.seg$loc.end - sample.seg$loc.start > 999999)  & 
    
      # Diplod signal is 1.  We ensure that Germline is diploid by subtracting weighted gmean verifying the value is within 10% cutoff range
      abs(sample.seg$GMean - weighted.gmean) < 0.1 &
      
      # Ensures single copy gain, copy neutral, and single copy loss samples
      sample.seg$seg.mean < 0.60 &  
      sample.seg$seg.mean > -0.60)
  } else {  # analysis.type == "full"
    id.seg = which(
    sample.seg$chrom < 23 &
      (
        # length.ratio is a measure of the quality of the read
        # Get either regions with at least 1 MB, ignorming small segments or
        (sample.seg$length.ratio > 0.9 &
         sample.seg$loc.end - sample.seg$loc.start > 999999)
         |
        # Get large segments (whole chromosome regions)
        (sample.seg$length.ratio > 0.75 &
         sample.seg$loc.end - sample.seg$loc.start > 9999999)
      )  &
      # Diplod signal is 1.  We ensure that Germline is diploid by subtracting weighted gmean verifying the value is within 10% cutoff range
      abs(sample.seg$GMean - weighted.gmean) < 0.1 &

      # Ensures single copy gain, copy neutral, and single copy loss samples
      sample.seg$seg.mean < 0.60 &
      sample.seg$seg.mean > -0.60)
  } 

  if( length(id.seg) == 0)
  {
    print("There are no segments that meet selection criteria. Purity cannot be estimated")
    print( "final_purity: NA")
    var.name <- paste("purity_", sample.name, sep="")
    estimation <- paste(var.name, ".txt", sep="")
    tmp.loc <- "Result/"
    purity.file.loc <- paste(tmp.loc, estimation, sep="")
    fileConn<- file(purity.file.loc)
    write("NA", fileConn)
    close(fileConn)
    quit(save = "no", status = 0, runLast = TRUE) 
  }

  if( print_id[1] == "all")
  {
    print_id <- id.seg
  }
  
  ############################################################################
  #                                                                          #
  #  Operation on Data: Iteration is over id.seg                             #
  #                                                                          #
  ############################################################################

  print ( c('total segments:', length(id.seg)))

  # initialize vector to hold purity estimates from CNV
  purity.list <- vector('list',length(id.seg)) 

  pos <- 1
  for (i in id.seg)
  {
    print ( c('Analyzing segment', i, pos))

    # 1 copy number gain cutoff
    if ( sample.seg$seg.mean[i] / sample.seg$GMean[i] > 0.5)
    {
      cna.mean = 0.5
    }
      
    # 1 copy number loss cutoff
    else if( sample.seg$seg.mean[i] / sample.seg$GMean[i] < -0.5)
    {
      cna.mean = -0.5
    }
    # region between 1 copy number gain and 1 copy number loss
    else
    {
      cna.mean = sample.seg$seg.mean[i] / sample.seg$GMean[i]
    }
      
    # if the cna.mean value falls below a threshold, set the sample to be diplod ( cna.mean = 0)
    if( abs(cna.mean) < 0.05)
    {
      cna.mean = 0
    }

    # count the total number of BAFS on a Concerting segment
    baf.on.segment = which(
      (sample.seg$chrom[i] == sample.baf$Chromosome) &
      (sample.baf$Position <  sample.seg$loc.end[i]) &
      (sample.baf$Position >= sample.seg$loc.start[i]))

    baf.total <- length(baf.on.segment)


    # IF total number of BAFS are less than baf.window (1000), skip to next iteration
    if(baf.total <= baf.window)
    {
      purity.baf.values    <- list(data.frame(purity = "NA", cna.mean = "NA", distance = "NA", chromosome = "NA", stringsAsFactors=FALSE)) 
      purity.list[[pos]]   <- purity.baf.values
      pos = pos + 1
      next
    }

    # baf.segments is a vector of lists, where each list contains baf.window size of BAF Ids
    baf.segments <- split_segment(baf.on.segment, baf.window)
    baf.segments.length <- sapply(baf.segments, length)
   
    # This piece of code handles the case of 1 item in the list
    if( length(baf.segments) == 1)
    {
      purity.baf.values <- list(data.frame(purity = "NA", cna.mean = "NA", distance = "NA", chromosome = "NA",  stringsAsFactors=FALSE)) 
    } else {
      purity.baf.values <- list(0, length(baf.segments))
    }
    
    total_segments <- length(baf.segments)
    baf_heatmap_plots <- vector(mode = "list", length = total_segments) 
    baf_density_plots <- vector(mode = "list", length = total_segments) 

    # iterate over baf.segments to compute a purity value based on each segment
    for (j in 1:total_segments)
    {
      # initialize a 100 x 100 matrix to 0
      heatmap_matrix <- matrix(0, nrow=100, ncol=100)
      
      # iterate over bafs in a segment to create the matrix
      segment.length <- baf.segments.length[j]
      for (k in 1:segment.length)
      {
        id <- baf.segments[[j]][k]
        baf_g = sample.baf$BAF_G[id]
        baf_d = sample.baf$BAF_D[id]
        
        # multiply decimal by matrix resolution
        baf_g_shift <- baf_g * 100
        baf_d_shift <- baf_d * 100
        
        # Round the number to the correct bin
        baf_g_bin = floor(baf_g_shift + 0.5)
        baf_d_bin = floor(baf_d_shift + 0.5)
 
        # increment identified bin
        value <- heatmap_matrix[baf_g_bin, baf_d_bin]
        heatmap_matrix[baf_g_bin, baf_d_bin] <- value + 1
      }
      # sum the columns of the heatmap_matrix
      heatmap_vector <- apply(heatmap_matrix, 2, sum)
      number <- sum(heatmap_matrix)
      
      interval <- find.mode(heatmap_vector)
      
      left_point <- interval$lp 
      right_point <- interval$rp
      area_diff <- interval$auc.dif

      purity <- "NA"
      distance <- right_point - left_point

      # QA Checks on BAF seperation
      if ( distance < seperation_threshold)
      {
        left_point <- 0.5
        right_point <- 0.5
        distance <- 0
      }
      

      # if the midpoint of the heatmap vector is greater than the peak
      # on the right or left interval, do not include sample in the
      # purity estimate set.
      if ( interval$ml < heatmap_vector[50] || interval$mr < heatmap_vector[50])
      {
        left_point <- 0.5
        right_point <- 0.5
        distance <- 0
      }
     
      
      # Case of copy number gain without BAF separation
      # This indicates that copy number gain derived from balanced amplification
      if( distance == 0 & cna.mean > 0 & area_diff < 1.5)
      {
        purity <- abs(cna.mean)
      }
      else if( distance > 0)
      {
        # Standard purity calculation
       
        # I need to understand the change to from the commented version to the live version
        # purity <- distance + 2 * abs(cna.mean) - abs( cna.mean / (1 + cna.mean))
        
         purity <- 2 * abs(cna.mean) +( distance - abs( cna.mean / (1 + cna.mean))) * ( 1 + cna.mean)

        
        # If purity estimation falls below concerting prediction, use concerting prediction.  
        # Theoretically, purity is equal to 2 * abs(cna.mean)
        if( purity < 2 * abs(cna.mean))
        {
          purity <- 2 * abs(cna.mean)
        }
      }

      # Remove low balanced amplification estimates
      if( cna.mean <= 0.25 & distance == 0 & analysis.type =="full")
      {
        purity <- "NA"
      }

      # put a ceiling on the purity calculation
      if( purity != "NA" & purity > 1.0)
      {
        purity <- 1.0
      }

      chromosome <- sample.seg$chrom[[i]]

      relevant_values.df <- data.frame(purity, cna.mean, distance, chromosome,  stringsAsFactors=FALSE)
      purity.baf.values[[j]] <- relevant_values.df
      
      if( print.density.plot & i %in% print_id)
      {
        title <- create_title(sample.name, i, pos, j, total_segments,  number, cna.mean, distance, purity, area_diff)
        density_plot <- create_density_plot( title, heatmap_vector, left_point, right_point)
        baf_density_plots[[j]] <- density_plot
      }
      if( print.heatmap.plot & i %in% print_id)
      {
        title <- create_title(sample.name, i, pos, j, total_segments, number, cna.mean, distance, purity, area_diff)
        heatmap_plot <- create_heatmap_plot( title, heatmap_matrix, left_point, right_point)
        baf_heatmap_plots[[j]] <- heatmap_plot
      }
    } # end for (j in baf.segments)
    if( print.density.plot & !print.heatmap.plot & i %in% print_id)
    {
      print.plots( sample.name,  baf_density_plots, 6, i, "d")
    }
    if( print.heatmap.plot & !print.density.plot & i %in% print_id)
    {
      print.plots( sample.name, baf_heatmap_plots, 6, i, "h")
    }
    if( print.heatmap.plot & print.density.plot & i %in% print_id)
    {
      combined_plots <- as.vector(rbind(baf_heatmap_plots, baf_density_plots))
      print.plots( sample.name, combined_plots, 6, i, "m")
    }
    
    purity.list[[pos]] <- purity.baf.values
    pos = pos + 1
  } # end for (i in id.seg)

  # Combine all of the purity estimates into one list
  purity.estimates <- unlist(purity.list, recursive = FALSE)
  purity.answers   <- do.call("rbind", purity.estimates)
  
  # Extract all values with purity not set to NA 
  purity.set       <- subset(purity.answers, purity != "NA")
  purity.set       <- subset(purity.set, purity >= 0.1)
  
  # Add new Type column to data frame based on cna.mean value
  purity.set       <- within(purity.set, {Type = ifelse(cna.mean >=0, "Gain", "Loss")})
  purity.set[purity.set$cna.mean == 0 & purity.set$purity > 0, "Type"] <- "CnLoh"
  purity.set[purity.set$cna.mean > 0  & purity.set$distance == 0,"Type"] <- "BalAmp" 

  purity.set <- subset(purity.set, Type != "BalAmp")

  # Must have a minimum of 20 points to estimate purity
  if ( dim(purity.set)[1] < 20)
  {
    print("The algorithm cannot predict tumor purity for this sample")
    num_points <- dim(purity.set)[1]
    print(c("The number of data points is: ", num_points))
    print( "final_purity: NA")
    var.name <- paste("purity_", sample.name, sep="")
    estimation <- paste(var.name, ".txt", sep="")
    tmp.loc <- "Result/"
    purity.file.loc <- paste(tmp.loc, estimation, sep="")
    fileConn<- file(purity.file.loc)
    write("NA", fileConn)
    close(fileConn)
    quit(save = "no", status = 0, runLast = TRUE)
  }
  
  
  # For Graphing add position variable to the data frame that is the length of the subset
  purity.pos       <- c(1:dim(purity.set)[1])
  purity.set$position <- purity.pos

  if( print.purity.plot)
  {
    data_plot <- create_purity_plot( purity.set)
    tmp.loc <- "Result/"
    print_purity_plot( data_plot, sample.name, tmp.loc)
  }
  
  # Retrieve the Purity values to send to MClust
  purity.value     <- as.numeric(purity.set$purity)
  
   # estimate gaussion mixture value
   # give bic for each G 1, 2, ..., 20
   # select the best model that minimizes bic
   # return the model based on the the G number
   # Maximum of 20 clusters
   # minimize the bic information criteria
   # return the model
  
  max_groupings <- max.clusters( length(purity.value), 5, 200)

  # This code handles the case where all the purity estimates are the same.
  # This hardcodes max_groupings to 1
  test_scalar = purity.value[1]
  result_test = lapply(purity.value, function(x) x==test_scalar)
  result.answers <- unlist(result_test, recursive = FALSE)
  if ( any(result.answers == FALSE)) {
    max_groupings <- max_groupings
  } else {
    max_groupings <- 1
  }


  if( max_groupings > 1) {
    purity.model <- try(Mclust(purity.value, G=1:max_groupings, modelNames="E"))
    while( class(purity.model) == "try-error"){
      max_groupings <- max_groupings - 1
      purity.model <- try(Mclust(purity.value, G=1:max_groupings, modelNames="E"))
    }

    print(c('max_groupings:', max_groupings))
    result <- rbind(purity.model$parameters$pro, purity.model$parameters$mean)
    print(c('result:' , result))
    print(c('result[2,]:', result[2,]))

    # purity is the value of the maximum cluster
    purity <- max(result[2,])
  }  else {
    purity <- test_scalar
  }

  print(c('purity:', purity))

  var.name <- paste("purity_", sample.name, sep="")
  estimation <- paste(var.name, ".txt", sep="")
  tmp.loc <- "Result/"
  purity.file.loc <- paste(tmp.loc, estimation, sep="")
  fileConn<- file(purity.file.loc)
  write(purity, fileConn)
  close(fileConn)
