#' Set genomic start position to accumulated start position.
#'
#' @param dataf the df must have 2 cols, "chromosome" character show as 1,2,...,"X","Y",
#' the "start" column must be numeric.
#'
#' @return A data with accumulated start position, called in bulkLOHplot func.
#' @export
#'
#' @examples \dontrun{
#' loh <- read.table(system.file("extdata/","p6LOH.txt",package = "stmut"), sep="\t",
#'  header = TRUE)
#'  loh <- loh[,c(1,2,3,10)]
#'  colnames(loh) <- c("chromosome", "start", "end","tumorshift")
#'  loh$chromosome <- substr(loh$chromosome,4,nchar(loh$chromosome))
#'  df8 <- accumStartPos(loh)
#' }
accumStartPos <- function(dataf){
  chromosome <- NULL
  chr <- unique(dataf$chromosome)[-1]
  new_ch <- dataf %>% filter(chromosome %in% "1")

  for (ch in chr){ # for loop to generate accumulated start position
    if (ch == "X"){
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "23"
      new_ch <- rbind(new_ch, df)
    } else if (ch == "Y") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "24"
      new_ch <- rbind(new_ch, df)
    } else if (ch == "MT") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "25"
      new_ch <- rbind(new_ch, df)
    } else {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- dataf %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      new_ch <- rbind(new_ch, df)
    }}
  new_ch$chromosome <- as.numeric(new_ch$chromosome) # this code is very important to plot the chromosomes in order.
  return(new_ch)
}


#' Get accumulated start position for cnr and cns file for arm-level CNV plotting.
#'
#' @param cnr the dataframe of cnr file from cnvkit output
#' @param seg the dataframe of cns file from cnvkit output
#'
#' @return A list of 2 dataframe, the new_ch for cnr file; and the newSeg for cns
#' @export
#'
#' @examples \dontrun{
#' cnr <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cnr",
#' package = "stmut"), sep="\t", header = TRUE)
#' cns <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cns",
#' package = "stmut"), sep="\t", header = TRUE)
#' df9 <- accStartCNR_CNS(cnr=cnr, seg=cns)
#' }
accStartCNR_CNS <- function(cnr,seg){
  chromosome <- NULL
  chr <- c(as.character(seq(22)), "X", "Y")[-1]
  new_ch <- cnr %>% filter(chromosome %in% "1")
  newSeg <- seg %>% filter(chromosome %in% "1")

  for (ch in chr){
    if (ch == "X"){
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- cnr %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "23"
      new_ch <- rbind(new_ch, df)
      # change cns table gene start position
      sg <- seg %>% filter(chromosome %in% ch)
      sg$start <- as.numeric(sg$start) + as.numeric(val)
      sg$end <- as.numeric(sg$end) + as.numeric(valEnd)
      sg$chromosome <- "23"
      newSeg <- rbind(newSeg, sg)
    } else if (ch == "Y") {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- cnr %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      df$chromosome <- "24"
      new_ch <- rbind(new_ch, df)
      # change cns table gene start position
      sg <- seg %>% filter(chromosome %in% ch)
      sg$start <- as.numeric(sg$start) + as.numeric(val)
      sg$end <- as.numeric(sg$end) + as.numeric(valEnd)
      sg$chromosome <- "24"
      newSeg <- rbind(newSeg, sg)
    } else {
      val <- max(new_ch$start)
      valEnd <- max(new_ch$end)
      df <- cnr %>% filter(chromosome %in% ch)
      df$start <- as.numeric(df$start) + as.numeric(val)
      df$end <- as.numeric(df$end) + as.numeric(valEnd)
      new_ch <- rbind(new_ch, df)

      # change cns table gene start position
      sg <- seg %>% filter(chromosome %in% ch)
      sg$start <- as.numeric(sg$start) + as.numeric(val)
      sg$end <- as.numeric(sg$end) + as.numeric(valEnd)
      newSeg <- rbind(newSeg, sg)
    }}

  new_ch$chromosome <- as.numeric(new_ch$chromosome)
  newSeg$chromosome <- as.numeric(newSeg$chromosome)

  returnedList <- list(new_ch, newSeg)
  return(returnedList)
}


#' Arm-level copy number variations of bulk DNA-seq data.
#'
#' @param centmere the simplified 3-column hg38 centromere data
#' @param dcnr cnr file which is cnvkit output
#' @param dcns segment file which is cnvkit output
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return An arm-level CNV plot, because the ref pool is male, so the plot also
#' has chrY data for female but the plot shows a deep deletion of chrY for female
#' @export
#'
#' @examples \dontrun{
#' centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
#' sep = "\t", header = FALSE)
#' cnr <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cnr",
#' package = "stmut"), sep="\t", header = TRUE)
#' cns <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cns",
#' package = "stmut"), sep="\t", header = TRUE)
#' bulkCNVs(centmere=centm, dcnr=cnr,dcns=cns)
#' }
bulkCNVs <- function(centmere, dcnr, dcns){
  chromosome <- NULL
  start <- NULL
  end <- NULL
  avg <- NULL
  colnames(centmere) <- c("chromosome", "start", "end")
  centmere$chromosome <- substr(centmere$chromosome,4,nchar(centmere$chromosome))
  centmere[which(centmere$chromosome=="X"),1] <- "23"
  centmere[which(centmere$chromosome=="Y"),1] <- "24"
  centmere$chromosome <- as.numeric(centmere$chromosome)

  # get cnr data
  dcnr$chromosome <- substr(dcnr$chromosome,4,nchar(dcnr$chromosome))
  dcns$chromosome <- substr(dcns$chromosome, 4, nchar(dcns$chromosome))

  r <- accStartCNR_CNS(dcnr,dcns)
  df <-r[[1]]
  df2 <-r[[2]]

  # add previous chromosome max start position to centromere start position
  nums <- centmere$chromosome[-1]
  for (i in nums){
    data <- df %>% filter(chromosome %in% (i-1))
    val1 <- max(data$start)
    val2 <- max(data$end)
    centmere[which(centmere$chromosome == i),2]  <- centmere[which(centmere$chromosome == i),2] + val1
    centmere[which(centmere$chromosome == i),3]  <- centmere[which(centmere$chromosome == i),3] + val2
  }

  # calculate mean column to centmere data
  centmere <- centmere[order(centmere$chromosome),]
  meanC <- rep(0,24)
  centmere <- cbind(centmere, meanC)
  for (i in 1:nrow(centmere)){
    centmere[i,4] <- mean(c(centmere[i,2], centmere[i,3]))
  }

  x <- c(centmere$start, centmere$end)
  cenMean <- mean(x)

  # set chromosome names as factor for later use as X-axis.
  if (length(unique(df$chromosome)) == 24){
    chrom <- factor(c(paste0("chr",seq(22)),"chrX","chrY"),levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
  }else{
    chrom <- factor(c(paste0("chr",seq(22)),"chrX"),levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
  }

  pos <- df %>%
    group_by(chromosome) %>%
    summarize(avg = round(mean(start))) %>% pull(avg)
  t1 <- df %>% group_by(chromosome) %>% summarize(avg=max(start)) %>% pull(avg)

  ggplot(df) +
    annotate("point", x=df$start, y=df$log2, size=0.00002, color="grey", pch = 1,alpha = 0.2) +
    annotate("segment", x= df2$start,xend = df2$end, y=df2$log2, yend = df2$log2, color = "orange",size=1)+ # add segment layer to the plot.
    annotate("segment", x=0, xend=max(df$start)*1.02, y=0, yend=0,color="black")+
    annotate("segment", x=t1, xend=t1, y=-1.5, yend = 1, color = "black")+
    geom_vline(data = centmere, mapping = aes(xintercept=end), linetype=10, color="black") + # centromere line
    scale_x_continuous(breaks = pos, labels=chrom, expand = c(0, 0)) +
    labs(x = "chromosome", y = "copy ratio (log2)") +
    ylim(c(-1.5, 1)) +
    theme_classic() + # use white color as background
    theme(axis.text.x = element_text(angle = 90))
}


#' Allelic Imbalance or loss of heterozygocity plot.
#'
#' @param centmere hg38 centromere data
#' @param alle_imbal tumorshif/shift, which is the variant allele frequency deviated from 50% VAF
#'
#' @return An arm-level LOH plot
#' @export
#'
#' @examples \dontrun{
#' centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
#' sep = "\t", header = FALSE)
#' loh <- read.table(system.file("extdata/","p6LOH.txt",package = "stmut"), sep="\t",
#'  header = TRUE)
#'  loh <- loh[,c(1,2,3,10)]
#'  colnames(loh) <- c("chromosome", "start", "end","tumorshift")
#'  loh$chromosome <- substr(loh$chromosome,4,nchar(loh$chromosome))
#'  bulkLOHplot(centmere = centm, alle_imbal = loh)
#' }
bulkLOHplot <- function(centmere, alle_imbal){
  chromosome <- NULL
  avg <- NULL
  start <- NULL
  median <- NULL
  colnames(centmere) <- c("chromosome", "start", "end")
  centmere$chromosome <- substr(centmere$chromosome,4,nchar(centmere$chromosome))
  centmere[which(centmere$chromosome=="X"),1] <- "23"
  centmere[which(centmere$chromosome=="Y"),1] <- "24"
  centmere$chromosome <- as.numeric(centmere$chromosome)
  centmere <- centmere[order(centmere$chromosome),] # adjust centmere data to be in chromosome order
  centmere$chromosome <- as.character(centmere$chromosome)

  # retrieve centromere row position of each chromosome of my sample
  pos_centM <- c()
  arm_median <- c()

  colIdx <- grep("chromosome", colnames(alle_imbal))

  alle_imbal[which(alle_imbal$chromosome == "X"),colIdx] <- "23" #
  alle_imbal[which(alle_imbal$chromosome == "Y"),colIdx] <- "24" # either deleting or keeping this code depends on male or female

  # create chromosome arm tumorshift median matrix
  armMedian <- data.frame(chromosome=character(),
                          arm = character(),
                          median = double(),
                          start = double(),
                          end = double(), stringsAsFactors = FALSE)

  chrs <- unique(alle_imbal$chromosome)
  for (ch in chrs){
    d1 <- alle_imbal %>% filter(chromosome %in% ch)
    d1$start <- as.numeric(d1$start)
    d2 <- centmere %>% filter(chromosome %in% ch)
    for (pos in d1$start){
      if (pos >= d2$start){
        val <- which(d1$start == pos) # retrieve row position where the centromere should locate
        tumshft <- d1$tumorshift

        if (val == 1){
          med <- round(median(tumshft[1:length(tumshft)]),4)
          newR <- c(ch,"arm",med, min(d1$start),max(d1$end))
          armMedian[nrow(armMedian)+1,] <- newR
        } else {
          med1 <- round(median(tumshft[1:(val-1)]),4)
          newR1 <- c(ch, "arm1", med1,min(d1[c(1:(val-1)),]$start),max(d1[c(1:(val-1)),]$end))
          armMedian[nrow(armMedian)+1,] <- newR1
          med2 <- round(median(tumshft[val:length(tumshft)]),4)
          newR2 <- c(ch, "arm2", med2, min(d1[c(val:length(tumshft)),]$start),max(d1[c(val:length(tumshft)),]$end))
          armMedian[nrow(armMedian)+1,] <- newR2
        }
        break
      }
    }
  }

  # adjust start position by calling accumStartPos
  alle_imbal[which(alle_imbal$chromosome == "23"),colIdx] <- "X" #
  alle_imbal[which(alle_imbal$chromosome == "24"),colIdx] <- "Y"

  data <- accumStartPos(alle_imbal)
  df <- data

  # adjust arm start, end position fo armMedian data for plotting
  armMedian$chromosome <- as.numeric(armMedian$chromosome)
  armMedian$median <- as.numeric(armMedian$median)
  armMedian$start <- as.numeric(armMedian$start)
  armMedian$end <- as.numeric(armMedian$end)

  new_armD <- armMedian %>% filter(chromosome %in% 1)
  new_CM <- centmere %>% filter(chromosome %in% 1)

  nums <- unique(armMedian$chromosome)[-1]
  for (i in nums){
    darM <- armMedian %>% filter(chromosome %in% i)
    data <- df %>% filter(chromosome %in% (i-1))
    val1 <- max(data$start)
    val2 <- max(data$end)
    darM$start <- darM$start + val1
    darM$end <- darM$end + val2
    new_armD <- rbind(new_armD, darM)

    # adjust centromere start, end position for plotting
    CMdata <- centmere %>% filter(chromosome %in% i)
    CMdata$start <- CMdata$start + val1
    CMdata$end <- CMdata$end +val2
    new_CM <- rbind(new_CM, CMdata)
  }

  # preparing for plot
  if (length(unique(df$chromosome)) == 24){
    chrom <- factor(c(paste0("chr",seq(22)),"chrX","chrY"),levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
  }else{
    chrom <- factor(c(paste0("chr",seq(22)),"chrX"),levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"))
  }

  pos <- df %>%
    group_by(chromosome) %>%
    summarize(avg = round(mean(start))) %>% pull(avg)
  t1 <- df %>% group_by(chromosome) %>% summarize(max=max(start)) %>% pull(max)

  ggplot(df) +
    annotate("point", x=df$start, y=df$tumorshift, size=0.00002, color="grey", pch = 0.2,alpha = 0.5) + # each row in data
    annotate("segment", x=0, xend=max(df$start)*1.02, y=0, yend=0,color="black") + # x-axis
    annotate("segment", x=t1, xend=t1, y=0, yend = 0.25, color = "black") +
    annotate("segment", x=new_armD$start, xend = new_armD$end, y=new_armD$median, yend = new_armD$median, color="orange")+ # add tumorshift median for each arm
    geom_vline(data = new_CM, mapping = aes(xintercept=start), linetype=10, color="black") + # centromere line
    scale_x_continuous(breaks = pos, labels=chrom, expand = c(0, 0)) +
    labs(x = "chromosome", y = "TumorShift") +
    ylim(c(0, 0.25)) + # if change ylim(), making sure also change the annotate() yend=?, otherwise, the chromosome boundaries will disappear
    theme_classic() + # use white color as background
    theme(axis.text.x = element_text(angle = 90))
}


#' Grouping spots from the same cluster to generate new spots to save spots with few reads.
#'
#' @param data A dataframe of spots to be grouped
#' @param Treads The min reads of the a new spot, by default Treads=5000
#' @param NumSpt NumSpt*NumSpt nearest spots will be included to determine which
#' spots to grouped together. By default, NumSpt=6
#'
#' @importFrom raster pointDistance
#'
#' @return A dataframe of new spots
#' @export
#'
#' @examples \dontrun{
#' clst1 <- read.csv(system.file("extdata/", "c1.csv", package = "stmut"),header = FALSE)
#' df10 <- groupSpots(clst1, Treads=5000, NumSpt=8)
#' }
groupSpots <- function(data,Treads=5000,NumSpt=8){
  array_row <- NULL
  barcode <- NULL
  # to store newly combined group
  spt2 = NumSpt*NumSpt
  d1 <- data # save to d1 for later use in merge(), the original data variable will be consumed.

  groups <- data.frame(barcode = character(),
                       groups = character(),
                       stringsAsFactors = FALSE)
  n=1 # initialize new group number
  if (sum(data$TotalRDs) < Treads){ # if total reads less than 5000, put all spots in one group
    g1 <- data.frame(barcode = data$barcode, groups = rep(paste0("group",n),dim(data)[1]))
    groups <- rbind(groups,g1)
    n = n + 1
  }else{ # if the total reads > given number, grouping spots

    for (j in 1:nrow(data)){ # loop each row to group spots
      if (dim(data)[1] == 0){
        break
      }else if(dim(data)[1] == 1){
        g1 <- data.frame(barcode = data$barcode, groups = paste0("group",n))
        groups <- rbind(groups,g1)
        break
      }else{
        newdata <- data.frame(barcode = character(), # empty dataframe to store tempor 36 spots for distance calculation
                              reads = integer(),
                              array_row = integer(),
                              array_col = integer(),
                              stringsAsFactors = FALSE)

        # first do QC to only keep 36 points for distance calculation by only keep 6 continuous x-axis spots, 6 continuous y-axis
        b = unique(data$array_row)
        for (i in 1:length(b)){ # loop each row/x-axis
          nam <- data %>% filter(array_row %in% b[i]) # dataframe of the same row.
          a <- dim(nam)[1] # number of spots- spots having the same x-axis but different y-axis

          # each x position, only keep 6 spots with different y positions
          if (a <= NumSpt){
            newdata <- rbind(newdata, nam)
          }else{ # only keep 6 continuous spots
            nam <- nam[-c((NumSpt+1):a),] # remove spots too far away
            newdata <- rbind(newdata, nam)
          }
        }
        ########################
        # subset the first 36 spots
        if (nrow(newdata)>= spt2 & sum(newdata$reads) >= Treads){
          newc4 <- newdata[c(1:spt2),]
        }else if(nrow(newdata) >= spt2 & sum(newdata$reads) >= Treads){
          print(paste0("The ", spt2, " spots has less than ", Treads, "reads; Suggest to assign a bigger number to variable NumSpt. "))
          break
        }else if(nrow(newdata) < spt2 & sum(data$TotalRDs) >= Treads){
          newc4 <- newdata
        }else{
          newc4 <- newdata
        }
        ############################
        dall <- c(0) # store distance of the 35 distance, because distance of the same spot is 0.
        p1 = newc4[1,c("array_row","array_col")]
        for (i in 2:nrow(newc4)){
          p2 = newc4[i,c("array_row","array_col")]
          a <- pointDistance(p1,p2,lonlat = FALSE)
          dall <- append(dall,a)
        }
        newc4 <- cbind(newc4, dall, row.names=NULL)
        newc5 <- newc4[order(newc4$dall),] # sort the df based on distance to the first spot
        # get fewest spots whose sum of reads is greater than 25000

        if (sum(newc5$TotalRDs) < Treads){
          g1 <- data.frame("barcode"=newc5$barcode, "groups"=rep(paste0("group",n), dim(newc5)[1])) # new group
          groups <- rbind(groups, g1)
          break
        }else{
          for (i in 1:nrow(newc5)){
            reads = newc5$TotalRDs
            readT <- sum(reads[1:i])

            if (readT >= Treads){
              k = i # to filter each group
              break
            }
          }
          rmCode <- newc5$barcode[1:k] # spot barcodes grouped as a new spot
          g1 <- data.frame("barcode"=rmCode, "groups"=rep(paste0("group",n),length(rmCode))) # new group
          groups <- rbind(groups, g1)

          data <- data %>% filter(! barcode %in% rmCode) # remove the grouped spots and start loop again.
          n = n +1
        }
      }
    }} # loop each row ends; else
  ##########
  final <- merge(d1,groups, by = "barcode")
  return(final)
}


#' Generating new_spots_feature_bc.csv after calling groupSpots function.
#'
#' @param df1 dataframe of spots to be grouped
#' @param Treads min reads of new spot, default is 5000
#' @param NumSpt number of spots to be considered, default is 8
#' @param data Original filtered_feature_bc.csv
#'
#' @importFrom stats na.omit
#'
#' @return A list of 2 dataframe. One is grouping spots with grouping-name; the other
#' is newspots(grouped_spots + non_group_spots) feature barcode dataframe
#' @export
#'
#' @examples \dontrun{
#' clst1 <- read.csv(system.file("extdata/", "c1.csv", package = "stmut"),header = FALSE)
#' c1BC <- read.csv(system.file("extdata/", "c1_feature_bc.csv", package = "stmut"),header = FALSE)
#' df11 <- newSptBC(df1=clst1, Treads=5000,NumSpt=8,data=c1BC)
#' }
newSptBC <- function(df1,Treads=5000,NumSpt=8,data){
  cluster <- NULL
  uniqueName <- c()
  fdf <- data.frame(matrix(ncol = 6, nrow = 0)) # the final fdf has same dimension as df1 but adding group-name within each cluster
  name <- c("barcode","array_row","array_col","cluster","TotalRDs","groups")
  colnames(fdf) <- name
  n <- length(unique(df1$cluster)) # number of clusters
  n <- unique(df1$cluster)
  for (i in n){ # loop clusters
    clst <- df1 %>% filter(cluster == i)
    clst <- clst[order(c(clst$array_row,clst$array_col)),]
    clst <- na.omit(clst)
    df <- groupSpots(data = clst,Treads = Treads, NumSpt = NumSpt) # call R function
    fdf <- rbind(fdf,df)
  }
  fdf1 <- fdf #

  # generate new new_feature_bc.csv
  fdf$barcode <- str_replace(fdf$barcode,"-1",".1")
  keep <- data[,!(names(data) %in% c(fdf$barcode))] # not included for spotcombining
  combSpts <- data[,c(fdf$barcode)] # spots for grouping

  k=0 # number of new spots
  for (cl in unique(fdf$cluster)){
    fdf2 <- fdf %>% filter(cluster==cl)
    len <- length(unique(fdf2$groups))
    k = k + len
    for(gp in unique(fdf2$groups)){
      d1 <- (fdf %>% filter(cluster == cl) %>% filter(groups == gp))

      if (dim(d1)[1] == 1){
        d1 <- d1$barcode
        df1 <- combSpts[,c(d1)]
        nam1 <- d1[1]
        newcol <- df1
      }else{
        d1 <- d1$barcode
        df1 <- combSpts[,c(d1)]
        nam1 <- d1[1]
        newcol <- rowSums(df1)
      }
      keep <- cbind(keep,newcol)
      name <- colnames(keep)
      colnames(keep) <- replace(name, name == "newcol",nam1 )
      uniqueName <- append(uniqueName, nam1)
    }
  }
  returnList <- list(fdf1, keep,uniqueName)
  return(returnList)
}




#' Calculate Chromosome Arm CNVs Weighted Median to Represent the CNVs of that arm.
#'
#' @param data The cnr dataframe of each spot.
#' @param centmere The centromere dataframe of hg38.
#'
#' @importFrom matrixStats weightedMedian
#'
#' @return A list of 2 dataframe, the weighted cnr file and a summary file counting
#' how many genes in each arm.
#' @export
#'
#' @examples \dontrun{
#' cnr <- read.table(system.file("extdata/", "spot1_rep1.cnr", package = "stmut"), header = TRUE)
#' centm <- read.csv(system.file("extdata/", "hg38_centromereSimple.bed", package = "stmut"),
#' sep = "\t", header = FALSE)
#' df5 <- wtArmMedianOne(data = cnr, centmere = centm)
#' }
wtArmMedianOne <- function(data,centmere){
  chromosome <- NULL
  # when calculate rolling median, rolling through each arm within of each chromosome.
  mt <- data %>% filter(chromosome %in% "MT")
  dataNew <- data %>% filter(!chromosome %in% "MT") # remove MT
  dataNew[which(dataNew$chromosome == "X"),1] <- "23" # replace X with "23"

  chr <- unique(dataNew[, which(colnames(data)=="chromosome")])
  colnames(centmere) <- c("chromosome", "start", "end")
  centmere$chromosome <- substr(centmere$chromosome,4,nchar(centmere$chromosome))
  centmere[which(centmere$chromosome=="X"),1] <- "23"
  centmere[which(centmere$chromosome=="Y"),1] <- "24"
  centmere$chromosome <- as.numeric(centmere$chromosome)
  centmere <- centmere[order(centmere$chromosome),] # adjust centmere data to be in chromosome order
  centmere$chromosome <- as.character(centmere$chromosome)

  qArm <- c() # store only qArms

  df <- data.frame(characters=character(),
                   Integers=integer(),
                   Integers=integer(),
                   Characters=character(),
                   Doubles=double(),
                   Doubles=double(),
                   Doubles=double(),
                   Integers=integer(),
                   Doubles=double(),
                   stringsAsFactors=FALSE)
  dfArms <- data.frame(matrix(nrow = length(chr), ncol = 5)) # store number genes of each arm for cdt sorting use
  colnames(dfArms) <- c("chromosome","p_Genes","q_Genes","pArmEnds","CM_row_pos")

  for (ch in chr){ # first loop over chromosome
    data1 <- dataNew %>% filter(chromosome %in% ch) # chr1
    ctm1 <- centmere %>% filter(chromosome %in% ch) # centromere of chr1
    nch <- length(data1[,which(colnames(data1)=="log2")])
    #to classify chromosomes into 2 groups, with and withour p-arms
    for (pos in data1$start){ # second loop over start pos of each chromosome to get break-points for centromere
      if (pos > ctm1$start){
        val <- which(data1$start == pos) # row pos of p arm ends,q arm starts
        if (val == 1){
          qArm <- append(qArm, ch)
        }
        dfArms[ch,"chromosome"] <- ch
        dfArms[ch,"p_Genes"] <- val-1
        dfArms[ch,"q_Genes"] <- nrow(data1)-val+1
        dfArms[ch,"pArmEnds"] <- pos
        dfArms[ch,"CM_row_pos"] <- val-1
        break
      }
    } # end of second for loop

    if (ch %in% qArm){ # weighted-median for chromosomes without p-arms
      dat <- data1[c(1:nch), which(colnames(data1)=="log2")]
      w <- data1[c(1:nch), which(colnames(data1)=="weight")]
      m <- weightedMedian(dat,w)
      data1[which(data1$chromosome == ch),which(colnames(data1)=="log2")] <- m
      df <- rbind(df, data1)

    } else { # chromosome with both arms
      dat1 <- data1[c(1:(val-1)), which(colnames(data1)=="log2")]
      w1 <- data1[c(1:(val-1)), which(colnames(data1)=="weight")]
      m1 <- weightedMedian(dat1,w1)
      data1[c(1:(val-1)),which(colnames(data1)=="log2")] <- m1
      dat2 <- data1[c(val:nch), which(colnames(data1)=="log2")]
      w2 <- data1[c(val:nch), which(colnames(data1)=="weight")]
      m2 <- weightedMedian(dat2,w2)
      data1[c(val:nch),which(colnames(data1)=="log2")] <- m2
      df <- rbind(df, data1)
    }}

  returnList <- list(df,dfArms)
  return(returnList)
} # function ending




#' Reformat Gene Summary of Each Arm for cdt_filt_sort function use in the next step.
#'
#' @param cdt the filtered cdt file.
#' @param data A summary file generated by wtArmMedianOne function.
#'
#' @return A dataframe prepared for sorting the cdt based on bulk CNVs information.
#' @export
#'
#' @examples \dontrun{
#' cdt <- read.table(system.file("extdata/", "cdt.cdt", package = "stmut"), header = TRUE)
#' data3 <- read.table(system.file("extdata/", "summary.txt", package = "stmut"),
#' sep = "\t", header = TRUE)
#' df6 <- CtArmGenes(cdt = cdt, data = data3)
#' }
CtArmGenes <- function(cdt, data){
  dfgene <- data.frame(matrix(nrow = 0,ncol = 3)) # store chromosome,number of genes in each arm, arm_start_pos...
  colnames(dfgene) <- c("chromosome","arm_genes","gene_row")

  genes <- cdt[,3]
  spt <- names(cdt)[3]
  col1 <- cdt[,spt]

  i = 1
  for (arm in unique(genes)){ # because each arm has the same median, that is why unique(genes)= number of arms
    #print(arm)
    d0 <- cdt %>% filter(!!as.symbol(spt) == arm)
    ch <- str_split(d0[1,2],":",simplify = TRUE)[1]
    numGene <- dim(d0)[1]
    rowN <- which(col1 == arm)[1]
    dfgene[i,"chromosome"] <- ch
    dfgene[i,"arm_genes"] <- numGene
    dfgene[i,"gene_row"] <- rowN
    i = i +1
  }

  # extract p-arm, q-arm information
  d3 <- data.frame(matrix(nrow = dim(data)[1]*2, ncol = 2))
  colnames(d3) <- c("arms","genes")
  for (i in 1:dim(data)[1]){
    d3[2*i-1,1] = paste0(i,"p")
    d3[2*i-1,2] <- data[i,2]
    d3[2*i,1] <- paste0(i,"q")
    d3[2*i,2] <- data[i,3]
  }

  d4 <- d3 %>% filter(genes!=0)
  final <- cbind(d4, dfgene)
  final <- final[,-4]
  return(final)
}




#' Rank the spots by their similarity to the DNA-seq copy number alterations.
#'
#' @param cdt cdt file to be ranked
#' @param genes Number of genes in each arm
#' @param gainLoss Gain or Loss info of each arm, 1 represent gain, -1 represent loss.
#' @param rs Representative row position of each arm.
#'
#' @return A spots dataframe sorted by CNVs
#' @export
#'
#' @examples \dontrun{
#' cdt <- read.table(system.file("extdata/", "cdt.cdt", package = "stmut"), header = TRUE)
#' data4 <- read.csv(system.file("extdata/", "CtArmGenSummary.csv", package = "stmut"), header = TRUE)
#' arm <- c("1p","3p","3q","4q","5q","8q","9q","10p","10q","11q","13p","13q",
#' "20p","20q","21q","14q","17q")
#' d4 <- data4 %>% filter(arms %in% arm)
#' genes <- d4[,"genes"]
#' rs <- d4[,"gene_row"]
#' gainLoss <- c(1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,1,1,-1,1,1)
#' df7 <- cdt_filt_sort(cdt = cdt,genes = genes,gainLoss = gainLoss,rs=rs)
#' }
cdt_filt_sort <- function(cdt,genes,gainLoss,rs){
  name <- cdt[,c(1,2)]
  cdt <- cdt[,-c(1,2)]

  # convert cdt values to be numeric
  cdt[] <- lapply(cdt, as.numeric) # convert the dataframe to numeric

  wt <- c() # normalize number of genes by divide the max number of genes
  for (n in genes){
    weit <- n/max(genes)
    wt <- append(wt, weit)
  }
  newR <- c() # newR created for sorting the dataframe
  nc <- dim(cdt)[2]
  for (i in c(1:nc)){
    val = 0
    for (k in 1:length(genes)){
      val1 <- gainLoss[k]*(cdt[rs[k],i])*wt[k]  #val <- (-1)*(cdt[rs[1],i])*wt[1]+cdt[rs[2],i]*wt[2]+(-1)*(cdt[rs[3],i])*wt[3] # the sum of weighted values, 3q and 8q has CN gains, 3p and 8p has CN loss, so val = 3q*wt - 3p*wt +8q*wt-8p*wt
      val = val + val1
    }
    newR <- append(newR, val)
  }

  cdt1 <- rbind(cdt,newR)
  cdt1 <- cdt1[,order(cdt1[nrow(cdt1),],decreasing = FALSE)] # reorder by the added row, order() has warning, but works.
  cdt2 <- cdt1[-nrow(cdt1),] # remove the added row
  cdt3 <- cbind(name,cdt2)
  return(cdt3)
}












