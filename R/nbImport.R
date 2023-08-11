#' Combine CNVs and SNVs
#' @description
#' Merges CNVs and SNVs into a single data.frame. Each variant is assigned to its corresponding copy number segment and status.
#' @param cnv CNV data from \code{\link{readCNV}}
#' @param snv SNV data from \code{\link{readVCF}}
#' @examples
#' snvs <- system.file("extdata", "NBE15/NBE15.snvs.dkfz.tsv.gz", package = "NBevolution")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15/NBE15.aceseq.tsv.gz", package = "NBevolution")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data)
#' @return a data.table
#' @export

nbImport <- function(cnv, snv){

  colnames(cnv)[1:3] <- c("chrom", "start", "end")
  data.table::setDT(x = cnv, key = c("chrom", "start", "end"))

  colnames(snv)[1:2] <- c("chrom", "start")
  snv[,end := start]
  data.table::setDT(x = snv, key = c("chrom", "start", "end"))

  sv <- data.table::foverlaps(x = snv, y = cnv, type = "within")
  colnames(sv)[which(colnames(sv) == "i.start")] <- "snv_start"
  colnames(sv)[which(colnames(sv) == "i.end")] <- "snv_end"
  attr(sv, "cnv") <- cnv
  sv
}


#' Plot VAF distribution per copy number
#' @description
#' Visualize results from  \code{\link{nbImport}}
#' @param nb output generated from \code{\link{nbImport}}
#' @param ref_build Default `hg19`. Can be `hg18`, `hg19` or `hg38`
#' @examples
#' snvs <- system.file("extdata", "NBE15/NBE15.snvs.dkfz.tsv.gz", package = "NBevolution")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15/NBE15.aceseq.tsv.gz", package = "NBevolution")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data)
#' plotNB(nb)
#' @export
#' @importFrom graphics abline axis box grid hist mtext par rect text title

plotNB <- function(nb, ref_build = "hg19"){

  segs <- attr(nb, "cnv")
  segs[order(chrom, start)]
  segs <- segs[order(chrom, start)]
  colnames(segs)[1:3] <- c("Chromosome", "Start_Position", "End_Position")
  segs <- .transformSegments(segmentedData = segs, build = "hg19")

  contig_lens <- cumsum(.getContigLens(build = ref_build))

  lo_mat <- matrix(data = c(rep(1, 3), 2:4), nrow = 2, byrow = TRUE)
  graphics::layout(mat = lo_mat, heights = c(3, 2))
  par(mar = c(3, 4, 4, 3))
  plot(NA, ylim = c(0, 4), xlim = c(0, max(contig_lens)), axes = FALSE, xlab = NA, ylab = NA)
  abline(h = 1:4, v = contig_lens, lty = 2, col = "gray", lwd = 0.4)
  rect(xleft = segs$Start_Position_updated, ybottom = segs$TCN-0.1, xright = segs$End_Position_updated, ytop = segs$TCN+0.1, col = ifelse(segs$TCN > 2, "#16a085", "#7f8c8d"), border = NA, lty = 3)
  axis(side = 1, at = contig_lens, labels = 1:24)
  axis(side = 2, at = 0:4, labels = 0:4, las = 2)
  mtext(text = "Total CN", side = 2, line = 2)
  title(main = attr(nb, "t.sample"))

  for(tcn in rev(split(nb, nb$TCN))){
    par(mar = c(3, 4, 3, 1))
    hist(tcn[,t_vaf], breaks = 100, xlim = c(0, 1), xlab = NA, ylab = NA,  border = NA, col = "#34495e", main = NA)
    title(main = paste0("Total CN:", tcn[1, TCN]), cex.main = 1.2)
    mtext(text = "No. of SNVs", side = 2, line = 2.5, cex = 0.7)
  }
}


#contig lengths for hg19, hg38 and hg18
.getContigLens <- function(build = "hg19"){

  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }

  chr.lens

}

#--- Change segment sizes into linear scale
.transformSegments <- function(segmentedData, build = 'hg19'){

  build.opts <- c('hg19', 'hg18', 'hg38')

  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }

  #Get chr lens
  chr.lens <- .getContigLens(build = build)

  segmentedData[,Start_Position := as.numeric(as.character(Start_Position))]
  segmentedData[,End_Position := as.numeric(as.character(End_Position))]

  #Replace chr x and y with numeric value (23 and 24) for better ordering
  segmentedData$Chromosome <- gsub(pattern = 'chr', replacement = '', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome <- gsub(pattern = 'X', replacement = '23', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome <- gsub(pattern = 'Y', replacement = '24', x = segmentedData$Chromosome, fixed = TRUE)

  segmentedData$Chromosome <- factor(x = segmentedData$Chromosome, levels = 1:24, labels = 1:24)

  segmentedData <- segmentedData[order(Chromosome, Start_Position, decreasing = FALSE)]

  seg.spl <- split(segmentedData, segmentedData$Chromosome)

  seg.spl.transformed <- seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated <- seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated <- seg.spl.transformed$End_Position
  }

  chr.lens.sumsum <- cumsum(chr.lens)

  for(i in 2:length(seg.spl)){

    x.seg <- seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated <- x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated <- x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed <- rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }

  return(seg.spl.transformed)
}

