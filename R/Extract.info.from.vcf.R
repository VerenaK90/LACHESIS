#' Extracts information from a vcf file
#'
#' @param vcf Mutation file in VCF representation as a list (as returned by read.vcf from package bedR).
#' @param info Which information is to be retrieved? Possible values are "readcounts", "depth", "VAF", "AA_change". Defaults to "readcounts".
#' @param type Is the VCF file reporting snvs or indels? Possible values are "snvs" and "indel". Defaults to "snvs".
#' @param mutationcaller which mutation caller was used to generate the vcf-files? Defaults to "DKFZ".
#' @return purity and ploidy.
#' Extract.info.from.vcf()
#' @export


Extract.info.from.vcf <- function(vcf, info="readcounts", type="snvs", mutationcaller="DKFZ"){

  if(names(vcf)[1]!="header"){
    vcf <- list(header=c(), vcf=vcf)
  }
  if(mutationcaller=="DKFZ"){
    if(type=="snvs"){
      if(info=="readcounts"){
        readcounts <- t(sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="DP4=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x <- as.numeric(strsplit(x, split=",")[[1]])
          x <- unname(x)
          c(sum(x[c(1,2)]), sum(x[c(3,4)]))
        }))
        rownames(readcounts) <- paste(vcf$vcf$GENE, vcf$vcf$POS, sep=".")
        colnames(readcounts) <- c("REF", "ALT")

        return(readcounts)
      }
      if(info=="depth"){
        depth <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="DP4=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x <- as.numeric(strsplit(x, split=",")[[1]])
          x <- unname(x)
          sum(x)
        })

        return(depth)
      }
      if(info=="VAF"){
        vaf <- t(sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="DP4=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
          x <- as.numeric(strsplit(x, split=",")[[1]])
          x <- unname(x)
          c(sum(x[c(3,4)])/ sum(x))
        }))
        return(vaf)
      }
      if(info=="AA_change"){
        aa_change <- sapply(vcf$vcf$ANNOVAR_TRANSCRIPTS, function(x){
          x <- strsplit(x, split="p.")[[1]][2]
          x <- strsplit(x, split=",")[[1]][1]
        })
        return(aa_change)
      }
      if(info=="clonality"){
        clonality <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="CLS=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(clonality)
      }
      if(info=="Allele"){
        allele <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="MutCN= ")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(Allele)
      }
      if(info=="MajCN"){
        MajCN <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="MajCN=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(MajCN)
      }
      if(info=="MinCN"){
        MinCN <- sapply(vcf$vcf$INFO, function(x){
          x <- strsplit(x, split="MinCN=")[[1]][2]
          x <- strsplit(x, split=";")[[1]][1]
        })
        return(MinCN)
      }
    }else if(type=="indel"){

      if(info=="readcounts"){
        readcounts <- t(sapply(vcf$vcf$FORMAT_INFO, function(x){
          x <- strsplit(x, split=":")[[1]][c(5,6)]
          x <- as.numeric(x)
          x <- unname(x)
          x
        }))
        rownames(readcounts) <- paste(vcf$vcf$GENE, vcf$vcf$POS, sep=".")
        colnames(readcounts) <- c("REF", "ALT")
        return(readcounts)
      }
      if(info=="depth"){
        depth <- sapply(vcf$vcf$FORMAT_INFO, function(x){
          x <- strsplit(x, split=":")[[1]][c(5,6)]
          x <- as.numeric(x)
          x <- unname(x)
          x
          sum(x)
        })

        return(depth)
      }
      if(info=="VAF"){
        vaf <- t(sapply(vcf$vcf$FORMAT_INFO, function(x){
          x <- strsplit(x, split=":")[[1]][c(5,6)]
          x <- as.numeric(x)
          x <- unname(x)
          x[6]/sum(x)
        }))
        return(vaf)
      }
      if(info=="AA_change"){
        aa_change <- sapply(vcf$vcf$ANNOVAR_TRANSCRIPTS, function(x){
          x <- strsplit(x, split="p.")[[1]][2]
          x <- strsplit(x, split=",")[[1]][1]
        })
        return(aa_change)
      }


    }
  }


}
