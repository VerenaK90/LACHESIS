#' Estimate mutation counts at MRCA and ECA.
#'
#' @param clonal.mutation.matrix.per.tumor Matrix with mutations sorted by copy number and allele as returned by Count.clonal.mutations; rows: CN/B state, columns: chromosomes
#' @param segment.length.matrix.per.tumor Matrix with segment lengths sorted by copy number and allele as returned by Count.clonal.mutations; rows: CN/B state, columns: chromosomes
#' @param chromosomes chromosomes of interest. Defaults to 1-22.
#' @param haploid.genome.length size of the genome under consideration. Defaults to 3.3*10^9 (human genome)
#' @param minimal.segment.size What is the minimal segment size to take into consideration? Defaults to 10^7bp.
#' @param max.CN maximal copy number to be considered. Defaults to 4.  Which information is to be retrieved? Possible values are "readcounts", "depth", "VAF", "AA_change". Defaults to "readcounts".
#' @param mutationcaller which mutation caller was used to generate the vcf-files? Defaults to "DKFZ".
#' @return A list containing the estimates for the mutation time at MRCA and ECA along with lower and upper boundaries according to 95% confidence intervals estimated by bootstrapping the segments 1000 times. Earliest.mutation.times returns the density of events earlier than a common ECA. Other list elements report the segments assigned to ECA and MRCA: gains.at.earliest.time - chromosomal gains assigned to the earliest datable time point (NA if all early gains map to a common ECA), gains.not.mapping.to.earliest.time - gains that do not map to MRCA, ECA or the earliest time point, other.gains.mapping.to.earliest.time - gains that can be mapped to the earliest time point or to ECA or to MRCA, gains.at.mrca - gains mapping to the MRCA (may also include gains that can also be mapped to ECA, i.e., not uniquely dateable gains. These pop up in gains.at.mrca.conforming.eca), gains.at.mrca.conforming.eca - gains mapping to the ECA or the MRCA, gains.uniquely.mapped.to.eca - gains mapping to ECA only, gains.not.maping.to.eca.or.mrca - gains that cannot be mapped to the ECA and cannot be mapped to the MRCA, monosomic.states.not.matching.mrca - monosomic states (i.e., non-amplified clonal segments) that do not conform to the MRCA, p.values.mrca - the p values for each segment comparing its density to the MRCA (one-sided negative binomial test), p.values.eca - the p values for each segment comparing its density to the ECA (one-sided negative binomial test)
#' MRCA.ECA.quantification()
#' @export


MRCA.ECA.quantification <- function(clonal.mutation.matrix.per.tumor, segment.length.matrix.per.tumor, chromosomes = c(1:22),
                                    haploid.genome.length = 3.3*10^9, minimal.segment.size = 10^7, max.CN=4){

  ## initialize output values
  mutation.time.mrca = c()
  mutation.time.mrca.lower = c()
  mutation.time.mrca.upper = c()

  earliest.mutation.time = c()
  earliest.mutation.time.lower = c()
  earliest.mutation.time.upper = c()
  gains.at.earliest.time = c()
  gains.not.mapping.to.earliest.time = c()
  other.gains.mapping.to.earliest.time = c()

  mutation.time.eca = c()
  mutation.time.eca.lower = c()
  mutation.time.eca.upper = c()
  gains.at.mrca = c()
  gains.uniquely.mapped.to.eca = c()
  gains.not.maping.to.eca.or.mrca = c()
  monosomic.states.not.matching.mrca = c()
  gains.at.mrca.conforming.eca = c()

  p.values.mrca <- c()
  p.values.eca <- c()

  ## indicator vectors for the copy numbers and B-alleles that were taken into account
  copy.number.indicator <- c(unlist(sapply(1:max.CN, function(x){rep(x, each=x)})))
  B.allele.indicator <- c(unlist(sapply(1:max.CN, function(x){1:x})))

  mrca.mutation.count <- c()
  mrca.genome.length <- c()

    for(j in chromosomes){
      ## Subset the respective chromosome
      tmp.mut.count <- Mutation.time.converter(clonal.mutation.matrix.per.tumor[,j])

      tmp.genome.length <- segment.length.matrix.per.tumor[,j]
      ## Take only lower-order peaks
      tmp.mut.count <- tmp.mut.count[which(B.allele.indicator==1)]
      tmp.genome.length <- tmp.genome.length[which(B.allele.indicator==1)]

      ## Restrict analysis to fragments > minimal.segment.size bp
      tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>minimal.segment.size)]
      tmp.genome.length <- tmp.genome.length[tmp.genome.length>minimal.segment.size]

      if(length(tmp.genome.length)==0){next}

      names(tmp.genome.length) <- paste("chr", j, names(tmp.genome.length), sep="_")

      mrca.mutation.count <- c(mrca.mutation.count, tmp.mut.count)
      mrca.genome.length <- c(mrca.genome.length, tmp.genome.length)

    }
    if(length(mrca.mutation.count)==0){next}

    names(mrca.mutation.count) <- names(mrca.genome.length)

    ## subtract the estimated error of clonal mutations due to spatial sampling
    mutation.time.mrca <- sum(mrca.mutation.count)*haploid.genome.length/sum(mrca.genome.length)*(1- mean(clonal.mutations.false.positives/clonal.mutations.all))

    ## bootstrap upper and lower limits of the mutation time; + subtract the number of false positive mutations due to spatial sampling
    bootstrapped.time <- sapply(1:1000, function(x){
      res <- sample(x=1:length(mrca.mutation.count), size=length(mrca.mutation.count), prob= mrca.genome.length, replace=T)
      res <- sum(mrca.mutation.count[res])*haploid.genome.length/sum(mrca.genome.length[res])
      res <- res - res*rnorm(n=1, mean=mean(clonal.mutations.false.positives/clonal.mutations.all), sd=sd(clonal.mutations.false.positives/clonal.mutations.all))
    })
    mutation.time.mrca.lower <- quantile(bootstrapped.time, 0.025)
    mutation.time.mrca.upper <- quantile(bootstrapped.time, 0.975)


    ## now, test whether the individual fragments match to a joint time point. Test with a negative binomial distribution to account for overdispersion of the
    ## local mutation rate:

    does.fragment.match <- apply(rbind(mrca.mutation.count, mrca.genome.length), 2, function(x){
      if(x[1] < round(sum(mrca.mutation.count))* x[2]/sum(mrca.genome.length)){
        test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)), prob = round(sum(mrca.mutation.count))/(round(sum(mrca.mutation.count))*(1+x[2]/sum(mrca.genome.length))))
      }else{
        test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)), prob = round(sum(mrca.mutation.count))/(round(sum(mrca.mutation.count))*(1+x[2]/sum(mrca.genome.length))), lower.tail = F)
      }
      test})
    haploid.fragment.names <- names(does.fragment.match)
    ## collect the p-values and adjust them later, once contributions from early clonal mutations are accounted for, too
    p.values <- does.fragment.match

    ## now, test whether mutation counts on amplified fragments are consistent with the mrca
    mean.fp <- mean(clonal.mutations.false.positives/clonal.mutations.all)
    eca.mutation.counts <- c()
    eca.genome.length <- c()
    eca.mutation.counts.at.mrca <- c()
    eca.genome.length.at.mrca <- c()
    aneuploid.fragment.name <- c()
    aneuploid.mutation.counts <- c()
    aneuploid.genome.lengths <- c()
    for(j in chromosomes){
      ## Subset the respective chromosome
      tmp.mut.count <- Mutation.time.converter(clonal.mutation.matrix.per.tumor[,j])
      tmp.genome.length <- segment.length.matrix.per.tumor[,j]

      tmp.mut.count <- tmp.mut.count[-which(B.allele.indicator==1)]
      tmp.genome.length <- tmp.genome.length[-which(B.allele.indicator==1)]

      ## Restrict analysis to fragments > minimal.segment.size bp
      tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>minimal.segment.size)]
      tmp.genome.length <- tmp.genome.length[tmp.genome.length>minimal.segment.size]

      if(length(tmp.genome.length)==0){next}

      names(tmp.genome.length) <- paste("chr", j, names(tmp.genome.length), sep="_")
      names(tmp.mut.count) <- names(tmp.genome.length)

      ## does the fragment conform to the mrca? Test with negative binomial distribution

      does.fragment.match <- apply(rbind(tmp.mut.count, tmp.genome.length), 2, function(x){
        if(x[1] <= round(sum(mrca.mutation.count))*(1-mean.fp)* x[2]/sum(mrca.genome.length)){
          test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)*(1-mean.fp)), prob = round(sum(mrca.mutation.count)*(1-mean.fp))/(round(sum(mrca.mutation.count)*(1-mean.fp))*(1+x[2]/sum(mrca.genome.length))))
        }else{
          test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)*(1-mean.fp)), prob = round(sum(mrca.mutation.count)*(1-mean.fp))/(round(sum(mrca.mutation.count)*(1-mean.fp))*(1+x[2]/sum(mrca.genome.length))), lower.tail = F)
        }
        test})
      aneuploid.fragment.name <- c(aneuploid.fragment.name, names(does.fragment.match))
      p.values <- c(p.values, does.fragment.match)

      aneuploid.genome.lengths <- c(aneuploid.genome.lengths, tmp.genome.length)
      aneuploid.mutation.counts <- c(aneuploid.mutation.counts, tmp.mut.count)

    }

    ## Adjust p values and cut off at 0.01. Fragments with p < 0.01 likely arose from a different time point
    adjusted.p.values <- p.adjust(p.values)
    p.values.mrca <- rbind(p.values.mrca, data.frame(p=p.values, adj.p = adjusted.p.values, row.names = names(adjusted.p.values)))
    monosomic.states.not.matching.mrca <- names(adjusted.p.values[haploid.fragment.names][adjusted.p.values[haploid.fragment.names]<0.01])
    gains.at.mrca <- names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]>=0.01])
    eca.mutation.counts.at.mrca <- aneuploid.mutation.counts[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]>=0.01])]
    eca.genome.length.at.mrca <- aneuploid.genome.lengths[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]>=0.01])]
    eca.mutation.counts <- aneuploid.mutation.counts[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]<0.01])]
    eca.genome.length <- aneuploid.genome.lengths[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]<0.01])]

    ## From the non-mapping fragments estimate a time of origin and test whether all of them conform to it.
    ## here, fragments with counts higher than the mrca should be excluded. OK as these are likely outliers
    not.used.for.quantification <- which(eca.mutation.counts*haploid.genome.length/eca.genome.length > mutation.time.mrca[length(mutation.time.mrca)])
    if(length(not.used.for.quantification) > 0){
      mutation.time.eca <- c(mutation.time.eca, sum(eca.mutation.counts[-not.used.for.quantification])*haploid.genome.length/sum(eca.genome.length[-not.used.for.quantification]))
      ## bootstrap upper and lower limits of the mutation time
      if(length(eca.mutation.counts[-not.used.for.quantification])>0){
        bootstrapped.time <- sapply(1:1000, function(x){
          res <- sample(x=1:length(eca.mutation.counts[-not.used.for.quantification]), size=length(eca.mutation.counts[-not.used.for.quantification]), prob= eca.genome.length[-not.used.for.quantification], replace=T)
          sum(eca.mutation.counts[-not.used.for.quantification][res])*haploid.genome.length/sum(eca.genome.length[-not.used.for.quantification][res])
        })
        mutation.time.eca.lower <- quantile(bootstrapped.time, 0.025)
        mutation.time.eca.upper <- quantile(bootstrapped.time, 0.975)
      }else{
        mutation.time.eca.lower <- NA
        mutation.time.eca.upper <- NA
      }
      gains.not.maping.to.eca.or.mrca <- names(not.used.for.quantification)
      eca.mutation.counts <- eca.mutation.counts[-not.used.for.quantification]
      eca.genome.length <- eca.genome.length[-not.used.for.quantification]
    }else{
      mutation.time.eca <- sum(eca.mutation.counts)*haploid.genome.length/sum(eca.genome.length)
      ## bootstrap upper and lower limits of the mutation time
      if(length(eca.mutation.counts)>0){
        bootstrapped.time <- sapply(1:1000, function(x){
          res <- sample(x=1:length(eca.mutation.counts), size=length(eca.mutation.counts), prob= eca.genome.length, replace=T)
          sum(eca.mutation.counts[res])*haploid.genome.length/sum(eca.genome.length[res])
        })
        mutation.time.eca.lower <- quantile(bootstrapped.time, 0.025)
        mutation.time.eca.upper <- quantile(bootstrapped.time, 0.975)

      }else{
        mutation.time.eca.lower <- NA
        mutation.time.eca.upper <- NA
      }
    }




    ## now, test whether the individual fragments match with the eca estimate:
    if(length(eca.genome.length)>0){
      does.fragment.match <- apply(rbind(eca.mutation.counts, eca.genome.length), 2, function(x){
        if(length(eca.mutation.counts)==1){
          x[1] <- round(x[1])
        }
        if(x[1]<1){
          x[1] <- ceiling(x[1])
        }
        if(x[1]==0 & sum(eca.mutation.counts)==0){
          return(1)
        }
        total.eca.counts <- if(sum(eca.mutation.counts) < 1){ceiling(sum(eca.mutation.counts))}else{sum(eca.mutation.counts)}
        if(x[1] <= total.eca.counts* x[2]/sum(eca.genome.length)){
          test <- pnbinom(q = x[1], size = total.eca.counts, prob = total.eca.counts/(total.eca.counts*(1+x[2]/sum(eca.genome.length))))
        }else{
          test <- pnbinom(q = x[1], size = total.eca.counts, prob = total.eca.counts/(total.eca.counts*(1+x[2]/sum(eca.genome.length))), lower.tail = F)
        }
        test
      })
      if(round(sum(eca.mutation.counts))==0){
        does.fragment.match <- rep(T, length(eca.mutation.counts))
        names(does.fragment.match) <- names(eca.mutation.counts)
        p.values.eca <- rbind(p.values.eca, data.frame(p=rep(NA, length(does.fragment.match)),
                                                       adj.p = rep(NA, length(does.fragment.match)), row.names = names(does.fragment.match)))

      }else{

        p.values <- does.fragment.match
        does.fragment.match <- p.adjust(does.fragment.match)
        p.values.eca <- rbind(p.values.eca, data.frame(p=p.values, adj.p = does.fragment.match, row.names = names(does.fragment.match)))

        does.fragment.match <- sapply(does.fragment.match, function(x){
          if(x < 0.01){
            F
          }else{
            T
          }
        })
      }

      ## Report fragments that conform and don't conform to a common origin
      gains.uniquely.mapped.to.eca <- names(does.fragment.match[does.fragment.match])
      gains.not.maping.to.eca.or.mrca <- names(does.fragment.match[!does.fragment.match])

      ## If some fragments do not match, check whether they request an even earlier time point and if so, determine
      if(length(does.fragment.match)>0 & any(!does.fragment.match)){
        fragments.before.eca <- which(eca.mutation.counts[!does.fragment.match]*haploid.genome.length/eca.genome.length[!does.fragment.match] < mutation.time.eca[length(mutation.time.eca)])

        if(length(fragments.before.eca)>0){

          earliest.mutation.count <- eca.mutation.counts[!does.fragment.match][fragments.before.eca]
          earliest.genome.length <- eca.genome.length[!does.fragment.match][fragments.before.eca]

          earliest.mutation.time <-  sum(eca.mutation.counts[!does.fragment.match][fragments.before.eca])*haploid.genome.length/
                                        sum(eca.genome.length[!does.fragment.match][fragments.before.eca])


          ## bootstrap upper and lower limits of the mutation time
          if(length(earliest.mutation.count)>0){
            bootstrapped.time <- sapply(1:1000, function(x){
              res <- sample(x=1:length(earliest.mutation.count), size=length(earliest.mutation.count), prob= earliest.genome.length, replace=T)
              sum(earliest.mutation.count[res])*haploid.genome.length/sum(earliest.genome.length[res])
            })
            earliest.mutation.time.lower <- quantile(bootstrapped.time, 0.025)
            earliest.mutation.time.upper <- quantile(bootstrapped.time, 0.975)
          }else{
            earliest.mutation.time.lower <- NA
            earliest.mutation.time.upper <- NA
          }


          does.fragment.match <- apply(rbind(earliest.mutation.count, earliest.genome.length), 2, function(x){
            if(length(earliest.mutation.count)==1){
              x[1] <- round(x[1])
            }
            if(x[1]==0 & sum(earliest.mutation.count)==0){
              return(1)
            }
            if(x[1] <= round(sum(earliest.mutation.count))* x[2]/sum(eca.genome.length)){
              test <- pnbinom(q = x[1], size = round(sum(earliest.mutation.count)), prob = round(sum(earliest.mutation.count))/(round(sum(earliest.mutation.count))*(1+x[2]/sum(earliest.genome.length))))
            }else{
              test <- pnbinom(q = x[1], size = round(sum(earliest.mutation.count)), prob = round(sum(earliest.mutation.count))/(round(sum(earliest.mutation.count))*(1+x[2]/sum(earliest.genome.length))), lower.tail = F)
            }
            test
          })
          does.fragment.match <- p.adjust(does.fragment.match)
          does.fragment.match <- sapply(does.fragment.match, function(x){
            if(x < 0.01){
              F
            }else{
              T
            }
          })
          ## Report fragments that conform and don't conform to a common origin
          gains.at.earliest.time <- names(does.fragment.match[does.fragment.match])
          gains.not.mapping.to.earliest.time <- names(does.fragment.match[!does.fragment.match])
        }



      }

      ## Finally, does a gain that matches to the MRCA also agree with eca?
      if(length(eca.mutation.counts.at.mrca)>0){
        does.fragment.match <- apply(rbind(eca.mutation.counts.at.mrca, eca.genome.length.at.mrca), 2, function(x){
          if(length(eca.mutation.counts.at.mrca)==1){
            x[1] <- round(x[1])
          }
          if(x[1] < 1){
            x[1] <- round(x[1])
          }
          if(x[1]==0 & round(sum(eca.mutation.counts) + x[1])==0){return(T)}
          if(x[1] <= round((sum(eca.mutation.counts) + x[1]))* x[2]/(sum(eca.genome.length)+x[2]) ){
            test <- pnbinom(q = x[1], size = round(sum(eca.mutation.counts) + x[1]), prob = round(sum(eca.mutation.counts) + x[1])/(round(sum(eca.mutation.counts) + x[1])*(1+x[2]/(sum(eca.genome.length) + x[2]))))
          }else{
            test <- pnbinom(q = x[1], size = round(sum(eca.mutation.counts) + x[1]), prob = round(sum(eca.mutation.counts) + x[1])/(round(sum(eca.mutation.counts) + x[1])*(1+x[2]/(sum(eca.genome.length) + x[2]))), lower.tail = F)
          }
          test
        })

        p.values <- does.fragment.match
        does.fragment.match <- p.adjust(does.fragment.match)
        p.values.eca <- rbind(p.values.eca, data.frame(p=p.values, adj.p = does.fragment.match, row.names = names(does.fragment.match)))

        does.fragment.match <- sapply(does.fragment.match, function(x){
          if(x < 0.01){
            F
          }else{
            T
          }
        })
        gains.at.mrca.conforming.eca <- names(does.fragment.match[does.fragment.match])
      }

    }

    return(list(mutation.time.mrca = mutation.time.mrca,
                mutation.time.mrca.lower = mutation.time.mrca.lower,
                mutation.time.mrca.upper = mutation.time.mrca.upper,

                earliest.mutation.time = earliest.mutation.time,
                earliest.mutation.time.lower = earliest.mutation.time.lower,
                earliest.mutation.time.upper = earliest.mutation.time.upper,
                gains.at.earliest.time = gains.at.earliest.time,
                gains.not.mapping.to.earliest.time = gains.not.mapping.to.earliest.time,
                other.gains.mapping.to.earliest.time = other.gains.mapping.to.earliest.time,

                mutation.time.eca = mutation.time.eca,
                mutation.time.eca.lower = mutation.time.eca.lower,
                mutation.time.eca.upper = mutation.time.eca.upper,
                gains.at.mrca = gains.at.mrca,
                gains.uniquely.mapped.to.eca = gains.uniquely.mapped.to.eca,
                gains.not.maping.to.eca.or.mrca = gains.not.maping.to.eca.or.mrca,
                monosomic.states.not.matching.mrca = monosomic.states.not.matching.mrca,
                gains.at.mrca.conforming.eca = gains.at.mrca.conforming.eca,

                p.values.mrca = p.values.mrca,
                p.values.eca = p.values.eca))

}
