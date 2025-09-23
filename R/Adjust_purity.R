##############################################################################################################################################
## source the settings
#source("./Nextcloud/NB_manuscript/Submission_NG/RevisionII//Plots_and_scripts/Custom_scripts/Settings.R")

##############################################################################################################################################
## For most tumors, the purity-ploidy estimate is okay, but for XI003_19079 it needs to be adjusted

# purity <- c()
# ploidy <- c()
#
# for(i in c(tumors.discovery, tumors.validation)){
#
#   if(i %in% tumors.discovery){
#     data.directory <- data.directory.discovery
#   }else{
#     data.directory <- data.directory.validation
#   }
#
#   aceseq.file.name. <- list.files(paste0(data.directory, i, "/ACEseq"), pattern="comb_pro_extra")[1]
#
#   purity. <- Extract.purity.ploidy.from.ACEseq(aceseq.file.name.)$purity
#   ploidy. <- Extract.purity.ploidy.from.ACEseq(aceseq.file.name.)$ploidy
#
#   if(i==tumor.to.adjust.purity){
#     ## re-estimate the purity from the VAF-distribution using normal mixture models
#
#     ## Find mutation file
#     snvs <- list.files(paste0(data.directory, "/", i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10", full.names = T)[1]
#     ## read in the mutation file
#     snvs <- read.vcf(snvs)
#     ## Get VAF, depth and number of variant readds
#     snvs$vcf$VAF <- Extract.info.from.vcf(snvs, info="VAF", type="snvs", mutationcaller="DKFZ")
#
#     snvs <- snvs$vcf
#
#     mixmdl = normalmixEM(snvs$VAF)
#     purity. <- max(mixmdl$mu)*2
#     while(purity.>1){
#       mixmdl = normalmixEM(snvs$VAF)
#       purity. <- max(mixmdl$mu)*2
#     }
#   }
#
#   purity[i] <- purity.
#   ploidy[i] <- ploidy.
# }

#save(purity, ploidy, file=paste0(rdata.directory, "/Purity_ploidy.RData"))
