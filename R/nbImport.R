#' Combine CNVs and SNVs
#' @description
#' Merges CNVs and SNVs into a single data.table. Each variant is assigned to its corresponding copy number segment and status.
#' @param cnv CNV data from \code{\link{readCNV}}.
#' @param snv SNV data from \code{\link{readVCF}}.
#' @param purity tumor cell content.
#' @param ploidy average copy number in the tumor sample.
#' @param sig.assign Logical. If TRUE, each variant will be assigned to a mutational signature.
#' @param assign.method Method to assign signatures: "max" to assign the signature with the highest probability, "sample" to randomly assign based on signature probabilities.
#' @param ID sample name.
#' @param sig.file File path to the SigAssignment output file, typically named "Decomposed_MutationType_Probabilities.txt".
#' @param sig.select A character vector of specific signatures to include in the analysis (e.g., c("SBS1", "SBS5", "SBS40") to focus on clock-like mutational processes).
#' @param min.p Numeric. The minimum probability threshold from the SigAssignment output that a variant must meet to be considered as matching a specific signature.
#' @param ref.build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or `hg38`.
#' @param seed Integer. Can be user-specified or an automatically generated random seed, it will be documented in the log file.
#'
#' @examples
#' # Example using all variants from vcf file
#' snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
#'
#' # Example using variants associated with specific SBS mutational signatures from vcf file
#' snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' sig.filepath <- system.file("extdata", "NBE15_Decomposed_MutationType_Probabilities.txt", package = "LACHESIS")
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51, sig.assign = TRUE, ID = "NBE15", sig.file = sig.filepath, sig.select = c("SBS1", "SBS5", "SBS40a", "SBS18"))
#' @seealso \code{\link{plotNB}}
#' @return a data.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Biostrings getSeq
#' @export

nbImport <- function(cnv = NULL, snv = NULL, purity = NULL, ploidy = NULL,
                     sig.assign = FALSE, assign.method = "sample", ID = NULL,
                     sig.file = NULL, sig.select = NULL, min.p = NULL,
                     ref.build = "hg19", seed = NULL) {
     end <- start <- sequence_context <- chrom <- i.end <- i.start <- TCN <- NULL

    if (any(is.null(cnv), is.null(snv))) {
        stop("Missing snv and cnv inputs!")
    }
    if (any(is.null(purity), is.null(ploidy))) {
        stop("Missing purity and ploidy inputs!")
    }

    if (is.null(seed)) {
        seed <- sample.int(.Machine$integer.max, 1)
    }

    colnames(cnv)[c(1, 2, 3)] <- c("chrom", "start", "end")
    data.table::setDT(x = cnv, key = c("chrom", "start", "end"))
    colnames(snv)[c(1, 2)] <- c("chrom", "start")
    snv[, end := start]
    data.table::setDT(x = snv, key = c("chrom", "start", "end"))

    sv <- data.table::foverlaps(x = snv, y = cnv, type = "within")

    if (nrow(sv) == 0) {
        stop("No overlapping SNVs found within CNV regions")
    }

    if (nrow(sv[is.na(start)])) {
        warning("Removed ", nrow(sv[is.na(start)]), " variants with no copy number overlaps")
        sv <- sv[!is.na(start)]
    }

    if (sig.assign == TRUE) {
        t.sample <- attributes(sv)$t.sample
        assign.result <- .assign_signatures(
            sv, sig.file, assign.method, ID,
            sig.select, min.p, ref.build, seed
        )
        sv <- assign.result$sv
        sig.colors <- assign.result$sig.colors
        attr(sv, "t.sample") <- t.sample
        attr(sv, "sig.colors") <- sig.colors
    }

    # Make columns more intuitive
    colnames(sv)[which(colnames(sv) == "i.start")] <- "snv_start"
    colnames(sv)[which(colnames(sv) == "i.end")] <- "snv_end"
    colnames(sv)[which(colnames(sv) == "start")] <- "cn_start"
    colnames(sv)[which(colnames(sv) == "end")] <- "cn_end"
    attr(sv, "cnv") <- cnv
    attr(sv, "purity") <- as.numeric(purity)
    attr(sv, "ploidy") <- as.numeric(ploidy)
    sv
}

.assign_signatures <- function(sv = NULL, sig.file = NULL,
                               assign.method = "sample", ID = NULL,
                               sig.select = NULL, min.p = NULL, ref.build = NULL,
                               seed = NULL) {
  strand <- ref <- sequence_context <- chrom <- i.start <- i.end <- Sample <- MutationType <- alt <- NULL

    if (is.null(sv)) {
        stop("Missing 'sv' input data!")
    }

    if (is.null(sig.file)) {
        stop("Missing 'SigAssignment' input data!")
    }

    sig.data <- data.table::fread(sig.file)

    data.table::setnames(sig.data, c("Sample Names"), c("Sample"))

    sbs.cols <- grep("^SBS", names(sig.data), value = TRUE)

    if (!"sequence_context" %in% colnames(sv)) {
        genome <- switch(ref.build,
            "hg18" = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
            "hg19" = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
            "hg38" = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        )

        # Mapping purine bases to reverse strand ("-")
        sv[, strand := ifelse(ref %in% c("A", "G"), "-", "+")]

        # Extracting 3-base sequence context (ref base on "-" will be reverse-complemented)
        sv[, sequence_context := as.character(Biostrings::getSeq(
            genome,
            names = paste0("chr", chrom),
            start = i.start - 1,
            end = i.end + 1,
            strand = strand
        ))]

        sv[, strand := NULL]
    }

    sv[, Sample := ID]

    # Constructing MutationType (alt base on "-" will be reverse-complemented)
    sv[, MutationType := {
        ctx <- sequence_context
        corrected.alt <- alt
        corrected.alt[alt == "A" & ref %in% c("A", "G")] <- "T"
        corrected.alt[alt == "C" & ref %in% c("A", "G")] <- "G"
        corrected.alt[alt == "G" & ref %in% c("A", "G")] <- "C"
        corrected.alt[alt == "T" & ref %in% c("A", "G")] <- "A"
        paste0(
            substr(ctx, 1, 1), "[", substr(ctx, 2, 2), ">", corrected.alt, "]",
            substr(ctx, 3, 3)
        )
    }]

    sv <- merge(sv, sig.data,
        by = c("Sample", "MutationType"), all.x = TRUE,
        all.y = FALSE
    )

    if (assign.method == "sample") {
        set.seed(seed)
        tmp <- sv[,
            {
                probs <- as.numeric(unlist(.SD, use.names = FALSE))

                if (sum(probs, na.rm = TRUE) == 0 || anyNA(probs)) {
                    Signature <- as.character(NA)
                    Probability <- as.numeric(NA)
                } else {
                    probs <- probs / sum(probs)
                    sampled.index <- sample(seq_along(probs), 1, prob = probs)
                    Signature <- sbs.cols[sampled.index]
                    Probability <- probs[sampled.index]
                }

                list(Signature = Signature, Probability = Probability)
            },
            .SDcols = sbs.cols,
            by = seq_len(nrow(sv)),
        ]
        sv[, `:=`(Signature = tmp$Signature, Probability = tmp$Probability)]
    } else {
        sv <- sv[,
            {
                if (any(is.na(unlist(.SD, use.names = FALSE)))) {
                    list(
                        Signature = as.character(NA),
                        Probability = as.numeric(NA)
                    )
                } else {
                    max.p.sig <- which.max(unlist(.SD, use.names = FALSE))
                    list(
                        Signature = sbs.cols[max.p.sig],
                        Probability = .SD[[max.p.sig]]
                    )
                }
            },
            .SDcols = sbs.cols,
            by = seq_len(nrow(sv))
        ]
    }
    sv <- sv[!is.na(Probability), ]

    if (!is.null(min.p)) {
        sv <- sv[Probability >= min.p]
    }

    if (!is.null(sig.select)) {
        sv <- sv[Signature %in% sig.select]
        sig.number <- length(sig.select)
        sig.colors <- setNames(.get_sig_colors(sig.number), sig.select)
    } else {
        sig.options <- unique(sv$Signature)
        sig.number <- length(sig.options)
        sig.colors <- setNames(.get_sig_colors(sig.number), sig.options)
    }


    sv[, "Sample" := NULL]

    return(list(sv = sv, sig.colors = sig.colors))
}

.get_sig_colors <- function(n, palette = "Set3", max.colors = 12) {
    base.colors <- RColorBrewer::brewer.pal(min(max.colors, n), palette)
    if (n > max.colors) {
        colorRampPalette(base.colors)(n)
    } else {
        base.colors
    }
}

#' Plot VAF distribution per copy number
#' @description
#' Visualizes results from  \code{\link{nbImport}}. Top plot, measured copy numbers along the genome; bottom plots, VAF histograms of SNVs stratified by copy number and minor/major allele count.
#' @param nb output generated from \code{\link{nbImport}}.
#' @param snvClonality output generated from \code{\link{estimateClonality}}.
#' @param ref.build Reference genome. Default `hg19`. Can be `hg18`, `hg19` or `hg38`.
#' @param min.cn maximum copy number to be included in the plotting. Defaults to 2.
#' @param max.cn maximum copy number to be included in the plotting. Defaults to 4.
#' @param nb.col.abline optional, the color code for the lines depicting clonality in the VAF histograms.
#' @param nb.col.cn.2 optional, the color code for tcn = 2 in the CNV plot.
#' @param nb.col.cn optional, the color code for other copy numbers in the CNV plot.
#' @param nb.col.hist optional, the color code for bars in the VAF histograms.
#' @param nb.border, optional, the line color in the VAF histograms.
#' @param nb.breaks optional, the number of bins in the histograms.
#' @param samp.name Sample name. Optional. Default NULL
#' @param output.file optional, will save the plot.
#' @param sig.show plot stratified VAF histogram with assigned mutational signatures.
#' @param ... further arguments and parameters passed to other LACHESIS functions.
#' @examples
#' # Example using all variants from vcf file
#' snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51)
#' cl_muts <- clonalMutationCounter(nb)
#' norm_muts <- normalizeCounts(cl_muts)
#' mrca <- MRCA(norm_muts)
#' snvClonality <- estimateClonality(nbObj = nb, mrcaObj = mrca, ID = "NBE15", purity = 1)
#' plotNB(nb = nb, snvClonality = snvClonality)
#'
#' # Example using variants assosciated with specific SBS mutational signatures from vcf file
#' snvs <- system.file("extdata", "NBE15", "snvs_NBE15_somatic_snvs_conf_8_to_10.vcf", package = "LACHESIS")
#' s_data <- readVCF(vcf = snvs, vcf.source = "dkfz")
#' aceseq_cn <- system.file("extdata", "NBE15", "NBE15_comb_pro_extra2.51_1.txt", package = "LACHESIS")
#' c_data <- readCNV(aceseq_cn)
#' sig.filepath <- system.file("extdata", "NBE15_Decomposed_MutationType_Probabilities.txt", package = "LACHESIS")
#' nb <- nbImport(cnv = c_data, snv = s_data, purity = 1, ploidy = 2.51, sig.assign = TRUE, ID = "NBE15", sig.file = sig.filepath)
#' cl_muts <- clonalMutationCounter(nb)
#' norm_muts <- normalizeCounts(cl_muts)
#' mrca <- MRCA(norm_muts)
#' snvClonality <- estimateClonality(nbObj = nb, mrcaObj = mrca, ID = "NBE15", purity = 1)
#' plotNB(nb = nb, snvClonality = snvClonality, sig.show = TRUE)
#'
#' @export
#' @importFrom graphics abline axis box grid hist mtext par rect text title
#' @import ggplot2
#' @import gridExtra
#' @import data.table

plotNB <- function(nb = NULL, snvClonality = NULL, ref.build = "hg19", min.cn = 2,
                   max.cn = 4, nb.col.abline = "gray70", nb.col.cn.2 = "#7f8c8d",
                   nb.col.cn = "#16a085", nb.col.hist = "#34495e", nb.border = NA,
                   nb.breaks = 100, samp.name = NULL, output.file = NULL,
                   sig.show = FALSE, ...) {
    chrom <- start <- t_vaf <- Clonality <- Signature <- End <- End_Position_updated <- Start <- Start_Position_updated <- TCN <- NULL

    if (is.null(nb)) {
        stop("Missing input. Please provide the output generated by nbImport")
    }

    if (is.null(snvClonality)) {
        stop("Missing input. Please provide the output generated by estimateClonality")
    }

    if (max.cn <= min.cn) {
        stop("max.cn must be larger than min.cn")
    }


    sig.colors <- attr(nb, "sig.colors")
    purity <- attr(nb, "purity")
    clonality_colors <- c(
        "Precnv" = "#66c2a5", "Postcnv" = "#fc8d62",
        "C" = "#8da0cb", "SC" = "#e78ac3"
    )

    segs <- attr(nb, "cnv")
    segs <- segs[order(chrom, start)]
    colnames(segs)[c(1, 2, 3)] <- c("Chromosome", "Start_Position", "End_Position")
    segs <- .transformSegments(segmentedData = segs, build = ref.build)

    contig_lens <- cumsum(.getContigLens(build = ref.build))
    contig_lens_dt <- data.table(
        Chromosome = c(seq_len(22), "X", "Y"),
        mid = c(
            contig_lens[1] / 2,
            (contig_lens[-length(contig_lens)] +
                contig_lens[-1]) / 2
        )
    )

    label_subset <- contig_lens_dt[seq(1, .N, by = 2)]

    cnv_plot <- ggplot() +
        geom_rect(
            data = segs, aes(
                xmin = Start_Position_updated,
                xmax = End_Position_updated, ymin = TCN - 0.1,
                ymax = TCN + 0.1,
                fill = factor(ifelse(TCN == 2, nb.col.cn.2, nb.col.cn))
            ),
            color = nb.border, linetype = "dotted"
        ) +
        scale_fill_identity() +
        geom_hline(
            yintercept = seq_len(max.cn), linetype = "dashed",
            color = nb.col.abline, size = 0.3
        ) +
        geom_vline(
            xintercept = contig_lens, linetype = "dashed",
            color = nb.col.abline, size = 0.3
        ) +
        scale_x_continuous(
            breaks = label_subset$mid,
            labels = label_subset$Chromosome, expand = c(0, 0)
        ) +
        scale_y_continuous(breaks = c(0, seq_len(max.cn))) +
        labs(
            x = "Chromosome", y = "Total CN",
            title = ifelse(is.null(samp.name), attr(nb, "t.sample"), samp.name)
        ) +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 9),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    # Clonality histograms
    snvClonality <- snvClonality[TCN >= min.cn & TCN <= max.cn]
    snvClonality[, TCN := factor(TCN, levels = seq_len(max.cn))]
    snvClonality_split_TCN <- split(snvClonality, by = "TCN")

    clonality_plots <- list()
    for (cn in names(snvClonality_split_TCN)) {
        snvClonality_split_TCN_B <- split(
            snvClonality_split_TCN[[cn]],
            snvClonality_split_TCN[[cn]]$B
        )
        for (b in names(snvClonality_split_TCN_B)) {
            tcn <- snvClonality_split_TCN_B[[b]]
            tcn$Clonality <- factor(tcn$Clonality,
                levels = c("Precnv", "Postcnv", "C", "SC")
            )
            if (nrow(tcn) == 0) next
            max_count <- max(hist(tcn$t_vaf, breaks = nb.breaks, plot = FALSE)$counts)
            p_clonality <- ggplot(tcn, aes(x = t_vaf, fill = Clonality)) +
                geom_histogram(
                    bins = nb.breaks, color = NA,
                    position = "stack", show.legend = T
                ) +
                scale_fill_manual(
                    values = clonality_colors,
                    labels = c(
                        "Precnv" = "Clonal\n- Pre-CNV",
                        "Postcnv" = "Clonal\n- Post-CNV",
                        "C" = "Clonal\n- NOS",
                        "SC" = "Subclonal"
                    ), drop = F
                ) +
                scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
                labs(
                    x = "VAF", y = "No. of SNVs",
                    title = paste0("CN:", cn, " (", as.numeric(cn) -
                        as.numeric(b), ":", b, ")")
                ) +
                theme_classic()
            if (!is.null(purity)) {
                expected_vafs <- .expectedClVAF(CN = as.numeric(cn), purity = purity)
                p_clonality <- p_clonality +
                    geom_vline(
                        xintercept = expected_vafs,
                        linetype = "dashed", color = "black"
                    )
            }
            clonality_plots[[paste0(cn, "_", b)]] <- p_clonality
        }
    }

    # Signature histograms
    signature_plots <- NULL
    if (sig.show) {
        signature_plots <- list()
        for (cn in names(snvClonality_split_TCN)) {
            snvClonality_split_TCN_B <- split(
                snvClonality_split_TCN[[cn]],
                snvClonality_split_TCN[[cn]]$B
            )
            for (b in names(snvClonality_split_TCN_B)) {
                tcn <- snvClonality_split_TCN_B[[b]]
                if (nrow(tcn) == 0) next
                p_signature <- ggplot(tcn, aes(x = t_vaf, fill = Signature)) +
                    geom_histogram(
                        bins = nb.breaks, color = NA,
                        position = "stack", show.legend = T
                    ) +
                    scale_fill_manual(values = sig.colors, drop = F) +
                    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
                    labs(
                        x = "VAF", y = "No. of SNVs",
                        title = paste0("CN:", cn, " (", as.numeric(cn) -
                            as.numeric(b), ":", b, ")")
                    ) +
                    theme_classic() +
                    theme(
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 9)
                    )
                if (!is.null(purity)) {
                    expected_vafs <- .expectedClVAF(CN = as.numeric(cn), purity = purity)
                    p_signature <- p_signature +
                        geom_vline(
                            xintercept = expected_vafs,
                            linetype = "dashed", color = "black"
                        )
                }
                signature_plots[[paste0(cn, "_", b)]] <- p_signature
            }
        }
    }

    if (!is.null(output.file)) {
        pdf(output.file, width = 7, height = 9)
    }

    # Copy number plot and clonality histograms
    clonality_plot <- do.call(
        .grid_arrange_shared_legend,
        c(clonality_plots,
            ncol = 2,
            nrow = ceiling(length(clonality_plots) / 2)
        )
    )
    first_page <- gridExtra::arrangeGrob(cnv_plot, clonality_plot,
        ncol = 1,
        heights = c(0.4, 0.6)
    )
    grid::grid.draw(first_page)

    # Optional signature histograms
    if (sig.show && length(signature_plots) > 0) {
        do.call(gridExtra::grid.arrange, c(signature_plots, ncol = 2))
    }

    if (!is.null(output.file)) {
        dev.off()
    }
}

# Contig lengths for hg19, hg38 and hg18
.getContigLens <- function(build = "hg19") {
    if (build == "hg19") {
        chr.lens <- c(
            249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 81195210, 78077248, 59128983, 63025520, 48129895,
            51304566, 155270560, 59373566
        )
    } else if (build == "hg18") {
        chr.lens <- c(
            247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
            158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
            114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
            63811651, 62435964, 46944323, 49691432, 154913754, 57772954
        )
    } else if (build == "hg38") { # hg38
        chr.lens <- c(
            248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
            159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
            114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
            58617616, 64444167, 46709983, 50818468, 156040895, 57227415
        )
    } else {
        stop("Available reference builds: hg18, hg19, hg38")
    }

    chr.lens
}

#--- Change segment sizes into linear scale
.transformSegments <- function(segmentedData, build = "hg19") {
    Start_Position <- End_Position <- Chromosome <- NULL

    build.opts <- c("hg19", "hg18", "hg38")

    if (!build %in% build.opts) {
        stop("Available reference builds: hg18, hg19, hg38")
    }

    # Get chr lens
    chr.lens <- .getContigLens(build = build)

    segmentedData[, Start_Position := as.numeric(as.character(Start_Position))]
    segmentedData[, End_Position := as.numeric(as.character(End_Position))]

    # Replace chr x and y with numeric value (23 and 24) for better ordering
    segmentedData$Chromosome <- gsub(
        pattern = "chr", replacement = "",
        x = segmentedData$Chromosome,
        fixed = TRUE
    )
    segmentedData$Chromosome <- gsub(
        pattern = "X", replacement = "23",
        x = segmentedData$Chromosome, fixed = TRUE
    )
    segmentedData$Chromosome <- gsub(
        pattern = "Y", replacement = "24",
        x = segmentedData$Chromosome, fixed = TRUE
    )

    segmentedData$Chromosome <- factor(
        x = segmentedData$Chromosome,
        levels = seq_len(24), labels = seq_len(24)
    )

    segmentedData <- segmentedData[order(Chromosome, Start_Position,
        decreasing = FALSE
    )]

    seg.spl <- split(segmentedData, segmentedData$Chromosome)

    seg.spl.transformed <- seg.spl[[1]]
    if (nrow(seg.spl.transformed) > 0) {
        seg.spl.transformed$Start_Position_updated <-
            seg.spl.transformed$Start_Position
        seg.spl.transformed$End_Position_updated <-
            seg.spl.transformed$End_Position
    }

    chr.lens.sumsum <- cumsum(chr.lens)

    for (i in seq(2, length(seg.spl))) {
        x.seg <- seg.spl[[i]]
        if (nrow(x.seg) > 0) {
            x.seg$Start_Position_updated <- x.seg$Start_Position + chr.lens.sumsum[i - 1]
            x.seg$End_Position_updated <- x.seg$End_Position + chr.lens.sumsum[i - 1]
        }
        seg.spl.transformed <- rbind(seg.spl.transformed, x.seg, fill = TRUE)
    }

    return(seg.spl.transformed)
}

# Expected clonal VAFs for copy number CN at a given purity on autosomes
.expectedClVAF <- function(CN, purity) {
    seq_len(CN) * purity / (purity * CN + 2 * (1 - purity))
}



# Plot shared legend, taken from https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
.grid_arrange_shared_legend <- function(..., ncol = length(list(...)),
                                        nrow = 1,
                                        position = c("bottom", "right")) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(vapply(g, function(x) x$name, character(1)) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(position,
        "bottom" = gridExtra::arrangeGrob(
            do.call(gridExtra::arrangeGrob, gl),
            legend,
            ncol = 1, heights = grid::unit.c(
                unit(1, "npc") - lheight, lheight
            )
        ),
        "right" = gridExtra::arrangeGrob(
            do.call(
                gridExtra::arrangeGrob, gl
            ),
            legend,
            ncol = 2, widths = grid::unit.c(unit(1, "npc") -
                lwidth, lwidth)
        )
    )

    # return gtable invisibly
    invisible(combined)
}
