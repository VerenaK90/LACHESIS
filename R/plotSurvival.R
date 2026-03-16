#' Correlate SNV density timing at MRCA with Survival
#' @description
#' Takes SNV density timing as computed by `LACHESIS` as input and compares
#' survival between tumors with high and low SNV densities
#' @param lachesis output generated from \code{\link{LACHESIS}}.
#' @param mrca.cutpoint optional, MRCA density value to be used for survival
#' stratification, will be computationally inferred to maximize survival
#' differences if not specified by user.
#' @param output.dir the directory to which the plot will be stored.
#' @param surv.time column name containing survival time; defaults to `OS.time`.
#' @param surv.event column name containing event; defaults to `OS`.
#' @param surv.time.breaks numeric value controlling time axis breaks; defaults
#' to `NULL`.
#' @param surv.time.scale numeric value by which survival time is to be divided
#' (e.g., 365 for converting days into years, 30 for months), defaults to `1`.
#' @param surv.palette color palette to be used. Allowed values include "hue"
#' for the default hue color scale; "grey" for grey color palettes; brewer
#' palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g.
#' c("blue", "red").
#' @param surv.title main title.
#' @param surv.ylab y-axis label, defaults to `Survival`.
#' @param surv.xlab x-axis label, defaults to `Time`.
#' @param output.dir link to directory in which output is to be stored.
#' @return survival graphs
#' @examples
#' # An example file with sample annotations and meta data
#' input.files <- system.file("extdata", "Sample_template.txt",
#'     package =
#'         "LACHESIS"
#' )
#' input.files <- data.table::fread(input.files)
#'
#' # cnv and snv files for example tumors
#' nbe11 <- list.files(system.file("extdata/NBE11/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe15 <- list.files(system.file("extdata/NBE15/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#' nbe26 <- list.files(system.file("extdata/NBE26/", package = "LACHESIS"),
#'     full.names = TRUE
#' )
#'
#' cnv.file <- c(nbe11[1], nbe15[1], nbe26[1])
#' snv.file <- c(nbe11[2], nbe15[2], nbe26[2])
#'
#' input.files$cnv.file <- cnv.file
#' input.files$snv.file <- snv.file
#'
#' # Make an example input file with paths to cnv and snv file along with other
#' # meta data
#' lachesis_input <- tempfile(
#'     pattern = "lachesis", tmpdir = tempdir(),
#'     fileext = ".tsv"
#' )
#' data.table::fwrite(x = input.files, file = lachesis_input, sep = "\t")
#'
#' # Example with template file with paths to multiple cnv/snv files as an input
#' lachesis <- LACHESIS(input.files = lachesis_input)
#' plotSurvival(lachesis,
#'     surv.time = "EFS.time", surv.event = "EFS",
#'     mrca.cutpoint = 0.05
#' )
#'
#' @export
#' @import ggplot2
#' @import survival
#' @import survminer
#' @import gridExtra
#' @importFrom stats pchisq

plotSurvival <- function(lachesis = NULL, mrca.cutpoint = NULL,
                         output.dir = NULL, surv.time = "OS.time",
                         surv.event = "OS", surv.palette = c(
                             "dodgerblue",
                             "dodgerblue4"
                         ),
                         surv.time.breaks = NULL, surv.time.scale = 1,
                         surv.title = "Survival Probability",
                         surv.ylab = "Survival", surv.xlab = "Time") {
    MRCA_time_mean <- ECA_time_mean <- NULL

    if (is.null(lachesis)) {
        stop("'lachesis' dataset must be provided.")
    }

    if (any(is.na(lachesis[["MRCA_time_mean"]]))) {
        tmp1 <- sum(is.na(lachesis[["MRCA_time_mean"]]))
        warning(sprintf(
            "Removing %s samples with missing MRCA density estimate.",
            tmp1
        ))
        lachesis <- lachesis[!is.na(MRCA_time_mean)]
    }

    if (!surv.time %in% colnames(lachesis)) {
        stop("Please provide a valid column name for `surv.time`.")
    }

    if (!surv.event %in% colnames(lachesis)) {
        stop("Please provide a valid column name for `surv.event`.")
    }

    if (any(is.na(lachesis[[surv.time]]))) {
        tmp1 <- sum(is.na(lachesis[[surv.time]]))
        warning(sprintf(
            "Removing %s samples with missing survival time.", tmp1
        ))
        lachesis <- lachesis[!is.na(get(surv.time))]
    }

    if (any(is.na(lachesis[[surv.event]]))) {
        tmp1 <- sum(is.na(lachesis[[surv.event]]))
        warning(sprintf(
            "Removing %s samples with missing survival event.", tmp1
        ))
        lachesis <- lachesis[!is.na(get(surv.event))]
    }

    if (nrow(lachesis) == 0) {
        warning(
            "No sample with MRCA density estimate provided. Returning zero."
        )
        return(NULL)
    }

    if (all(lachesis[[surv.event]] == 0)) {
        warning("No survival events in cohort Returning zero.")
        return(NULL)
    }

    mrca.cutpoint.obj <- NULL

    # Calculating MRCA cutpoint
    if (is.null(mrca.cutpoint)) {
        mrca.cutpoint.obj <- survminer::surv_cutpoint(
            lachesis,
            time = surv.time,
            event = surv.event,
            variables = c("MRCA_time_mean")
        )

        mrca.cutpoint <- as.numeric(
            mrca.cutpoint.obj$cutpoint[
                "MRCA_time_mean",
                "cutpoint"
            ]
        )
    }

    # Categorizing according to MRCA
    lachesis.categorized <- lachesis
    lachesis.categorized[[surv.time]] <-
        as.numeric(lachesis.categorized[[surv.time]]) / surv.time.scale
    lachesis.categorized$MRCA_timing <-
        ifelse(lachesis.categorized$MRCA_time_mean < mrca.cutpoint, "Early MRCA",
            "Late MRCA"
        )
    lachesis.categorized$MRCA_timing <-
        factor(lachesis.categorized$MRCA_timing, levels = c("Early MRCA", "Late MRCA"))

    # Survival analysis
    survival.fit <- survival::survfit(
        Surv(
            time  = lachesis.categorized[[surv.time]],
            event = lachesis[[surv.event]]
        ) ~ MRCA_timing,
        data = lachesis.categorized
    )
    survival.diff <- survival::survdiff(
        Surv(
            time  = lachesis.categorized[[surv.time]],
            event = lachesis[[surv.event]]
        ) ~ MRCA_timing,
        data = lachesis.categorized
    )

    p_value <- 1 - stats::pchisq(survival.diff$chisq, length(survival.diff$n) -
        1)
    p.value.pos <- max(survival.fit$time) * (1 / 6)

    survival.fit.plot <- survminer::ggsurvplot_df(
        surv_summary(survival.fit, data = lachesis.categorized),
        title = surv.title, conf.int = TRUE, color = "strata",
        censor.shape = 124,
        palette = surv.palette, xlab = surv.xlab, ylab = surv.ylab,
        legend.labs = c("Early MRCA", "Late MRCA"),
        legend.title = " ",
        break.time.by = surv.time.breaks,
        font.main = c(15, "bold", "black"),
        font.x = 12,
        font.y = 12,
        font.tickslab = 10,
        font.legend = 10
    ) +
        annotate("text",
            x = p.value.pos, y = 0.2,
            label = paste0("p = ", ifelse(p_value > 0 &
                p_value < 0.0001, "< 0.0001",
            formatC(p_value,
                format = "f",
                digits = 4
            )
            )), size = 4
        ) +
        theme(plot.title = element_text(hjust = 0.5))

    survival.fit.risk.table <- survminer::ggrisktable(
        survival.fit,
        data = lachesis.categorized,
        colour = "strata",
        legend.labs = c("Early MRCA", "Late MRCA"),
        legend.title = NULL,
        break.time.by = surv.time.breaks,
        fontsize = 4,
        ylab = NULL,
        xlab = surv.xlab,
        ggtheme = theme_classic() +
            theme(
                axis.title.y = element_blank(),
                axis.text.y = element_text(size = 10),
                axis.text.x = element_text(size = 10)
            )
    )

    # Printing cutpoint txt
    if (!is.null(output.dir)) {
        if (!is.null(mrca.cutpoint.obj)) {
            mrca.cutpoint.dt <- data.table::data.table(
                cutpoint = mrca.cutpoint,
                statistic = as.numeric(mrca.cutpoint.obj$cutpoint[
                    "MRCA_time_mean",
                    "statistic"
                ]),
                p_value = p_value
            )

            mrca.cutpoint.rounded <- formatC(mrca.cutpoint,
                format = "f",
                digits = 2
            )
            data.table::fwrite(mrca.cutpoint.dt,
                file = file.path(output.dir, paste0(
                    "cutpoint_estimate_",
                    surv.event, "_",
                    mrca.cutpoint.rounded,
                    ".txt"
                )), sep = "\t"
            )
        }
    }

    # Printing pdf
    if (!is.null(output.dir)) {
        pdf(paste0(output.dir, "/Stratified_", surv.event, ".pdf"),
            width = 9,
            height = 8
        )
    }
    gridExtra::grid.arrange(survival.fit.plot, survival.fit.risk.table,
        ncol = 1, nrow = 2,
        widths = c(1),
        heights = c(3, 1)
    )
    if (!is.null(output.dir)) {
        dev.off()
    }
}
