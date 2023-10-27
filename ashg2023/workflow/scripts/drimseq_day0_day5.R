run_drimseq_day0_day5 <- function(input_path, output_path) {
    set.seed(42)
    Sys.setenv("TZDIR"=paste0(Sys.getenv("CONDA_PREFIX"), "/share/zoneinfo"))
    #print("YUPSTER?")
    #library(edgeR)
    library(DRIMSeq)
    #print("EVEN GETTTIUNG HERE?")
    library(vroom)
    #print("NOT HERE THO")
    library(dplyr)
    #print("OR MAYBE?")

    counts <- data.frame(vroom::vroom(input_path), check.names = FALSE)

    counts <- (counts[, c(1:5, 11:13)])
    colnames(counts) <- c("feature_id", "gene_id", "day0-rep1", "day0-rep2", "day0-rep3", "day5-rep1", "day5-rep2", "day5-rep3")

    md <- data.frame(
    sample_id = c("day0-rep1", "day0-rep2", "day0-rep3", "day5-rep1", "day5-rep2", "day5-rep3"),
    group = c(rep("day0", 3), rep("day5", 3))
    
    )

    d <- dmDSdata(counts = counts, samples = md)

    d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_expr = 3,
                min_gene_expr = 100, min_feature_expr = 20)
    design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
    design_full
    d <- dmPrecision(d, design = design_full)

    d <- dmFit(d, design = design_full, verbose = 1)

    d <- dmTest(d, coef = "groupday5", verbose = 1)

    data.frame(results(d, level = "gene"), check.names = FALSE) %>% vroom::vroom_write(
        paste0(output_path, "_gene_level_results.tsv")
    )

    data.frame(results(d, level = "feature"), check.names=FALSE) %>% vroom::vroom_write(
        paste0(output_path, "_transcript_level_results.tsv")
    )



}

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("--input_path")
parser$add_argument("--output_path")
args <- parser$parse_args()
run_drimseq_day0_day5(input_path=args$input_path,
                    output_path=args$output_path
)