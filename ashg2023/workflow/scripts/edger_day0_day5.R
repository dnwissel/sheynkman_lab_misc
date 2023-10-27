run_edger_day0_day5 <- function(input_path, output_path) {
    set.seed(42)
    Sys.setenv("TZDIR"=paste0(Sys.getenv("CONDA_PREFIX"), "/share/zoneinfo"))
    library(edgeR)
    library(vroom)
    library(dplyr)
    
    samples <- c("day0-rep1", "day0-rep2", "day0-rep3", "day5-rep1", "day5-rep2", "day5-rep3")
    counts <- data.frame(vroom::vroom(input_path), check.names=FALSE)
    if (ncol(counts) == 13) {
        counts <- counts[, -2]
    }

    dge <- edgeR::DGEList(counts=counts[, unname(sapply(samples, function(x) grep(x, colnames(counts))))], group = c(paste0("day0-rep", 1:3), paste0("day5-rep", 1:3)), genes = counts[, 1])

    keep <- edgeR::filterByExpr(dge)
    dge <- dge[keep, ]
    cps <- edgeR::cpm(dge)
    rs <- rowSums(cps > 1) >= 3

    dge <- dge[rs, ]

    ge <- calcNormFactors(dge)
    ss <- sapply(strsplit(colnames(dge), ".", fixed=TRUE), function(x) x[[1]])

    day <- sapply(ss, function(u) ifelse(grepl("2023",u[1]), u[2], u[1]))

    colnames(dge) <- day

    (grp <- sapply(strsplit(day, "-", fixed=TRUE), .subset, 1))
    (mm <- model.matrix(~0+grp))

    cps <- edgeR::cpm(dge)
    cps <- cps[,order(colnames(cps))]
    dge$genes <- data.frame(dge$genes, round(cps, 2))

    head(dge$genes)

    # estimate dispersion
    dge <- estimateDisp(dge,mm)

    f <- glmQLFit(dge,mm)
    mc <- makeContrasts(grpday0-grpday5, levels=colnames(mm))
    qlf <- glmQLFTest(f,contrast = mc)


    tt <- topTags(qlf, n = Inf, sort.by = "none")
    data.frame(tt, check.names=FALSE) %>% vroom::vroom_write(output_path)

}

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("--input_path")
parser$add_argument("--output_path")
args <- parser$parse_args()
run_edger_day0_day5(input_path=args$input_path,
                    output_path=args$output_path
)