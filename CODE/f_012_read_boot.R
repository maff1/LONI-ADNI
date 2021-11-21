################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F012 READ BOOTSTRAP MODEL COEFFICIENTS
# ---------------------------------------------------------------------------- #

rm(list = ls())

pheno <-c("BPDIA", "BPSYS", "MMSE", "DX", "EXF", "LAN", "MEM")
lsPheno = lapply(pheno, 
       function(x) list.files(path = "./RESULTS/BOOT", pattern = x, full.names = TRUE)
); names(lsPheno) <- pheno

lsRes = lapply(lsPheno$BPDIA, function(x) readRDS(x))
bpdia = do.call(rbind, lsRes)
bpdia["metabolite"] <- rep(
  gsub("_BPDIA_boot_1000.rds", "", basename(lsPheno$BPDIA)),
  each = 67)
# ---------------------------------------------------------------------------- #

f_make_tbl <- function(x) {
  lsdf <- lapply(x, function(dat) readRDS(dat))
  df = do.call(rbind, lsdf)
  df["metabolite"] <- rep(
    gsub("_[A-Z]{2,5}_boot_1000.rds", "", basename(x)), 
    each = nrow(lsdf[[1]]))
  return(df)
}

lRes <- list()
for(p in 1:length(pheno)) {
  lRes[[p]] <- do.call(rbind, lapply(lsPheno[[p]], function(x) f_make_tbl(x)))
}
names(lRes) <- pheno
# ---------------------------------------------------------------------------- #

res = do.call(rbind, lRes)
resF = data.frame(subset(res,
              coef == "BPDIA" |
                coef == "BPSYS" |
                coef == "MMSE" |
                coef == "DXMCI" |
                coef == "DXALZ" |
                coef == "EXF" |
                coef == "LAN" |
                coef == "MEM"
              ), row.names = NULL
)
saveRDS(resF, "./RESULTS/bootBetas.rds")
# ---------------------------------------------------------------------------- #