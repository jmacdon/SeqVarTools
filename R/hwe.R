
.f <- function(counts) {
    nhet <- counts$nAa
    ntot <- rowSums(counts)
    afreq <- (2*counts$nAA + counts$nAa) / (2*ntot)
    exp.het <- 2*afreq*(1-afreq)*ntot
    1 - (nhet/exp.het)
}

setMethod("hwe",
          "SeqVarGDSClass",
          function(gdsobj, permute=FALSE, parallel=FALSE) {
              counts <- .countGenotypes(gdsobj, permute, parallel=parallel)
              counts2 <- counts
              names(counts2) <- c("AA","AB","BB")
              afreq <- (2*counts$nAA + counts$nAa) / (2*rowSums(counts))
              p <- HWExactStats(counts2)
              f <- .f(counts)

              ## set non-biallelic or monomorphic variants to NA
              sel <- nAlleles(gdsobj) != 2L |
                  counts$nAa + counts$naa == 0L |
                      counts$nAa + counts$nAA == 0L
              p[sel] <- NA
              f[sel] <- NA

              variant.id <- seqGetData(gdsobj, "variant.id")
              cbind(variant.id, counts, afreq, p, f)
          })
