#' AMRfinder: Main function for segmentation analysis
#'
#' @param intput_dat The input data matrix.
#' @param y The response variable vector.
#' @param cov.mod A data frame of covariates.
#' @param controllist A list of control parameters
#' @return A data frame with segmentation results including p-values, coefficients, etc.
#' @export
AMRfinder <- function(intput_dat, y, cov.mod = NULL, controlist = list(
                           maxdist = 300, maxseg = -1, mincpgs = 5, threads = 1, mode = 1, mtc = 1, name = "sample", trend = 0.6,
                           minNo = -1, minFactor = 0.8, valley = 0.7, minMethDist = 0.1, randomseed = 26061981
                       )) {
    nfo <- controlist
    id <- grep("sample", colnames(intput_dat))
    if (nfo$minNo < 0) nfo$minNo <- floor(length(id) * nfo$minFactor)


    nfo[["outputList"]] <- NULL
    chr.ids <- unique(intput_dat$chr)

    for (i in 1:length(chr.ids)) {
        d1 <- intput_dat[which(intput_dat$chr == chr.ids[i]), ]
        id1 <- which(diff(d1$pos) > nfo$maxdist)
        if (length(id1) > 0) {
            start1 <- c(1, id1 + 1)
            end1 <- c(id1, nrow(d1))
            for (j in 1:(length(start1))) {
                d2 <- d1[start1[j]:end1[j], ]
                if (nfo$maxseg > 0 && nrow(d2) > nfo$maxseg) {
                    al2 <- c(1:nrow(d2))
                    id2 <- which(al2 %% nfo$maxseg == 0)
                    start2 <- c(1, id2 + 1)
                    end2 <- c(id, nrow(d2))
                    for (k in 1:length(start2)) {
                        d3 <- d2[start2[k]:end2[k], ]
                        out <- segmentation(d3, y, cov.mod, chr.ids[i], nfo$mincpgs, nfo$trend, nfo$valley)
                        nfo$outputList <- rbind(nfo$outputList, as.data.frame(out))
                    }
                } else {
                    out <- segmentation(d2, y, cov.mod, chr.ids[i], nfo$mincpgs, nfo$trend, nfo$valley)
                    nfo$outputList <- rbind(nfo$outputList, as.data.frame(out))
                }
            }
        } else {
            if (nfo$maxseg > 0 && nrow(d1) > nfo$maxseg) {
                al4 <- c(1:nrow(d1))
                id4 <- which(al4 %% nfo$maxseg == 0)
                start4 <- c(1, id4 + 1)
                end4 <- c(id4, nrow(d1))
                for (k in 1:length(start4)) {
                    d4 <- d1[start4[k]:end4[k], ]
                    out <- segmentation(d4, y, cov.mod, chr.ids[i], nfo$mincpgs, nfo$trend, nfo$valley)
                    nfo$outputList <- rbind(nfo$outputList, as.data.frame(out))
                }
            } else {
                out <- segmentation(d1, y, cov.mod, chr.ids[i], nfo$mincpgs, nfo$trend, nfo$valley)
                nfo$outputList <- rbind(nfo$outputList, as.data.frame(out))
            }
        }
    }
    colnames(nfo$outputList)[-4] <- c("chr", "start", "end", "N.CpGs", "cor_est", "coef_lm", "p_value", "methX", "methY")
    nfo$outputList$FDR <- p.adjust(nfo$outputList$p_value, method = "BH")
    return(nfo$outputList[, -4])
}
