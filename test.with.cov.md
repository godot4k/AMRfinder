# dmr.with.cov test

## Test data

The `dmr.with.cov` branch was tested with the longest single chromosome from the common-site methylation matrix, not with the built-in demo data.

Input file:

```text
../data_input/merged_cln_cov/methylation_matrix_common_sites.tsv
```

Chromosome row counts showed chromosome `1` was the longest single chromosome:

```text
1 90858
10 61036
3 59896
5 59137
7 53557
```

The test subset was therefore:

```r
dat <- read.delim(
  "../data_input/merged_cln_cov/methylation_matrix_common_sites.tsv",
  check.names = FALSE
)
dat <- dat[as.character(dat[["chr"]]) == "1", ]
names(dat)[2] <- "pos"
```

Phenotype and covariate model:

```r
y <- data.frame(y = c(0, 0, 0, 0, 1, 1, 1, 1))
cov.mod <- data.frame(batch = rep(c(0, 1), 4))
```

## Method

This branch keeps the AMRfinder segmentation and k-means region summarization steps.

Region-level testing uses logistic regression:

```r
glm(y ~ x + batch, data = lm.dat, family = binomial())
```

The output columns include:

```text
cor_est coef_glm p_value
```

where `p_value` is the p value for the methylation beta term `x`.

## Command

```bash
Rscript -e ".libPaths(c('r-lib', .libPaths())); library(AMRfinder); dat <- read.delim('../data_input/merged_cln_cov/methylation_matrix_common_sites.tsv', check.names=FALSE); dat <- dat[as.character(dat[['chr']]) == '1', ]; names(dat)[2] <- 'pos'; y <- data.frame(y=c(0,0,0,0,1,1,1,1)); cov.mod <- data.frame(batch=rep(c(0,1), 4)); ctl <- list(maxdist=300, method='pearson', maxseg=-1, mincpgs=5, threads=1, mode=1, mtc=1, name='sample', trend=0, minNo=-1, minFactor=0.8, valley=0, minMethDist=0.1, randomseed=26061981); t0 <- proc.time(); res <- AMRfinder(dat, y, cov.mod, ctl); elapsed <- proc.time() - t0; print(dim(dat)); print(dim(res)); print(head(res)); print(tail(res)); print(elapsed); write.table(res, 'package-r-naive-logistic-cov-common-site-chr1-output.tsv', sep='\t', quote=FALSE, row.names=FALSE);"
```

## Result

The test completed successfully.

Input dimensions:

```text
[1] 90858    10
```

Output dimensions:

```text
[1] 5091   10
```

Runtime:

```text
user    50.157
system   1.808
elapsed 52.025 seconds
```

The retained k-means step produced three warnings:

```text
did not converge in 10 iterations
```

These warnings did not stop the run.

First output rows:

```text
chr start   end     N.CpGs cor_est      coef_glm     p_value   methX       methY FDR
1   41543   41709   5      0.47989913    19.2745354  0.2276438 0.82616064  0.5   0.8110661
1   626393  626413  5     -0.07051005    -0.3565804  0.7777305 0.21484533  0.5   0.9998751
1   626418  626447  5     -0.46300351   -11.8947168  0.2039615 0.55944617  0.5   0.8110661
1   626447  626567 13      0.32132627    19.7457185  0.3812653 0.30458655  0.5   0.8110661
```

Full local output file:

```text
package-r-naive-logistic-cov-common-site-chr1-output.tsv
```
