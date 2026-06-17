findMaxZ<-function(s,t,mincpgs,XS){
  max<-0
  a_max_id<-0
  b_max_id<-1
  S <-function(s,t,XS){
    dat<-as.numeric(XS[,2])
    Sst<-dat[t]
    if (s>1) Sst <- Sst - dat[s=1]
    Sst <- abs(Sst)
    return(Sst)
  }
  if(s<=t){
    for(a in s:t){
      if(a+mincpgs-1<=t){
        for(b in (a+mincpgs-1):t){
          if ((b-a-mincpgs == 0) || ((a ==s) && (b==t)) || (a==b)){
            Zab <- 0
          }else{
            lab <- b - a + 1
            lst <- t - s + 1
            u <- (S(a, b, XS) - (lab * S(s, t, XS)/lst))^2
            Zab <- u / (lab * (1-lab/lst))
          }
          if((is.finite(Zab)==TRUE)&&((a!=s)||(b!=t))&&(max<Zab)) {
            max<-Zab
            a_max_id<-(a-s+1)
            b_max_id<-(b-s+1)
          }
        }
      }
    }
  }
  max_id<-c(a_max_id,b_max_id)
  return(max_id)
}

calcSingleTrendAbs<-function(S,s,t){
  if(s==1){
    trend<-S[t,3]/(t-1)
    trend<-abs(trend)
  }
  else{
    trend<-(S[t,3]-S[s-1,3])/(t-s+1)
    trend<-abs(trend)
  }
  return (trend)
}

noValley<-function(S,s,t,mincpgs,valley){
  if((t-s+1)<mincpgs) {return (1)}
  Sst<-S[t,2]
  minlength<-ifelse(mincpgs>10,mincpgs,10)
  if(s>1){Sst<-Sst-S[s-1,2]}
  mean<-abs(Sst/(t-s+1))
  S0<-as.numeric(S[,2])
  len1<-nrow(S)+1
  S1<-append(S0,rep(0,minlength),0)
  len2<-length(S1)
  S2<-S1[-c(len1:len2)]
  imean<-abs((S0-S2)/minlength)
  imean<-imean[(s+minlength-1):t]
  y<-which(abs(imean)<mean*valley)
  if(length(y)>0){return (0)}
  return (1)
}

is_contiguous <- function(x) {
  rle_x <- rle(x)
  length(unique(rle_x$values)) == length(rle_x$values)
}

buildPhenotypeModelData <- function(meth, y, cov.mod = NULL, NR = 1) {
  y_df <- as.data.frame(y)
  phenotype <- as.numeric(y_df[, 1])
  if (!all(is.finite(phenotype))) {
    stop("For phenotype model mode, the first column of y must be a numeric phenotype.")
  }
  if (length(meth) != length(phenotype) * NR) {
    stop("The methylation summary length must equal nrow(y) * NR.")
  }

  adjustment <- NULL
  if (ncol(y_df) > 1) {
    adjustment <- y_df[, -1, drop = FALSE]
  }
  if (!is.null(cov.mod)) {
    cov_df <- as.data.frame(cov.mod)
    if (nrow(cov_df) != length(phenotype)) {
      stop("cov.mod must have one row per sample.")
    }
    adjustment <- cbind(adjustment, cov_df)
  }

  model_dat <- data.frame(
    phenotype = rep(phenotype, each = NR),
    meth = meth
  )
  if (!is.null(adjustment) && ncol(adjustment) > 0) {
    names(adjustment) <- make.unique(make.names(names(adjustment)))
    keep <- vapply(adjustment, function(z) length(unique(z[!is.na(z)])) > 1, logical(1))
    adjustment <- adjustment[, keep, drop = FALSE]
    if (ncol(adjustment) > 0) {
      model_dat <- cbind(model_dat, adjustment[rep(seq_len(nrow(adjustment)), each = NR), , drop = FALSE])
    }
  }
  char_id <- vapply(model_dat, is.character, logical(1))
  model_dat[char_id] <- lapply(model_dat[char_id], factor)
  model_dat
}

selectRandomEffect <- function(model_dat) {
  candidates <- c(
    "subject", "subject_id", "Subject", "SubjectID", "subject.ID",
    "id", "ID", "patient", "patient_id", "individual", "individual_id"
  )
  candidates <- make.names(candidates)
  candidates <- candidates[candidates %in% names(model_dat)]
  for (candidate in candidates) {
    value <- model_dat[[candidate]]
    n_level <- length(unique(value[!is.na(value)]))
    if (n_level >= 2 && any(table(value) > 1)) {
      return(candidate)
    }
  }
  NULL
}

extractMethCoefficient <- function(fit) {
  coef_table <- tryCatch(summary(fit)$coef, error = function(e) NULL)
  if (is.null(coef_table) || !"meth" %in% rownames(coef_table)) {
    return(c(p_value = 1, coef_meth = 0))
  }
  p_col <- grep("^Pr\\(", colnames(coef_table), value = TRUE)
  if (length(p_col) == 0) {
    return(c(p_value = 1, coef_meth = unname(coef_table["meth", 1])))
  }
  p_value <- unname(coef_table["meth", p_col[1]])
  coef_meth <- unname(coef_table["meth", 1])
  if (!is.finite(p_value)) p_value <- 1
  if (!is.finite(coef_meth)) coef_meth <- 0
  c(p_value = p_value, coef_meth = coef_meth)
}

fitPhenotypeMethModel <- function(meth, y, cov.mod = NULL, NR = 1, method = "pearson") {
  model_dat <- buildPhenotypeModelData(meth, y, cov.mod, NR)
  model_dat <- model_dat[complete.cases(model_dat), , drop = FALSE]
  if (nrow(model_dat) < 3 ||
      length(unique(model_dat$phenotype)) < 2 ||
      length(unique(model_dat$meth)) < 2) {
    return(c(p_value = 1, coef_meth = 0, cor_est = 0))
  }

  random_effect <- selectRandomEffect(model_dat)
  fixed_terms <- setdiff(names(model_dat), c("phenotype", "meth", random_effect))
  rhs <- c("meth", fixed_terms)
  fixed_formula <- as.formula(paste("phenotype ~", paste(rhs, collapse = " + ")))

  fit <- NULL
  if (!is.null(random_effect) && requireNamespace("lmerTest", quietly = TRUE)) {
    mixed_formula <- as.formula(
      paste("phenotype ~", paste(rhs, collapse = " + "), "+ (1 |", random_effect, ")")
    )
    model_dat[[random_effect]] <- factor(model_dat[[random_effect]])
    fit <- tryCatch(
      suppressWarnings(lmerTest::lmer(mixed_formula, data = model_dat, REML = FALSE)),
      error = function(e) NULL
    )
  }

  if (is.null(fit)) {
    fit <- tryCatch(
      suppressWarnings(lm(fixed_formula, data = model_dat)),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) {
    return(c(p_value = 1, coef_meth = 0, cor_est = 0))
  }

  meth_coef <- extractMethCoefficient(fit)
  cor_dat <- data.frame(y = model_dat$phenotype, x = model_dat$meth)
  cor_est <- cor(cor_dat$y, cor_dat$x, use = "complete.obs", method = method)
  if (!is.finite(cor_est)) cor_est <- 0
  p_value <- unname(meth_coef["p_value"])
  coef_meth <- unname(meth_coef["coef_meth"])
  if (!is.finite(p_value)) p_value <- 1
  if (!is.finite(coef_meth)) coef_meth <- 0
  c(p_value = p_value, coef_meth = coef_meth, cor_est = cor_est)
}

cortest <- function(intput_dat, y, method = "pearson", cov.mod = NULL, a, b) {
  set.seed(123)
  aa <- intput_dat[a:b, ]
  if (nrow(aa) <= 1) {
    op.num <- 1
  } else {
    cls.num <- NULL
    unique_points <- nrow(aa[, -c(1, 2)])
    max_possible_k <- min(unique_points - 1, 4)
    possible_ks <- 1:max_possible_k
    possible_ks <- possible_ks[possible_ks > 0 & possible_ks <= unique_points]
    for (k in possible_ks) {
      cls <- tryCatch(
        {
          tmp_data <- aa[, -c(1, 2)]
          tmp_data <- t(apply(tmp_data, 1, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); if(all(is.na(x))) x <- rep(0, length(x)); x }))
          kmeans(tmp_data, centers = k)
        },
        error = function(e) NULL
      )
      if (!is.null(cls) && is_contiguous(as.numeric(cls$cluster))) {
        cls.num <- c(cls.num, k)
      }
    }
    op.num <- ifelse(length(cls.num) > 0, max(cls.num), 1)
  }
  if (op.num > 1) {
    tmp_data_final <- aa[, -c(1, 2)]
    tmp_data_final <- t(apply(tmp_data_final, 1, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); if(all(is.na(x))) x <- rep(0, length(x)); x }))
    cls.op <- kmeans(tmp_data_final, centers = op.num)
    x.mean <- cls.op$centers
    NR <- nrow(x.mean)
    x <- as.numeric(unlist(x.mean))
  } else {
    x <- as.numeric(colMeans(aa[, -c(1, 2)], na.rm = TRUE))
    NR <- 1
  }
  model_stat <- fitPhenotypeMethModel(x, y, cov.mod, NR, method)
  model_stat[!is.finite(model_stat)] <- c(p_value = 1, coef_meth = 0, cor_est = 0)[!is.finite(model_stat)]
  return(c(unname(model_stat["p_value"]), unname(model_stat["coef_meth"]), unname(model_stat["cor_est"])))
}

calcSingleDiffSum<-function(intput_dat,y){
  phenotype <- as.numeric(as.data.frame(y)[, 1])
  if (!all(is.finite(phenotype)) || length(unique(phenotype)) < 2) {
    stop("The first column of y must be a numeric phenotype with at least two values.")
  }
  calcCR<-function(x){
    x<-as.numeric(x)
    res <- cor(x, phenotype, method = "pearson", use = "complete.obs")
    if(is.na(res)) res <- 0
    return(res)
  }
  correlation<-apply(intput_dat[,-c(1,2)],1,calcCR)
  smean<-cumsum(as.numeric(correlation))
  absmean<-abs(smean)
  sigm<-sign(correlation)
  sigsum<-cumsum(sigm)
  S<-cbind(absmean,smean,sigsum)
  return(S)
}

segment_pSTKopt<-function(intput_dat,y,cov.mod,XS,a,b,chr,mincpgs,trend,valley,KS,method){
  stacks<-NULL
  breaks<-NULL
  child<-0
  ab<-c(-1,0)
  ks1<-c(2,2)
  ks2<-c(2,2)
  ks3<-c(2,2)
  while(length(stacks)||(a!=-1)){
    if((a!=-1)&&(child<=2)){
      if(ab[1]==-1){
        ab<-c(0,0)
        ks1<-c(2,2)
        ks2<-c(2,2)
        ks3<-c(2,2)
        max_id<-findMaxZ(a,b,mincpgs,XS)
        ab<-(a-1)+max_id
        n<-a
        m<-ab[1]-1
        if((ab[1]>1)&&(m-n+1>=mincpgs)&&calcSingleTrendAbs(XS,n,m)>trend&&noValley(XS,n,m,mincpgs,valley)){
          ks1<-cortest(intput_dat,y, method, cov.mod,n,m)
        }
        n<-ab[1]
        m<-ab[2]
        if(m-n+1>=mincpgs&&calcSingleTrendAbs(XS,n,m)>trend&&noValley(XS,n,m,mincpgs,valley)){
          ks2<-cortest(intput_dat,y, method, cov.mod,n,m)
        }
        n<-ab[2]+1
        m<-b
        if(m-n+1>=mincpgs&&calcSingleTrendAbs(XS,n,m)>trend&&noValley(XS,n,m,mincpgs,valley)){
          ks3<-cortest(intput_dat,y,method, cov.mod,n,m)
        }
      }
      newp<-min(ks1[1],ks2[1],ks3[1], na.rm = TRUE)
      if(!is.finite(newp)) newp <- 2
      if((newp<KS[1]||(newp>1&&KS[1]>1))&&(b-a>=mincpgs)){
        stack<-list(a=a,b=b,ab=ab,child=child+1,ks1=ks1,ks2=ks2,ks3=ks3,KS=KS)
        stacks<-append(stacks,list(stack))
        if(child==0){
          n<-a
          m<-ab[1]-1
          a<--1
          if(ab[1]>1&&n<=m){
            if(m-n>=mincpgs){
              a<-n
              b<-m
              KS<-ks1
              child<-0
              ab[1]<--1
            }else{
              bre<-list(chr=chr,start=n,stop=m,p_value=ks1[1],coef_meth=ks1[2],cor_est=ks1[3])
              breaks<-rbind(breaks,data.frame(bre))
            }
          }
        }
        if(child==1){
          n<-ab[1]
          m<-ab[2]
          a<--1
          if(n<=m){
            if(m-n>=mincpgs){
              a<-n
              b<-m
              KS<-ks2
              child<-0
              ab[1]<--1
            }else{
              bre<-list(chr=chr,start=n,stop=m,p_value=ks2[1],coef_meth=ks2[2],cor_est=ks2[3])
              breaks<-rbind(breaks,data.frame(bre))
            }
          }
        }
        if(child==2){
          n<-ab[2]+1
          m<-b
          a<--1
          if(n<=m){
            if(m-n>=mincpgs){
              a<-n
              b<-m
              KS<-ks3
              child<-0
              ab[1]<--1
            }else{
              bre<-list(chr=chr,start=n,stop=m,p_value=ks3[1],coef_meth=ks3[2],cor_est=ks3[3])
              breaks<-rbind(breaks,data.frame(bre))
            }
          }
        }
      }else{
        if(child==0){
          bre<-list(chr=chr,start=a,stop=b,p_value=KS[1],coef_meth=KS[2],cor_est=KS[3])
          breaks<-rbind(breaks,data.frame(bre))
        }
        a<--1
        ab[1]<--1
      }
    }else{
      len<-length(stacks)
      pop<-stacks[[len]]
      stacks[[length(stacks)]]<-NULL
      a<-pop[[1]]
      b<-pop[[2]]
      ab<-pop[[3]]
      child<-pop[[4]]
      ks1<-pop[[5]]
      ks2<-pop[[6]]
      ks3<-pop[[7]]
      KS<-pop[[8]]
      if(child==3){
        a<--1
      }
    }
  }
  return (breaks)
}

segmenterSTK<-function(intput_dat,y,cov.mod,XS,a,b,chr,mincpgs,trend,valley,KS,method){
  stacks<-NULL
  globalbreaks<-NULL
  while(length(stacks)>0||a!=-1){
    if(a!=-1){
      bre<-segment_pSTKopt(intput_dat,y,cov.mod,XS,a,b,chr,mincpgs,trend,valley,KS,method)
      i<-which(bre$p_value==min(bre$p_value))
      max<-bre[i,]
      max<-max[1,]
      stack<-list(a=a,b=b,max=max)
      stacks<-append(stacks,list(stack))
      n<-a
      m<-max$start-1
      if(max$start>1&&n<=m){
        ks<-c(2,2)
        if(m-n+1>=mincpgs&&calcSingleTrendAbs(XS,n,m)>trend&&noValley(XS,n,m,mincpgs,valley)){
          ks<-cortest(intput_dat,y, method, cov.mod,n,m)
        }
        a<-n
        b<-m
        KS<-ks
      }else{
        a<--1
      }
    }
    else{
      pop<-stacks[[length(stacks)]]
      stacks[[length(stacks)]]<-NULL
      a<-pop$a
      b<-pop$b
      max<-pop$max
      globalbreaks<-rbind(globalbreaks,as.data.frame(max))
      n<-max$stop+1
      m<-b
      if(n<=m){
        ks<-c(2,2)
        if(m-n+1>=mincpgs&&calcSingleTrendAbs(XS,n,m)>trend&&noValley(XS,n,m,mincpgs,valley)){
          ks<-cortest(intput_dat,y, method, cov.mod,n,m)
        }
        a<-n
        b<-m
        KS<-ks
      }else{
        a<--1
      }
    }
  }
  return(globalbreaks)
}

output<-function(intput_dat,y,cov.mod,XS,global,chr,mincpgs,trend,valley,method){
  tmp<-NULL
  outputList<-NULL
  collapseTmp<-function(tmp){
    if(!is.null(tmp)){
      tmp<-as.data.frame(tmp)
      start<-min(as.numeric(unlist(tmp$start)))
      stop<-max(as.numeric(unlist(tmp$stop)))
      tmp<-tmp[1,]
      tmp$start<-start
      tmp$stop<-stop
    }
    return(tmp)
  }
  nbreaks<-nrow(global)
  if(nbreaks>0){
    for(i in 1:nbreaks){
      b<-global[i,]
      if(b$p_value>1){
        if(is.null(tmp)){
          tmp<-b
          methX <- mean(as.numeric(as.matrix(intput_dat[tmp$start:tmp$stop, -c(1, 2)])), na.rm = TRUE)
          methY <- mean(as.numeric(as.data.frame(y)[, 1]), na.rm = TRUE)
          tmp$methX<-methX
          tmp$methY<-methY
        }else{
          tmp$stop<-b$stop
        }
      }else{
        if(!is.null(tmp)){
          tmp<-collapseTmp(tmp)
          tmp_start<-as.integer(as.numeric(unlist(tmp$start))[1])
          tmp_stop<-as.integer(as.numeric(unlist(tmp$stop))[1])
          ks<-c(2,2)
          if(isTRUE((tmp_stop-tmp_start+1>=mincpgs)[1])&&isTRUE((calcSingleTrendAbs(XS,tmp_start,tmp_stop)>trend)[1])&&isTRUE((noValley(XS,tmp_start,tmp_stop,mincpgs,valley)==1)[1])){
            ks<-cortest(intput_dat,y, method, cov.mod,tmp_start,tmp_stop)
          }
          if(ks[1]<2){
            out<-data.frame(chr=chr,start=intput_dat$pos[tmp_start]-1,stop=intput_dat$pos[tmp_stop],q=-1,length=tmp_stop-tmp_start+1,cor_est = ks[3],coef_meth=ks[2],p_value=ks[1])
            methX <- mean(as.numeric(as.matrix(intput_dat[tmp_start:tmp_stop, -c(1, 2)])), na.rm = TRUE)
            methY <- mean(as.numeric(as.data.frame(y)[, 1]), na.rm = TRUE)
            out$methX<-methX
            out$methY<-methY
            outputList<-rbind(outputList,as.data.frame(out))
          }
          tmp<-NULL
        }
        out<-data.frame(chr=chr,start=intput_dat$pos[b$start]-1,stop=intput_dat$pos[b$stop],q=-1,length=b$stop-b$start+1,cor_est =b$cor_est,coef_meth=b$coef_meth,p_value=b$p_value)
        methX <- mean(as.numeric(as.matrix(intput_dat[b$start:b$stop, -c(1, 2)])), na.rm = TRUE)
        methY <- mean(as.numeric(as.data.frame(y)[, 1]), na.rm = TRUE)
        out$methX<-methX
        out$methY<-methY
        outputList<-rbind(outputList,as.data.frame(out))
      }
    }
  }
  if(!is.null(tmp)){
    tmp<-collapseTmp(tmp)
    tmp_start<-as.integer(as.numeric(unlist(tmp$start))[1])
    tmp_stop<-as.integer(as.numeric(unlist(tmp$stop))[1])
    ks<-c(2,2)
    if(isTRUE((tmp_stop-tmp_start+1>=mincpgs)[1])&&isTRUE((calcSingleTrendAbs(XS,tmp_start,tmp_stop)>trend)[1])&&isTRUE((noValley(XS,tmp_start,tmp_stop,mincpgs,valley)==1)[1])){
      ks<-cortest(intput_dat,y, method, cov.mod,tmp_start,tmp_stop)
    }
    if(ks[1]<2){
      out<-data.frame(chr=chr,start=intput_dat$pos[tmp_start]-1,stop=intput_dat$pos[tmp_stop],q=-1,length=tmp_stop-tmp_start+1,cor_est = ks[3],coef_meth=ks[2],p_value=ks[1])
      methX <- mean(as.numeric(as.matrix(intput_dat[tmp_start:tmp_stop, -c(1, 2)])), na.rm = TRUE)
      methY <- mean(as.numeric(as.data.frame(y)[, 1]), na.rm = TRUE)
      out$methX<-methX
      out$methY<-methY
      outputList<-rbind(outputList,as.data.frame(out))
    }
    tmp<-NULL
  }
  return(outputList)
}

segmentation<-function(intput_dat,y,cov.mod,chr,mincpgs,trend,valley,method){
  ks<-c(2,2)
  len<-nrow(intput_dat)
  XS<-calcSingleDiffSum(intput_dat,y)
  if(len-1>=mincpgs&&calcSingleTrendAbs(XS,1,len)>trend&&noValley(XS,1,len,mincpgs,valley)){
    ks<-cortest(intput_dat,y, method, cov.mod,1,len)
  }
  globalbreaks<-segmenterSTK(intput_dat,y,cov.mod,XS,1,len,chr,mincpgs,trend,valley,ks,method)
  outputList<-output(intput_dat,y,cov.mod,XS,globalbreaks,chr,mincpgs,trend,valley,method)
  return (outputList)
}
