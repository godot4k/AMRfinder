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

metilene_kscdf <- function(x) {
  if (!is.finite(x) || x < 0) return(1)
  sum_value <- 0
  old <- 0
  coeff <- 1
  base <- -2 * x * x
  for (k in seq_len(100)) {
    tmp <- exp(base * k * k)
    sum_value <- sum_value + coeff * tmp
    if (tmp <= 1e-3 * old || tmp <= 1e-8 * sum_value) {
      return(max(0, min(1, 2 * sum_value)))
    }
    coeff <- -coeff
    old <- tmp
  }
  1
}

metilene_calc_max <- function(l1, l2, c1, c2, a, b, m, n) {
  l1a <- l1[a + 1]
  l2a <- l2[a + 1]
  l1b <- l1[b + 1]
  l2b <- l2[b + 1]
  max(
    abs((c1 / m) - (c2 / n)),
    abs(((c1 + l1[1]) / m) - ((c2 + l2[1]) / n)),
    abs(((c1 + l1a) / m) - ((c2 + l2a) / n)),
    abs(((c1 + l1b) / m) - ((c2 + l2b) / n)),
    abs(((c1 + l1[1] + l1a) / m) - ((c2 + l2[1] + l2a) / n)),
    abs(((c1 + l1[1] + l1b) / m) - ((c2 + l2[1] + l2b) / n)),
    abs(((c1 + l1b + l1a) / m) - ((c2 + l2b + l2a) / n)),
    abs(((c1 + l1b + l1a + l1[1]) / m) - ((c2 + l2b + l2a + l2[1]) / n))
  )
}

metilene_counter <- function(x, y, x0, x1, y0, y1) {
  m <- length(x0)
  n <- length(y0)
  l_control <- c(
    sum(x0 == x & x1 == y),
    sum(x0 == x & x1 < y),
    sum(x0 == x & x1 > y),
    sum(x0 < x & x1 == y),
    sum(x0 > x & x1 == y)
  )
  l_test <- c(
    sum(y0 == x & y1 == y),
    sum(y0 == x & y1 < y),
    sum(y0 == x & y1 > y),
    sum(y0 < x & y1 == y),
    sum(y0 > x & y1 == y)
  )
  c_control <- c(
    sum(x0 > x & x1 > y),
    sum(x0 > x & x1 < y),
    sum(x0 < x & x1 > y),
    sum(x0 < x & x1 < y)
  )
  c_test <- c(
    sum(y0 > x & y1 > y),
    sum(y0 > x & y1 < y),
    sum(y0 < x & y1 > y),
    sum(y0 < x & y1 < y)
  )
  max(
    metilene_calc_max(l_control, l_test, c_control[1], c_test[1], 2, 4, m, n),
    metilene_calc_max(l_control, l_test, c_control[2], c_test[2], 1, 4, m, n),
    metilene_calc_max(l_control, l_test, c_control[3], c_test[3], 2, 3, m, n),
    metilene_calc_max(l_control, l_test, c_control[4], c_test[4], 1, 3, m, n)
  )
}

safe_cor <- function(x, y) {
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) return(0)
  value <- suppressWarnings(cor(x, y))
  ifelse(is.finite(value), value, 0)
}

metilene_ks2d <- function(control_values, control_pos, test_values, test_pos) {
  m <- length(control_values)
  n <- length(test_values)
  if (m == 0 || n == 0) return(c(p_value = 1, statistic = 0))
  d_control <- 0
  for (j in seq_len(m)) {
    d_control <- max(
      d_control,
      metilene_counter(control_values[j], control_pos[j], control_values, control_pos, test_values, test_pos)
    )
  }
  d_test <- 0
  for (j in seq_len(n)) {
    d_test <- max(
      d_test,
      metilene_counter(test_values[j], test_pos[j], control_values, control_pos, test_values, test_pos)
    )
  }
  d_stat <- (d_control + d_test) * 0.5
  s <- sqrt(m * n / (m + n))
  cor_control <- safe_cor(control_values, control_pos)
  cor_test <- safe_cor(test_values, test_pos)
  denom <- 1 + sqrt(max(0, 1 - 0.5 * (cor_control * cor_control + cor_test * cor_test))) * (0.25 - 0.75 / s)
  if (!is.finite(denom) || denom == 0) return(c(p_value = 1, statistic = d_stat))
  p_value <- metilene_kscdf(d_stat * s / denom)
  if (!is.finite(p_value)) p_value <- 1
  c(p_value = p_value, statistic = d_stat)
}

cortest <- function(intput_dat, y, method = "pearson", cov.mod = NULL, a, b) {
  aa <- intput_dat[a:b, ]
  y_group <- as.numeric(y[, 1])
  if (!all(y_group %in% c(0, 1))) {
    stop("For control/test mode, y must be coded as 0 for control and 1 for test.")
  }
  methylation <- aa[, -c(1, 2), drop = FALSE]
  if (ncol(methylation) != length(y_group)) {
    stop("The number of methylation sample columns must match the length of y.")
  }
  control_id <- which(y_group == 0)
  test_id <- which(y_group == 1)
  if (length(control_id) == 0 || length(test_id) == 0) {
    stop("Both control (0) and test (1) samples are required for 2D-KS/MWU testing.")
  }

  pos <- as.numeric(aa$pos)
  control_values <- as.numeric(t(as.matrix(methylation[, control_id, drop = FALSE])))
  test_values <- as.numeric(t(as.matrix(methylation[, test_id, drop = FALSE])))
  control_pos <- rep(pos, each = length(control_id))
  test_pos <- rep(pos, each = length(test_id))
  control_ok <- is.finite(control_values) & is.finite(control_pos)
  test_ok <- is.finite(test_values) & is.finite(test_pos)
  control_values <- control_values[control_ok]
  control_pos <- control_pos[control_ok]
  test_values <- test_values[test_ok]
  test_pos <- test_pos[test_ok]
  if (length(control_values) == 0 || length(test_values) == 0) {
    return(c(1, 0, 0, 1))
  }

  ks2d <- metilene_ks2d(control_values, control_pos, test_values, test_pos)
  ks2d_p <- as.numeric(ks2d[["p_value"]])
  if (!is.finite(ks2d_p)) ks2d_p <- 1
  ks2d_stat <- as.numeric(ks2d[["statistic"]])
  if (!is.finite(ks2d_stat)) ks2d_stat <- 0
  mwu_p <- suppressWarnings(wilcox.test(control_values, test_values, exact = FALSE)$p.value)
  if (!is.finite(mwu_p)) mwu_p <- 1
  mean_diff <- mean(test_values, na.rm = TRUE) - mean(control_values, na.rm = TRUE)
  if (!is.finite(mean_diff)) mean_diff <- 0
  return(c(ks2d_p, mean_diff, ks2d_stat, mwu_p))
}

calcEValue <- function(intput_dat, y, a, b) {
  y_group <- as.numeric(y[, 1])
  if (!all(y_group %in% c(0, 1))) {
    stop("For e-value calculation, y must be coded as 0 for control and 1 for test.")
  }
  if (b < a) return(1)
  region_dat <- intput_dat[a:b, -c(1, 2), drop = FALSE]
  if (nrow(region_dat) == 0) return(1)
  if (ncol(region_dat) != length(y_group)) {
    stop("The number of methylation sample columns must match the length of y.")
  }

  control_id <- which(y_group == 0)
  test_id <- which(y_group == 1)
  if (length(control_id) == 0 || length(test_id) == 0) {
    stop("Both control (0) and test (1) samples are required for e-value calculation.")
  }

  density_log <- function(x, mu, sigma) {
    vector_temp <- na.omit(as.numeric(x))
    n <- length(vector_temp)
    if (n == 0 || !is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
      return(NA_real_)
    }
    value <- mean(vector_temp)
    if (!is.finite(value)) return(NA_real_)
    d <- dnorm(x = value, mean = mu, sd = sigma / sqrt(n), log = TRUE)
    if (!is.finite(d)) return(NA_real_)
    d
  }

  control_values <- unlist(region_dat[, control_id, drop = FALSE])
  test_values <- unlist(region_dat[, test_id, drop = FALSE])
  all_values <- c(control_values, test_values)
  mu_control <- mean(control_values, na.rm = TRUE)
  sigma_control <- sd(control_values, na.rm = TRUE)
  mu_test <- mean(test_values, na.rm = TRUE)
  sigma_test <- sd(test_values, na.rm = TRUE)
  mu_pooled <- mean(all_values, na.rm = TRUE)
  sigma_pooled <- sd(all_values, na.rm = TRUE)

  log_up_control <- apply(region_dat[, control_id, drop = FALSE], 2, density_log, mu = mu_control, sigma = sigma_control)
  log_down_control <- apply(region_dat[, control_id, drop = FALSE], 2, density_log, mu = mu_pooled, sigma = sigma_pooled)
  log_up_test <- apply(region_dat[, test_id, drop = FALSE], 2, density_log, mu = mu_test, sigma = sigma_test)
  log_down_test <- apply(region_dat[, test_id, drop = FALSE], 2, density_log, mu = mu_pooled, sigma = sigma_pooled)

  log_e_value <- sum(c(log_up_control, log_up_test), na.rm = TRUE) -
    sum(c(log_down_control, log_down_test), na.rm = TRUE)
  if (!is.finite(log_e_value) || log_e_value <= 0) return(1)
  if (log_e_value >= log(.Machine$double.xmax)) return(Inf)
  exp(log_e_value)
}

adjustEValueBH <- function(e_value) {
  p_value <- rep(1, length(e_value))
  finite_id <- is.finite(e_value) & e_value > 0
  p_value[finite_id] <- 1 / e_value[finite_id]
  p_value[is.infinite(e_value) & e_value > 0] <- 0
  adjusted_p <- p.adjust(p_value, method = "BH")
  e_adjust <- rep(1, length(adjusted_p))
  zero_id <- adjusted_p == 0
  positive_id <- adjusted_p > 0
  e_adjust[zero_id] <- Inf
  e_adjust[positive_id] <- 1 / adjusted_p[positive_id]
  e_adjust[!is.finite(e_adjust) & !zero_id] <- 1
  e_adjust
}

calcSingleDiffSum<-function(intput_dat,y){
  y_group <- as.numeric(y[, 1])
  if (!all(y_group %in% c(0, 1))) {
    stop("For control/test mode, y must be coded as 0 for control and 1 for test.")
  }
  control_id <- which(y_group == 0)
  test_id <- which(y_group == 1)
  if (length(control_id) == 0 || length(test_id) == 0) {
    stop("Both control (0) and test (1) samples are required.")
  }
  calcDiff<-function(x){
    x<-as.numeric(x)
    res <- mean(x[test_id], na.rm = TRUE) - mean(x[control_id], na.rm = TRUE)
    if(is.na(res)) res <- 0
    return(res)
  }
  mean_difference<-apply(intput_dat[,-c(1,2)],1,calcDiff)
  smean<-cumsum(as.numeric(mean_difference))
  absmean<-abs(smean)
  sigm<-sign(mean_difference)
  sigsum<-cumsum(sigm)
  S<-cbind(absmean,smean,sigsum)
  return(S)
}

segment_pSTKopt<-function(intput_dat,y,cov.mod,XS,a,b,chr,mincpgs,trend,valley,KS,method){
  stacks<-NULL
  breaks<-NULL
  child<-0
  ab<-c(-1,0)
  ks1<-c(2,0,0,2)
  ks2<-c(2,0,0,2)
  ks3<-c(2,0,0,2)
  while(length(stacks)||(a!=-1)){
    if((a!=-1)&&(child<=2)){
      if(ab[1]==-1){
        ab<-c(0,0)
        ks1<-c(2,0,0,2)
        ks2<-c(2,0,0,2)
        ks3<-c(2,0,0,2)
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
      newp<-min(ks1[1],ks2[1],ks3[1])
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
              bre<-list(chr=chr,start=n,stop=m,segment_p=ks1[1],p_value=ks1[4],mean_diff=ks1[2],ks_stat=ks1[3])
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
              bre<-list(chr=chr,start=n,stop=m,segment_p=ks2[1],p_value=ks2[4],mean_diff=ks2[2],ks_stat=ks2[3])
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
              bre<-list(chr=chr,start=n,stop=m,segment_p=ks3[1],p_value=ks3[4],mean_diff=ks3[2],ks_stat=ks3[3])
              breaks<-rbind(breaks,data.frame(bre))
            }
          }
        }
      }else{
        if(child==0){
          bre<-list(chr=chr,start=a,stop=b,segment_p=KS[1],p_value=KS[4],mean_diff=KS[2],ks_stat=KS[3])
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
      i<-which(bre$segment_p==min(bre$segment_p))
      max<-bre[i,]
      max<-max[1,]
      stack<-list(a=a,b=b,max=max)
      stacks<-append(stacks,list(stack))
      n<-a
      m<-max$start-1
      if(max$start>1&&n<=m){
        ks<-c(2,0,0,2)
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
        ks<-c(2,0,0,2)
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
          methY <- mean(as.numeric(as.matrix(y)), na.rm = TRUE)
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
          ks<-c(2,0,0,2)
          if(isTRUE((tmp_stop-tmp_start+1>=mincpgs)[1])&&isTRUE((calcSingleTrendAbs(XS,tmp_start,tmp_stop)>trend)[1])&&isTRUE((noValley(XS,tmp_start,tmp_stop,mincpgs,valley)==1)[1])){
            ks<-cortest(intput_dat,y, method, cov.mod,tmp_start,tmp_stop)
          }
          if(ks[1]<2){
            out<-data.frame(chr=chr,start=intput_dat$pos[tmp_start]-1,stop=intput_dat$pos[tmp_stop],q=-1,length=tmp_stop-tmp_start+1,ks_stat = ks[3],mean_diff=ks[2],p_value=ks[4])
            methX <- mean(as.numeric(as.matrix(intput_dat[tmp_start:tmp_stop, -c(1, 2)])), na.rm = TRUE)
            methY <- mean(as.numeric(as.matrix(y)), na.rm = TRUE)
            out$methX<-methX
            out$methY<-methY
            out$e_value<-calcEValue(intput_dat,y,tmp_start,tmp_stop)
            outputList<-rbind(outputList,as.data.frame(out))
          }
          tmp<-NULL
        }
        out<-data.frame(chr=chr,start=intput_dat$pos[b$start]-1,stop=intput_dat$pos[b$stop],q=-1,length=b$stop-b$start+1,ks_stat =b$ks_stat,mean_diff=b$mean_diff,p_value=b$p_value)
        methX <- mean(as.numeric(as.matrix(intput_dat[b$start:b$stop, -c(1, 2)])), na.rm = TRUE)
        methY <- mean(as.numeric(as.matrix(y)), na.rm = TRUE)
        out$methX<-methX
        out$methY<-methY
        out$e_value<-calcEValue(intput_dat,y,b$start,b$stop)
        outputList<-rbind(outputList,as.data.frame(out))
      }
    }
  }
  if(!is.null(tmp)){
    tmp<-collapseTmp(tmp)
    tmp_start<-as.integer(as.numeric(unlist(tmp$start))[1])
    tmp_stop<-as.integer(as.numeric(unlist(tmp$stop))[1])
    ks<-c(2,0,0,2)
    if(isTRUE((tmp_stop-tmp_start+1>=mincpgs)[1])&&isTRUE((calcSingleTrendAbs(XS,tmp_start,tmp_stop)>trend)[1])&&isTRUE((noValley(XS,tmp_start,tmp_stop,mincpgs,valley)==1)[1])){
      ks<-cortest(intput_dat,y, method, cov.mod,tmp_start,tmp_stop)
    }
    if(ks[1]<2){
      out<-data.frame(chr=chr,start=intput_dat$pos[tmp_start]-1,stop=intput_dat$pos[tmp_stop],q=-1,length=tmp_stop-tmp_start+1,ks_stat = ks[3],mean_diff=ks[2],p_value=ks[4])
      methX <- mean(as.numeric(as.matrix(intput_dat[tmp_start:tmp_stop, -c(1, 2)])), na.rm = TRUE)
      methY <- mean(as.numeric(as.matrix(y)), na.rm = TRUE)
      out$methX<-methX
      out$methY<-methY
      out$e_value<-calcEValue(intput_dat,y,tmp_start,tmp_stop)
      outputList<-rbind(outputList,as.data.frame(out))
    }
    tmp<-NULL
  }
  return(outputList)
}

segmentation<-function(intput_dat,y,cov.mod,chr,mincpgs,trend,valley,method){
  ks<-c(2,0,0,2)
  len<-nrow(intput_dat)
  XS<-calcSingleDiffSum(intput_dat,y)
  if(len-1>=mincpgs&&calcSingleTrendAbs(XS,1,len)>trend&&noValley(XS,1,len,mincpgs,valley)){
    ks<-cortest(intput_dat,y, method, cov.mod,1,len)
  }
  globalbreaks<-segmenterSTK(intput_dat,y,cov.mod,XS,1,len,chr,mincpgs,trend,valley,ks,method)
  outputList<-output(intput_dat,y,cov.mod,XS,globalbreaks,chr,mincpgs,trend,valley,method)
  return (outputList)
}
