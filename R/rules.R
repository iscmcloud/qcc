# Nelson rules & Western Electric rules & AIAG rules
#
# With Western Electric rules
# A process is out of control if either
# 1. One point plots outside 3-sigma control limits.
# 2. Two of three consecutive points plot beyond a 2-sigma limit.
# 3. Four of five consecutive points plot beyond a 1-sigma limit.
# 4. Eight consecutive points plot on one side of the center line.
#
# With Nelson rules
# A process is out of control if either
# 11. One point is more than 3 standard deviations from the mean.
# 12. Nine (or more) points in a row are on the same side of the mean.
# 13. Six (or more) points in a row are continually increasing (or decreasing).
# 14. Fourteen (or more) points in a row alternate in direction, increasing then decreasing.
# 15. Two (or three) out of three points in a row are more than 2 standard deviations from the mean in the same direction.
# 16. Four (or five) out of five points in a row are more than 1 standard deviation from the mean in the same direction.
# 17. Fifteen points in a row are all within 1 standard deviation of the mean on either side of the mean.
# 18. Eight points in a row exist, but none within 1 standard deviation of the mean, and the points are in both directions from the mean.
#
# With AIAG rules
# A process is out of control if either
# 21. One point is more than 3 standard deviations from the mean.
# 22. Seven (or more) points in a row are on the same side of the mean.
# 23. Six (or more) points in a row are continually increasing (or decreasing).
# 24. Fourteen (or more) points in a row alternate in direction, increasing then decreasing.
# 25. Two (or three) out of three points in a row are more than 2 standard deviations from the mean in the same direction.
# 26. Four (or five) out of five points in a row are more than 1 standard deviation from the mean in the same direction.
# 27. Fifteen points in a row are all within 1 standard deviation of the mean on either side of the mean.
# 28. Eight points in a row exist, but none within 1 standard deviation of the mean, and the points are in both directions from the mean.

qccRules <- function(object, rules = object$rules)
{
  # Return a vector of indices for cases (statistics & new.statistics) 
  # in object violating specified rules (NA if no rule is violated)
  rules <- as.numeric(rules)
  if(!(inherits(object, "qcc") | inherits(object, "mqcc")))
    stop("input object must be of class 'qcc' or 'mqcc'")
  stats <- c(object$statistics, object$newstats)
  out <- rep(NA, length(stats))
  # Western Electric rules 
  if(any(rules == 4)) 
  {  
    wer <- qccRulesViolatingWER4(object)
    out[wer] <- 4
  }
  if(any(rules == 3)) 
  {  
    wer <- qccRulesViolatingWER3(object)
    out[wer] <- 3
  }
  if(any(rules == 2)) 
  {  
    wer <- qccRulesViolatingWER2(object)
    out[wer] <- 2
  }
  if(any(rules == 1)) 
  {  
    wer <- qccRulesViolatingWER1(object)
    out[wer] <- 1
  }
  # Nelson rules
  if(any(rules == 18)) 
  {  
    wer <- qccRulesViolatingWER18(object)
    out[wer] <- 18
  }
  if(any(rules == 17)) 
  {  
    wer <- qccRulesViolatingWER17(object)
    out[wer] <- 17
  }
  if(any(rules == 16)) 
  {  
    wer <- qccRulesViolatingWER16(object)
    out[wer] <- 16
  }
  if(any(rules == 15)) 
  {  
    wer <- qccRulesViolatingWER15(object)
    out[wer] <- 15
  }
  if(any(rules == 14)) 
  {  
    wer <- qccRulesViolatingWER14(object)
    out[wer] <- 14
  }
  if(any(rules == 13)) 
  {  
    wer <- qccRulesViolatingWER13(object)
    out[wer] <- 13
  }
  if(any(rules == 12)) 
  {  
    wer <- qccRulesViolatingWER12(object)
    out[wer] <- 12
  }
  if(any(rules == 11)) 
  {  
    wer <- qccRulesViolatingWER11(object)
    out[wer] <- 11
  }
  attr(out, "WesternElectricRules") <- rules
  return(out)
}

qccRulesViolatingWER1 <- function(object, limits = object$limits)
{
  # Return cases beyond control limits (WER #1)
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}

qccRulesViolatingWER2 <- function(object, 
                                  run.points = 2,
                                  run.length = 3,
                                  k = object$nsigmas*2/3)
{
  # Return indices of points violating runs
  center     <- object$center
  statistics <- c(object$statistics, object$newstats)
  limits     <- paste("limits.", object$type, sep = "")
  limits     <- do.call(limits, list(center = object$center, 
                                     std.dev = object$std.dev,
                                     sizes = c(object$sizes, object$newsizes),
                                     nsigmas = k))
  i <- if(nrow(limits) > 1) seq(run.length, length(statistics)) else 1
  viol.above <- embed(statistics, run.length) > limits[i,2]
  viol.above <- which(apply(viol.above, 1, sum) >= run.points & viol.above[,1])
  viol.above <- viol.above + (run.length-1)
  viol.below <- embed(statistics, run.length) < limits[i,1]
  viol.below <- which(apply(viol.below, 1, sum) >= run.points & viol.below[,1])
  viol.below <- viol.below + (run.length-1)
  return(c(viol.above, viol.below))
}

qccRulesViolatingWER3 <- function(object, ...)
{
  qccRulesViolatingWER2(object, 
                        run.points = 4,
                        run.length = 5,
                        k = object$nsigmas*1/3)
}  

qccRulesViolatingWER4 <- function(object)
{
  # Return indices of points violating runs (WER #4)
  run.length <- 8
  center <- object$center
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  diffs <- statistics - center
  diffs[diffs > 0] <- 1
  diffs[diffs < 0] <- -1
  runs <- rle(diffs)
  vruns <- rep(runs$lengths >= run.length, runs$lengths)
  vruns.above <- (vruns & (diffs > 0))
  vruns.below <- (vruns & (diffs < 0))
  rvruns.above <- rle(vruns.above)
  rvruns.below <- rle(vruns.below)
  vbeg.above <- cumsum(rvruns.above$lengths)[rvruns.above$values] -
    (rvruns.above$lengths - run.length)[rvruns.above$values]
  vend.above <- cumsum(rvruns.above$lengths)[rvruns.above$values]
  vbeg.below <- cumsum(rvruns.below$lengths)[rvruns.below$values] -
    (rvruns.below$lengths - run.length)[rvruns.below$values]
  vend.below <- cumsum(rvruns.below$lengths)[rvruns.below$values]
  violators <- numeric()
  if (length(vbeg.above)) 
  { for (i in 1:length(vbeg.above))
    violators <- c(violators, vbeg.above[i]:vend.above[i]) }
  if (length(vbeg.below)) 
  { for (i in 1:length(vbeg.below))
    violators <- c(violators, vbeg.below[i]:vend.below[i]) }
  return(violators)
}



qccRulesViolatingWER15<- function(object, 
                                  run.points = 2,
                                  run.length = 3,
                                  k = object$nsigmas*2/3)
{
  # Return indices of points violating runs
  center     <- object$center
  statistics <- c(object$statistics, object$newstats)
  limits     <- paste("limits.", object$type, sep = "")
  limits     <- do.call(limits, list(center = object$center, 
                                     std.dev = object$std.dev,
                                     sizes = c(object$sizes, object$newsizes),
                                     nsigmas = k))
  i <- if(nrow(limits) > 1) seq(run.length, length(statistics)) else 1
  viol.above <- embed(statistics, run.length) > limits[i,2]
  viol.above <- which(apply(viol.above, 1, sum) >= run.points & viol.above[,1])
  viol.above <- viol.above + (run.length-1)
  viol.below <- embed(statistics, run.length) < limits[i,1]
  viol.below <- which(apply(viol.below, 1, sum) >= run.points & viol.below[,1])
  viol.below <- viol.below + (run.length-1)
  return(c(viol.above, viol.below))
}

qccRulesViolatingWER16 <- function(object, ...)
{
  qccRulesViolatingWER2(object, 
                        run.points = 4,
                        run.length = 5,
                        k = object$nsigmas*1/3)
} 


qccRulesViolatingWER12 <- function(object)
{
  # Return indices of points violating runs (WER #12)
  run.length <- 9
  center <- object$center
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  diffs <- statistics - center
  diffs[diffs > 0] <- 1
  diffs[diffs < 0] <- -1
  runs <- rle(diffs)
  vruns <- rep(runs$lengths >= run.length, runs$lengths)
  vruns.above <- (vruns & (diffs > 0))
  vruns.below <- (vruns & (diffs < 0))
  rvruns.above <- rle(vruns.above)
  rvruns.below <- rle(vruns.below)
  vbeg.above <- cumsum(rvruns.above$lengths)[rvruns.above$values] -
    (rvruns.above$lengths - run.length)[rvruns.above$values]
  vend.above <- cumsum(rvruns.above$lengths)[rvruns.above$values]
  vbeg.below <- cumsum(rvruns.below$lengths)[rvruns.below$values] -
    (rvruns.below$lengths - run.length)[rvruns.below$values]
  vend.below <- cumsum(rvruns.below$lengths)[rvruns.below$values]
  violators <- numeric()
  if (length(vbeg.above)) 
  { for (i in 1:length(vbeg.above))
    violators <- c(violators, vbeg.above[i]:vend.above[i]) }
  if (length(vbeg.below)) 
  { for (i in 1:length(vbeg.below))
    violators <- c(violators, vbeg.below[i]:vend.below[i]) }
  return(violators)
}

qccRulesViolatingWER18 <- function(object, ...)
{
  qccRulesViolatingWER2(object, 
                        run.points = 8,
                        run.length = 8,
                        k = object$nsigmas*1/3)
} 

qccRulesViolatingWER14 <- function(object) {
  # Return indices of points violating runs (WER #14)
  run.length <- 14
  statistics <- c(object$statistics, object$newstats)

  values <- statistics
  current_state <- ifelse(values[1] > values[2], "increasing", ifelse(values[1] < values[2], "decreasing", "unchanged"))
  counter_seq <- c(1)
  result_seq <- c()
  # Check all points
  for (i in 2:(length(values) - 1)) {
    next_state <- ifelse(values[i] > values[i + 1], "increasing", ifelse(values[i] < values[i + 1], "decreasing", "unchanged"))
    if (next_state != current_state & next_state != "unchanged") {
      counter_seq <- c(counter_seq, i)  
    } else {
      if (length(counter_seq) >= run.length ) {
        result_seq <- c(result_seq, counter_seq)
      }
      counter_seq <- c()	  
    }
    current_state <- next_state
  }

  violators <- numeric()
  violators <- c(result_seq)
  return(violators)
}

qccRulesViolatingWER13 <- function(object) {
  # Return indices of points violating runs (WER #13)
  run.length <- 6
  statistics <- c(object$statistics, object$newstats)

  values <- statistics
  current_state <- ifelse(values[1] > values[2], "increasing", ifelse(values[1] < values[2], "decreasing", "unchanged"))
  counter_seq <- c(1)
  result_seq <- c()
  # Check all points
  for (i in 2:(length(values) - 1)) {
    next_state <- ifelse(values[i] > values[i + 1], "increasing", ifelse(values[i] < values[i + 1], "decreasing", "unchanged"))
    if (next_state == current_state & next_state != "unchanged") {
      counter_seq <- c(counter_seq, i)  
    } else {
      if (length(counter_seq) >= run.length ) {
        result_seq <- c(result_seq, counter_seq)
      }
      counter_seq <- c(i)	  
    }
    current_state <- next_state
  }

  violators <- numeric()
  violators <- c(result_seq)
  return(violators)
}

qccRulesViolatingWER17 <- function(object, limits = object$limits)
{
  # Return cases beyond control limits (WER #17)
  violators <- numeric()
  return(violators)
}

qccRulesViolatingWER11 <- function(object, limits = object$limits)
{
  # Return cases beyond control limits (WER #11)
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}
