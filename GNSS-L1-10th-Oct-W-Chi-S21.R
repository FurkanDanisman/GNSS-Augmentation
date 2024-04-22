# Required Libraries
library(rcompanion) 
library(isotree)
library(dbscan)
library(scales)

#setwd("C://Users//-//Downloads")

dft <- read.csv("PR_L1_prefit-10thOct2022_S21.csv",header=T)

df <- read.csv("NSSE_PR_L1_prefit-10thOct2022_S21.csv",header=T)

# Chi-Square Test 

p <- 1
Chi_sqr <- matrix(0,nrow = nrow(dft),ncol = 1)
for (i in 1:nrow(dft)) {
  row_values <- dft[i, ]
  e_values <- row_values[!is.na(row_values)]
  if (length(e_values)<=2) {
    next
  }
  n <- length(e_values)
  th <- qchisq(0.05,n-p) / (n-p)
  if (df[i,]>=th) {
    Chi_sqr[i,] <- 1 
  }
}
wh <- which(Chi_sqr==1)

# Method-1 : Chi-square -> t-test #  (Manual Threshold ~ 0.90 - 0.97)

to_be_tested <- dft[wh,]
updated_values <- list()
dft2 <- dft
dft5 <- dft
g <- 0
h_t_ch <- list()

for (i in 1:nrow(to_be_tested)) {
  
  row_values <- to_be_tested[i, ]
  e_values <- row_values[!is.na(row_values)]
  
  if (length(e_values)==0) {
    next
  }
  e_SE <- sqrt(var(e_values) * (length(e_values) - 1) / (length(e_values) - 2))
  e_mean <- mean(e_values)
  n <- length(e_values)
  c0 <- qt(0.90, n - 1)
  c1 <- qt(0.99, n - 1)
  for (j in 1:length(e_values)) {
    e_Tj <- abs(e_values[j] - e_mean) / e_SE
    
    # t-test boundries 
    
    if (e_Tj > c0 & e_Tj < c1) {
      g[j] <- 1 
      alpha <- c0 / e_Tj * ((c1 - e_Tj) / (c1 - c0))^2
      dft2[which(dft2 == e_values[j], arr.ind = TRUE)] <- alpha * e_values[j]
      dft5[which(dft5 == e_values[j], arr.ind = TRUE)] <- 9999
      e_values[j] <- alpha * e_values[j]
      t_ch_s_N_t <- which(dft5 == 9999, arr.ind = TRUE)
    } else if (e_Tj >= c1) {
      g[j] <- 2
      dft2[which(dft2 == e_values[j], arr.ind = TRUE)] <- 9999
      e_values[j] <- 9999
    } else {
      g[j] <- 3
      next
    } 
  }
  
  h_t_ch[[i]] <- g
  g <- 0
  # Update the row in the original data frame
  updated_values[[i]] <- e_values
  
}

# Updated values
upd_res <- unlist(updated_values)

# Number of outliers
outlier_number_t_ch <- length(upd_res[upd_res==9999]) 

# Number of NA's
combined_column3 <- unlist(to_be_tested)
length(combined_column3) - length(upd_res)


upd_res <- upd_res[upd_res!=9999]

# Transformed Data Plots

par(mfrow=c(1,2))
plotNormalHistogram( combined_column3, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( upd_res, breaks=100,col="lightblue",prob = FALSE, 
                     main = "After the transformation", 
                     length = 100000,linecol = "red")


# Transformation Plots for Whole Data

t_ch_s_N <- which(dft2 == 9999, arr.ind = TRUE)
t_ch_s_N_t <- which(dft5 == 9999, arr.ind = TRUE)
dft_t_ch <- dft2[dft2!=9999]
combined_column <- unlist(dft)
combined_column2 <- unlist(dft_t_ch)

par(mfrow=c(1,2))

plotNormalHistogram( combined_column, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 

plotNormalHistogram( dft_t_ch, breaks=100,col="lightblue",prob = FALSE, 
                     main = "t-test w/Chi-Square Transformation", 
                     length = 100000,linecol = "red") 


# Transformation QQ-plots for Whole Data
qqnorm(combined_column)
abline(a=0,b=1,col="red")

qqnorm(dft_t_ch)
abline(a=0,b=1,col="red")

# Absolute Values 
abs_combined <- abs(combined_column)
abs_dft_t_ch <- abs(dft_t_ch)

hist(abs_combined,breaks = 100,col='lightblue',main = 'Before the transformation')
hist(abs_dft_t_ch,breaks = 100,col = 'orange',add=F,main = 'After the transformation')

# Example Comparison
par(mfrow=c(1,2))

#head(wh,100)

row_values <- dft[92, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "Passed Chi-Square Test | Epoch 92",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")

row_values <- dft[93, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "Passed Chi-Square Test | Epoch 93",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")

#which(h_t_ch[[1]]==1)
#which(h_t_ch[[80]]==2)
row_values <- to_be_tested[1, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "t-test Process at Epoch 1",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")
points(y = e_values[5], x = 5,col = "orange", pch = 16) 
text(y=e_values[5],x=5 ,labels = "Transformed", pos = 1, col = "orange")

row_values <- to_be_tested[80, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "t-test Process at Epoch 80",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")
points(y=e_values[9], x = 9,col = "navyblue", pch = 16) 
text(y=e_values[9],x=9, labels = "Outlier", pos = 3, col = "navyblue")


# Bonferroni Correction
# Error Rate: 1-(1-alpha)^C : C is the number of tests
(1-0.05)*(1-0.07) # 0.97
(1-0.05)*(1-0.01) # 0.99

# Each test is need to be done with around alpha / 2 to get overall 95% Confidence Level.
1-(1-0.025)^2


# Method-2 : Chi-square -> LOF # (Manual Threshold ~ 2.5 - 4.5)
to_be_tested <- dft[wh,]
updated_values <- list()
dft2 <- dft
dft5 <- dft
g <- 0
h_lof_ch <- list()
lof_points <- list()

for (i in 1:nrow(to_be_tested)) {
  row_values <- to_be_tested[i, ,drop=FALSE]
  e_values <- row_values[!is.na(row_values)]
  e <- as.matrix(e_values)
  
  if (is.null(nrow(e))) {
    next
  }
  if (nrow(e)<=2) {
    next
  }else if(nrow(e)<=4) {
    lof_score <- lof(e,minPts = nrow(e))
  }else{
    lof_score <- lof(e)
    
  }
  lof_points[[i]] <- lof_score
  for (j in 1:length(e_values)) {
    
    if (lof_score[j]>2.5 & lof_score[j]<=4) {
      g[j] <- 1
      alpha <- 2.5 / lof_score[j] * ((4.5 - lof_score[j]) / (4.5 - 2.5))^2
      dft2[which(dft2 == e_values[j], arr.ind = TRUE)] <- e_values[j]*alpha
      dft5[which(dft5 == e_values[j], arr.ind = TRUE)] <- 9999
      e_values[j] <- e_values[j]*alpha
    } else if (lof_score[j]>4) {
      g[j] <- 2
      dft2[which(dft2 == e_values[j], arr.ind = TRUE)] <- 9999
      e_values[j] <- 9999
    } else {
      g[j] <- 3
      next
    } 
  }
  h_lof_ch[[i]] <- g
  g <- 0
  # Update the row in the original data frame
  updated_values[[i]] <- e_values
  
}

# LOF Plot 
lof_points_unlist <- unlist(lof_points) 
hist(lof_points_unlist)
quantile(lof_points_unlist,0.97)

# Updated values
upd_res <- unlist(updated_values)

# Number of outliers
outlier_number_lof_ch <- length(upd_res[upd_res==9999]) 

# Number of NA's
combined_column3 <- unlist(to_be_tested)

length(combined_column3) - length(upd_res)


upd_res <- upd_res[upd_res!=9999]

# Transformed Data Plots

par(mfrow=c(1,2))
plotNormalHistogram( combined_column3, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( upd_res, breaks=100,col="lightblue",prob = FALSE, 
                     main = "After the transformation", 
                     length = 100000,linecol = "red")


# Transformation Plots for Whole Data 
lof_ch_s_N <- which(dft2 == 9999, arr.ind = TRUE)
lof_ch_s_N_t <- which(dft5 == 9999, arr.ind = TRUE)

dft_lof_ch <- dft2[dft2!=9999]
combined_column <- unlist(dft)
combined_column2 <- unlist(dft_lof_ch)

par(mfrow=c(1,2))

plotNormalHistogram( combined_column, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( dft_lof_ch, breaks=100,col="lightblue",prob = FALSE, 
                     main = "LOF Transformation", 
                     length = 100000,linecol = "red") 

# Transformation QQ-plots for Whole Data
qqnorm(combined_column)
abline(a=0,b=1,col="red")

qqnorm(dft_lof_ch)
abline(a=0,b=1,col="red")

# Absolute Values 
abs_combined <- abs(combined_column)
abs_dft_lof_ch <- abs(dft_lof_ch)

hist(abs_combined,breaks = 100,col='lightblue',main = 'Before the transformation')
hist(abs_dft_lof_ch,breaks = 100,col = 'orange',add=F,main = 'After the transformation')



# Example Comparison
par(mfrow=c(1,2))

# h_lof_ch

row_values <- to_be_tested[1, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "LOF Process at Epoch 1",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")


row_values <- to_be_tested[80, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "LOF Process at Epoch 80",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")
points(y = e_values[c(9,10,11)], x = c(9,10,11),col = "orange", pch = 16) 
points(y = e_values[c(12)], x = c(12),col = "navyblue", pch = 16) 
text(y=e_values[c(9,10,11)],x=c(9,10,11) ,labels = "Transformed", pos = 3, col = "orange")
text(y=e_values[c(12)],x=c(12) ,labels = "Outlier", pos = 1, col = "navyblue")

# Method-3 : Chi-square -> Isolation Forest (Manual Threshold ~0.6)

to_be_tested <- dft[wh,]
updated_values <- list()
dft2 <- dft
g <- 0
h_if_ch <- list()
if_scores <- list()


for (i in 1:nrow(to_be_tested)) {
  row_values <- to_be_tested[i, ,drop=FALSE]
  e_values <- row_values[!is.na(row_values)]
  e <- as.matrix(e_values)
  
  if (is.null(nrow(e))) {
    next
  }
  model <- isolation.forest(e,ntrees = 100)
  if_score <- predict.isolation_forest(model,e)
  if_scores[[i]] <- if_score
  
  for (j in 1:length(e_values)) {
    
    if (if_score[j]>0.6) {
      g[j] <- 1
      dft2[which(dft2 == e_values[j], arr.ind = TRUE)] <- 9999
      e_values[j] <- 9999
    } else {
      g[j] <- 2
      next
    } 
  }
  h_if_ch[[i]] <- g
  g <- 0
  # Update the row in the original data frame
  updated_values[[i]] <- e_values
  
}

# IF Plot 
if_scores_unlist <- unlist(if_scores)
hist(if_scores_unlist)
quantile(if_scores_unlist,0.95)

# Updated values
upd_res <- unlist(updated_values)

# Number of outliers
outlier_number_if_ch <- length(upd_res[upd_res==9999]) 

# Number of NA's
combined_column3 <- unlist(to_be_tested)
length(combined_column3) - length(upd_res)


upd_res <- upd_res[upd_res!=9999]

# Transformed Data Plots

par(mfrow=c(1,2))
plotNormalHistogram( combined_column3, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( upd_res, breaks=100,col="lightblue",prob = FALSE, 
                     main = "After the transformation", 
                     length = 100000,linecol = "red")


# Transformation Plots for Whole Data 
if_ch_s_N <- which(dft2 == 9999, arr.ind = TRUE)
dft_if_ch <- dft2[dft2!=9999]
combined_column <- unlist(dft)
combined_column2 <- unlist(dft_if_ch)

par(mfrow=c(1,2))

plotNormalHistogram( combined_column, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( dft_if_ch, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Isolation Forest transformation", 
                     length = 100000,linecol = "red") 

# Transformation QQ-plots for Whole Data
qqnorm(combined_column)
abline(a=0,b=1,col="red")

qqnorm(dft_if_ch)
abline(a=0,b=1,col="red")


# Absolute Values 
abs_combined <- abs(combined_column)
abs_dft_if_ch <- abs(dft_if_ch)

hist(abs_combined,breaks = 100,col='lightblue',main = 'Before the transformation')
hist(abs_dft_if_ch,breaks = 100,col = 'orange',add=F,main = 'After the transformation')



# Example Comparison
par(mfrow=c(1,2))

# h_if_ch
#which(h_if_ch[[1]]==1)
row_values <- to_be_tested[1, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "IF Process at Epoch 1",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")
points(y=e_values[c(5)], x = c(5),col = "navyblue", pch = 16) 
text(y=e_values[c(5)],x=c(5), labels = "Outlier", pos = 3, col = "navyblue")

#which(h_if_ch[[80]]==1)
row_values <- to_be_tested[80, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "IF Process at Epoch 80",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")
points(y = e_values[c(9,12)], x = c(9,12),col = "navyblue", pch = 16) 
text(y=e_values[c(9,12)],x=c(9,12) ,labels = "Outlier", pos = 3, col = "navyblue")



# Method-4 : Chi-square -> Isolation Forest -> LOF # (Manual Threshold ~ 0.6 - 0.4 -> 2.2 - 3.15)

to_be_tested <- dft[wh,]
updated_values <- list()
dft2 <- dft
dft5 <- dft
g <- 0
h_if_lof_ch <- list()
if_lof_scores <- list()

for (i in 1:nrow(to_be_tested)) {
  row_values <- to_be_tested[i, ,drop=FALSE]
  e_values <- row_values[!is.na(row_values)]
  e <- as.matrix(e_values)
  
  if (is.null(nrow(e))) {
    next
  }
  
  # First Process
  model <- isolation.forest(e,ntrees = 100)
  if_score <- predict.isolation_forest(model,e)
  # Second Process
  l_e <- e[if_score>0.39,,drop=FALSE]
  if (nrow(l_e)<=2) {
    for (j in 1:length(e_values)) {
      
      if (if_score[j]>0.6) {
        g[j] <- 1
        dft2[which(dft2 == e_values[j], arr.ind = TRUE)] <- 9999
        e_values[j] <- 9999
      } else {
        g[j] <- 2
      } 
    }
    h_if_lof_ch[[i]] <- g
    g <- 0
    updated_values[[i]] <- e_values
    if_lof_scores[[i]] <- NA
    next
  }else if(nrow(l_e)<=4) {
    lof_score <- lof(l_e,minPts = nrow(l_e))
    if_lof_scores[[i]] <- lof_score
  }else{
    lof_score <- lof(l_e)
    if_lof_scores[[i]] <- lof_score
  }
  for (j in 1:length(l_e)) {
    if (lof_score[j]>2.2 & 3.15  >= lof_score[j]) {
      g[j] <- 1
      alpha <- 2.2 / lof_score[j] * ((3.15 - lof_score[j]) / (3.15 - 2.2))^2
      value <- e[which(e == l_e[j], arr.ind = TRUE)]
      dft2[which(dft2 == value, arr.ind = TRUE)] <- value*alpha
      dft5[which(dft5 == value, arr.ind = TRUE)] <- 9999
      e_values[which(e_values == l_e[j], arr.ind = TRUE)] <- value*alpha
    }else if (lof_score[j]>3.15) {
      g[j] <- 2
      value <- e[which(e == l_e[j], arr.ind = TRUE)]
      dft2[which(dft2 == value, arr.ind = TRUE)] <- 9999
      e_values[which(e_values == l_e[j], arr.ind = TRUE)] <- 9999
    } else {
      g[j] <- 3
    }
  }
  
  h_if_lof_ch[[i]] <- g
  g <- 0
  # Update the row in the original data frame
  updated_values[[i]] <- e_values
  
}

# IF-LOF Scores
if_lof_scores_unlist <- unlist(if_lof_scores)
if_lof_scores_unlist <- na.omit(if_lof_scores_unlist)
length(if_lof_scores_unlist)

# Updated values
upd_res_if_lof <- unlist(updated_values)


# Number of outliers
outlier_number_if_lof_ch <- length(upd_res_if_lof[upd_res_if_lof==9999]) 

upd_res_if_lof <- upd_res_if_lof[upd_res_if_lof!=9999]

# Transformed Data Plots

par(mfrow=c(1,2))
plotNormalHistogram( combined_column3, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( upd_res_if_lof, breaks=100,col="lightblue",prob = FALSE, 
                     main = "After the transformation", 
                     length = 100000,linecol = "red")


# Transformation Plots for Whole Data
if_lof_ch_s_N <- which(dft2 == 9999, arr.ind = TRUE)
if_lof_ch_s_N_t <- which(dft5 == 9999, arr.ind = TRUE)
dft_if__lof_ch <- dft2[dft2!=9999]
combined_column <- unlist(dft)
combined_column2 <- unlist(dft_if__lof_ch)

par(mfrow=c(1,2))

plotNormalHistogram( combined_column, breaks=100,col="lightblue",prob = FALSE, 
                     main = "Before the transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( dft_if__lof_ch, breaks=100,col="lightblue",prob = FALSE, 
                     main = "IF-LOF-transformation", 
                     length = 100000,linecol = "red") 


# Transformation QQ-plots for Whole Data
qqnorm(combined_column)
abline(a=0,b=1,col="red")

qqnorm(dft_if__lof_ch)
abline(a=0,b=1,col="red")

# Absolute Values 
abs_combined <- abs(combined_column)
abs_dft_if__lof_ch <- abs(dft_if__lof_ch)

hist(abs_combined,breaks = 100,col='lightblue',main = 'Before the transformation')
hist(abs_dft_if__lof_ch,breaks = 100,col = 'orange',add=F,main = 'After the transformation')



# Example Comparison
par(mfrow=c(1,2))

# h_if_lof_ch

row_values <- to_be_tested[1, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "IF-LOF Process at Epoch 1",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")


#which(h_if_lof_ch[[80]]==2)
row_values <- to_be_tested[80, ]
e_values <- row_values[!is.na(row_values)]
plot(e_values,main = "IF Process at Epoch 80",ylim = c(-30,30),xlab="Satellites",ylab="Standardized Residuals")
abline(h=0,col="red")
points(y = e_values[c(9,10,11,12)], x = c(9,10,11,12),col = "navyblue", pch = 16) 
text(y=e_values[c(9,10,11,12)],x=c(9,10,11,12) ,labels = "Outlier", pos = 3, col = "navyblue")



# Number of Outlier for each method # 
combined_column <- unlist(dft)
wo_na <- na.omit(combined_column)
Total_obs <- length(wo_na)
uht <- unlist(h_t_ch)
uhl <- unlist(h_lof_ch)
uhif <- unlist(h_if_ch)
uhif_lof <- unlist(h_if_lof_ch)
count_t_t <- sum(uht == 1)
count_t_lof <- sum(uhl == 1)
count_t_iflof <- sum(uhif_lof == 1)
outlier <- c(outlier_number_t_ch,outlier_number_lof_ch,outlier_number_if_ch,outlier_number_if_lof_ch)
outlier_df <- data.frame("Detected Outlier"=outlier)
outlier_df$Transformed <- c(count_t_t,count_t_lof,NA,count_t_iflof)
rownames(outlier_df) <- c("t-test","LOF","IF", "IF-LOF")
outlier_df <- cbind(outlier_df,rep(Total_obs,4))
colnames(outlier_df) <- c("Detected Outlier","Transformed" ,"Total Observation")
outlier_df$Outlier_Ratio <- round(outlier_df$`Detected Outlier` / outlier_df$`Total Observation`,4)
outlier_df$Transformed_Ratio <- round(outlier_df$Transformed / outlier_df$`Total Observation`,4)
outlier_df$Outlier_Percentage <- outlier_df$Outlier_Ratio * 100
outlier_df$Transformed_Percentage <- outlier_df$Transformed_Ratio * 100
outlier_df <- outlier_df[,-c(4:5)]
outlier_df

# Plot Comparison #

par(mfrow=c(2,2))

plotNormalHistogram( dft_t_ch, breaks=200,col="lightblue",prob = FALSE, 
                     main = "t-test w/Chi-Square Transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( dft_lof_ch, breaks=200,col="lightblue",prob = FALSE, 
                     main = "LOF w/Chi-Square Transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( dft_if_ch, breaks=200,col="lightblue",prob = FALSE, 
                     main = "IF w/Chi-Square Transformation", 
                     length = 100000,linecol = "red") 
plotNormalHistogram( dft_if__lof_ch, breaks=200,col="lightblue",prob = FALSE, 
                     main = "IF-LOF w/Chi-Square Transformation", 
                     length = 100000,linecol = "red") 


## QQ-Plots ## 

par(mfrow=c(2,2))
qqnorm(dft_t_ch,main = "t-test")
qqline(dft_t_ch,col="red")
qqnorm(dft_lof_ch,main = "LOF")
qqline(dft_lof_ch,col="red")
qqnorm(dft_if_ch,main = "IF")
qqline(dft_if_ch,col="red")
qqnorm(dft_if__lof_ch,main = "IF-LOF")
qqline(dft_if__lof_ch,col="red")

outlier_df

# OUTLIER

# t-test Satellite Detector 

t_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==t_ch_s_N[,1])) {
    a <- t_ch_s_N[i==t_ch_s_N[,1],,drop=FALSE]
    t_s_m[[i]] <- a[,2]
  }else{t_s_m[[i]] <- 0}
}

s_i1 <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(t_s_m[[i]]%in%j)){
      s_i1[i,j] <- j
    }
  }
}

non_zero_cols <- apply(s_i1 != 0, 2, any)

# Extract columns that are not all zeros
result_t <- s_i1[, non_zero_cols]

# Lof Satellite Detector
l_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==lof_ch_s_N[,1])) {
    a <- lof_ch_s_N[i==lof_ch_s_N[,1],,drop=FALSE]
    l_s_m[[i]] <- a[,2]
  }else{l_s_m[[i]] <- 0}
}

s_i <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(l_s_m[[i]]%in%j)){
      s_i[i,j] <- j
    }
  }
}

non_zero_cols <- apply(s_i != 0, 2, any)

# Extract columns that are not all zeros
result_lof <- s_i[, non_zero_cols]

plot(x = 1, y = result_lof[1,1], ylim = c(0,140),xlim=c(0,1700) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="t-test vs LOF: Outlier Detection of the Satellites for Each Epoch  ")

for (col in 1:max(ncol(result_t),ncol(result_lof))) {
  
  x <- 1:nrow(s_i1)
  
  if (col > ncol(result_t)) {
    y2=result_lof[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)}
    }
  }else if(col > ncol(result_lof)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[k], y = y[k],col=alpha("orange",0.25),pch=16)}
    }
  }else{
    y = result_t[,col]
    y2=result_lof[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16)
      }else{
        points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16) 
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "LOF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))

# With only first 20 Epoch 

plot(x = 1, y = result_lof[1,1], ylim = c(0,140),xlim=c(0,25) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="Outlier Detection on the Satellites for first 20 Epoch")

for (col in 1:max(ncol(result_t),ncol(result_lof))) {
  
  x <- 1:20
  
  if (col > ncol(result_t)) {
    y2=result_lof[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)
        text(paste("",y2[k]),x = x[k], y = y2[k],col="navyblue",pos = 3)
        text(paste("",x[k]),x = x[k], y = 0,col="red",pos = 1)
        abline(v=x[k],col="red",lty=2)}
    }
  }else if(col > ncol(result_lof)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[j], y = y[k],col="orange",pch=16)
        text(paste("",y[k]),x = x[j], y = y[k],col="orange",pos = 2)
        text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
        abline(v=x[j],col="red",lty=2)}
    }
  }else{
    y = result_t[,col]
    y2=result_lof[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col="navyblue",pch=15)
          text(paste("",y2[j]),x = x[j], y = y2[j],col="navyblue",pos = 3)
          text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
          abline(v=x[j],col="red",lty=2)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col="orange",pch=16)
        text(paste("",y[j]),x = x[j], y = y[j],col="orange",pos = 2)
        text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
        abline(v=x[j],col="red",lty=2)
      }else{
        points(x = x[j], y = y2[j],col="navyblue",pch=15)
        text(paste("",y2[j]),x = x[j], y = y2[j],col="navyblue",pos = 3)
        points(x = x[j], y = y[j],col="orange",pch=16)
        text(paste("",y[j]),x = x[j], y = y[j],col="orange",pos = 2)
        text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
        abline(v=x[j],col="red",lty=2)
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "LOF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))


# If Satellite Detector 
if_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==if_ch_s_N[,1])) {
    a <- if_ch_s_N[i==if_ch_s_N[,1],,drop=FALSE]
    if_s_m[[i]] <- a[,2]
  }else{if_s_m[[i]] <- 0}
}
s_i <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(if_s_m[[i]]%in%j)){
      s_i[i,j] <- j
    }
  }
}


non_zero_cols <- apply(s_i != 0, 2, any)

# Extract columns that are not all zeros
result_if <- s_i[, non_zero_cols]

plot(x = 1, y = result_if[1,1], ylim = c(0,140),xlim=c(0,1700) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="t-test vs IF: Outlier Detection of the Satellites for Each Epoch ")

for (col in 1:max(ncol(result_t),ncol(result_if))) {
  
  x <- 1:nrow(s_i1)
  
  if (col > ncol(result_t)) {
    y2=result_if[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)}
    }
  }else if(col > ncol(result_if)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[k], y = y[k],col=alpha("orange",0.25),pch=16)}
    }
  }else{
    y = result_t[,col]
    y2=result_if[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16)
      }else{
        points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16) 
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "IF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))

# If-LOF Satellite Detector

if_lof_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==if_lof_ch_s_N[,1])) {
    a <- if_lof_ch_s_N[i==if_lof_ch_s_N[,1],,drop=FALSE]
    if_lof_s_m[[i]] <- a[,2]
  }else{if_lof_s_m[[i]] <- 0}
}

s_i <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(if_lof_s_m[[i]]%in%j)){
      s_i[i,j] <- j
    }
  }
}

non_zero_cols <- apply(s_i != 0, 2, any)

# Extract columns that are not all zeros
result_if_lof <- s_i[, non_zero_cols]

plot(x = 1, y = result_if_lof[1,1], ylim = c(0,140),xlim=c(0,1700) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="t-test vs IF-LOF: Outlier Detection of the Satellites for Each Epoch")

for (col in 1:max(ncol(result_t),ncol(result_if_lof))) {
  
  x <- 1:nrow(s_i1)
  
  if (col > ncol(result_t)) {
    y2=result_if_lof[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)}
    }
  }else if(col > ncol(result_if_lof)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[k], y = y[k],col=alpha("orange",0.25),pch=16)}
    }
  }else{
    y = result_t[,col]
    y2=result_if_lof[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16)
      }else{
        points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16) 
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "IF-LOF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))


# Transformed

# t-test Satellite Detector 

t_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==t_ch_s_N_t[,1])) {
    a <- t_ch_s_N_t[i==t_ch_s_N_t[,1],,drop=FALSE]
    t_s_m[[i]] <- a[,2]
  }else{t_s_m[[i]] <- 0}
}

s_i1 <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(t_s_m[[i]]%in%j)){
      s_i1[i,j] <- j
    }
  }
}

non_zero_cols <- apply(s_i1 != 0, 2, any)

# Extract columns that are not all zeros
result_t <- s_i1[, non_zero_cols]

# Lof Satellite Detector
l_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==lof_ch_s_N_t[,1])) {
    a <- lof_ch_s_N_t[i==lof_ch_s_N_t[,1],,drop=FALSE]
    l_s_m[[i]] <- a[,2]
  }else{l_s_m[[i]] <- 0}
}

s_i <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(l_s_m[[i]]%in%j)){
      s_i[i,j] <- j
    }
  }
}

non_zero_cols <- apply(s_i != 0, 2, any)

# Extract columns that are not all zeros
result_lof <- s_i[, non_zero_cols]

plot(x = 1, y = result_lof[1,1], ylim = c(0,140),xlim=c(0,1700) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="t-test vs LOF Transformed Detection of the Satellites for Each Epoch")

for (col in 1:max(ncol(result_t),ncol(result_lof))) {
  
  x <- 1:nrow(s_i1)
  
  if (col > ncol(result_t)) {
    y2=result_lof[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)}
    }
  }else if(col > ncol(result_lof)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[k], y = y[k],col=alpha("orange",0.25),pch=16)}
    }
  }else{
    y = result_t[,col]
    y2=result_lof[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16)
      }else{
        points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16) 
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "LOF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))

# With only first 20 Epoch 

plot(x = 1, y = result_lof[1,1], ylim = c(0,140),xlim=c(0,25) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="Transformed Detection on the Satellites for first 20 Epoch")

for (col in 1:max(ncol(result_t),ncol(result_lof))) {
  
  x <- 1:20
  
  if (col > ncol(result_t)) {
    y2=result_lof[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)
        text(paste("",y2[k]),x = x[k], y = y2[k],col="navyblue",pos = 3)
        text(paste("",x[k]),x = x[k], y = 0,col="red",pos = 1)
        abline(v=x[k],col="red",lty=2)}
    }
  }else if(col > ncol(result_lof)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[k], y = y[k],col="orange",pch=16)
        text(paste("",y[k]),x = x[k], y = y[k],col="orange",pos = 2)
        text(paste("",x[k]),x = x[k], y = 0,col="red",pos = 1)
        abline(v=x[k],col="red",lty=2)}
    }
  }else{
    y = result_t[,col]
    y2=result_lof[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col="navyblue",pch=15)
          text(paste("",y2[j]),x = x[j], y = y2[j],col="navyblue",pos = 3)
          text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
          abline(v=x[j],col="red",lty=2)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col="orange",pch=16)
        text(paste("",y[j]),x = x[j], y = y[j],col="orange",pos = 2)
        text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
        abline(v=x[j],col="red",lty=2)
      }else{
        points(x = x[j], y = y2[j],col="navyblue",pch=15)
        text(paste("",y2[j]),x = x[j], y = y2[j],col="navyblue",pos = 3)
        points(x = x[j], y = y[j],col="orange",pch=16)
        text(paste("",y[j]),x = x[j], y = y[j],col="orange",pos = 2)
        text(paste("",x[j]),x = x[j], y = 0,col="red",pos = 1)
        abline(v=x[j],col="red",lty=2)
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "LOF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))


# If-LOF Satellite Detector

if_lof_s_m <- list()
for (i in 1:nrow(dft)) {
  if (any(i==if_lof_ch_s_N_t[,1])) {
    a <- if_lof_ch_s_N_t[i==if_lof_ch_s_N_t[,1],,drop=FALSE]
    if_lof_s_m[[i]] <- a[,2]
  }else{if_lof_s_m[[i]] <- 0}
}

s_i <- matrix(0,nrow = nrow(dft),ncol = ncol(dft2))
for (i in 1:nrow(dft)) {
  for (j in 1:ncol(dft)) {
    if(any(if_lof_s_m[[i]]%in%j)){
      s_i[i,j] <- j
    }
  }
}

non_zero_cols <- apply(s_i != 0, 2, any)

# Extract columns that are not all zeros
result_if_lof <- s_i[, non_zero_cols]

plot(x = 1, y = result_if_lof[1,1], ylim = c(0,140),xlim=c(0,1700) ,type = "p", col = "white", xlab = "Epoch", ylab = "Satellite",
     main="t-test vs IF-LOF Transformed Detection of the Satellites for Each Epoch")

for (col in 1:max(ncol(result_t),ncol(result_if_lof))) {
  
  x <- 1:nrow(s_i1)
  
  if (col > ncol(result_t)) {
    y2=result_if_lof[,col]
    for (k in 1:length(y2)) {
      if(y2[k]==0) next
      else{points(x = x[k], y = y2[k],col=alpha("navyblue",0.2),pch=15)}
    }
  }else if(col > ncol(result_if_lof)){
    y=result_t[,col]
    for (k in 1:length(y)) {
      if(y[k]==0) next
      else{points(x = x[k], y = y[k],col=alpha("orange",0.25),pch=16)}
    }
  }else{
    y = result_t[,col]
    y2=result_if_lof[,col]
    for (j in 1:length(y)) {
      if (y[j]==0) {
        if (y2[j]==0) {
          next
        }else{points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)}
      }else if(y2[j]==0){
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16)
      }else{
        points(x = x[j], y = y2[j],col=alpha("navyblue",0.2),pch=15)
        points(x = x[j], y = y[j],col=alpha("orange",0.25),pch=16) 
      }
    }
    
  }
}

legend("topleft", 
       inset = c(0,0), 
       cex = 1, 
       bty = "n", 
       legend = c("t-test", "IF-LOF"), 
       text.col = c("orange", "navyblue"),
       col = c("orange", "navyblue"), 
       pch = c(16,15))
