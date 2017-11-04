library(R.matlab)
library(stats)
library(huge)

# load data
n <- as.integer(param[1])
p <- as.integer(param[2])
exNum <- as.integer(param[3])
indir <- param[4]
outdir <- param[5]

#data_file <- paste(indir,"/data_",n,"_",p,".mat", sep = "", collapse = NULL)
data_file <- paste(indir,"/ecoli_",n,"_splits.mat", sep = "", collapse = NULL)
data <- readMat(data_file)


#sigmas <- list()
#mus <- list()
lambda <- c(1e-5, 0.0001, 0.001, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,  0.5, 0.7, 0.8, 0.9, 1, 2, 4)
lambda <- rev(lambda)

precision <- matrix(0, exNum, length(lambda))
recall <- matrix(0, exNum, length(lambda))
fpr <- matrix(0, exNum, length(lambda))

for(i in 1:exNum)
{
  xTrain <- data$xTrain[[i]][[1]]
  # gaussianize the data
  xTrain_G <- huge.npn(x = xTrain)
  output <- huge(x = xTrain_G, lambda = lambda, method = "glasso")
  
  for(j in 1:length(lambda))
  {
    adj_pred = output$icov[[j]]
    adj_pred[which(adj_pred != 0)] = 1; diag(adj_pred) = 0;
    adj <- data$adj
    precision[i, j] = length(which(adj != 0 & adj_pred != 0))/length(which(adj_pred != 0))
    recall[i, j] = length(which(adj != 0 & adj_pred != 0))/length(which(adj != 0))
    fpr[i, j] = length(which(adj == 0 & adj_pred != 0))/(length(which(adj == 0))-p) 
  }
}

# to find the best lambda, perform CV
# t_folds <- 5
# for(i in 1:exNum)
# {
#   x <- data$xTrain[[i]][[1]]
#   x <- huge.npn(x = x)
#   cv_score <- matrix(0,t_folds, length(lambda))
#   for(folds in 1:t_folds)
#   {
#     n5 <- round(n/5)
#     tst_idx <- ((folds-1)*n5+1) : (folds * n5)
#     tr_idx <- setdiff(1:n, tst_idx)
#     xtr <- x[tr_idx,]
#     xtst <- x[tst_idx,]
#     output_cv <- huge(x = xtr, lambda = lambda, method = "glasso")
#     mu = colMeans(xtst) # ideally use mu estimated from train data
#     xtst_c <- t(t(xtst)-mu)
#     S <- t(xtst_c)%*%xtst_c/length(tst_idx)
#     for(j in 1:length(lambda))
#     {
#       cv_score[folds, j] <- log(det(output_cv$icov[[j]])) - sum(S*output_cv$icov[[j]])
#     }
#   }
#   cv_score <- colMeans(cv_score, na.rm = TRUE)
#   best_idx <- which.max(cv_score)
#   best_lambda <- lambda[best_idx]
#   best_icov = output$icov[[best_lambda]]
# }


#outfile <- paste(outdir,"/results_",n,"_",p,".mat", sep = "", collapse = NULL)
outfile <- paste(outdir,"/results_",n,".mat", sep = "", collapse = NULL)
#writeMat(outfile, precision = precision, recall = recall, fpr = fpr, icov = best_icov)
writeMat(outfile, precision = precision, recall = recall, fpr = fpr)
