##imputdata: pro=protein_data, s=sex, a=age, lab=label

res_func<-function(a,s,lab,pro){
  lm_model <-lm(lab ~ a + s)
  Y.res <- residuals(lm_model)
  ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,savePredictions = "final")
  pls_model <- train(Y.res ~ ., data = data.frame(pro, Y.res),method = "pls",tuneLength = 5, trControl = ctrl,scale = TRUE)
  calculate_vip <- function(pls_model) {
    w <- pls_model$loading.weights
    s <- pls_model$Yloadings
    ss <- colSums(s^2)
    vip <- sqrt(nrow(w) * rowSums((w^2) %*% diag(ss)) / sum(ss))
    names(vip) <- colnames(pls_model$model$X)
    return(vip)
  }
  cv_folds <- pls_model$pred$Resample  # 获取交叉验证的分组信息
  vip_list <- list()
  for (fold in unique(cv_folds)) {
    fold_data <- pls_model$pred[cv_folds == fold, ]
    model_fold <- plsr(Y.res ~ ., data = data.frame(pro_facial[,tolower(pro_id$protein)], Y.res)[fold_data$rowIndex, ], ncomp = pls_model$bestTune$ncomp)
    vip_list[[fold]] <- calculate_vip(model_fold)
  }
  vip_df <- do.call(rbind, vip_list) |>
    as.data.frame() |>
    tibble::rownames_to_column("Fold") |>
    tidyr::pivot_longer(cols = -Fold, names_to = "Feature", values_to = "VIP")
    vip_df$Feature<-rep(pro_id$protein,50)
}
