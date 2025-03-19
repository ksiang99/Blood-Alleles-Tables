Results of Random Forest using randomForest package and the following parameters:
```r
"random_forest" = {
  class_weights <- 1 / table(df_train[[target_column]])
  class_weights <- setNames(as.numeric(class_weights), names(class_weights))
  randomForest(x = df_train[predictor_columns], y = df_train[[target_column]], 
  ntree = 500,
  classwt = class_weights)
