Results of Random Forest using the randomForest package and the following parameters:
```r
"random_forest" = {
min_class_size <- min(table(df_train[[target_column]]))
balanced_sampsize <- rep(min_class_size, length(unique(df_train[[target_column]])))

randomForest(x = df_train[predictor_columns], y = df_train[[target_column]], 
        ntree = 500,
        sampsize = balanced_sampsize,
        strata = df_train[[target_column]])
}
