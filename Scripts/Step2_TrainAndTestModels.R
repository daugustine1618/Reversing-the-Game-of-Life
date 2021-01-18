
# Train Naive Bayes model

p_y_hat <- train_data %>% dplyr::summarize(across(start_0:start_624, mean))

S_stop.N.indices <- list(); p_x.y_hat <- list(); train_data.subset <- list(); X.test <- list(); p_y.x_hat.Bayes.test <- list()

for(k in 1:625){
  # As was done in function definition for define_aux_matrix(), determine indices (i,j) and (I,J) relating index for each cell and indices for cells in its neighborhood

  j <- ifelse(k %% 25 == 0, 25, k %% 25)
  i <- ((k - j) + 25)/25
  J <- ifelse((j - 1):(j + 1) %% 25 == 0, 25, (j - 1):(j + 1) %% 25)
  I <- ifelse((i - 1):(i + 1) %% 25 == 0, 25, (i - 1):(i + 1) %% 25)

  # Define vector of indices for neighboring cells of kth cell

  S_stop.N.indices[[k]] <- plyr::mdply(expand.grid(x = I, y = J), function(x, y){(x - 1)*25 + y}) %>% .[3] %>% pull()

  # Define subset of data corresponding to kth cell that includes the observations of the kth starting grid cell value, delta, the kth stopping grid cell value, and the corresponding
  # neighborhood sum of cells in the stopping grid

  train_data.subset[[k]] <- train_data %>% dplyr::select(c(delta, 2 + k, 627 + S_stop.N.indices[[k]]))
  train_data.subset[[k]] <- train_data.subset[[k]] %>% mutate(X_3 = train_data.subset[[k]] %>% as.matrix() %>% .[, 3:11] %>% rowSums()) %>% dplyr::select(c(2, 1, 7, 12)) %>%
    setNames(c("Y", "X_1", "X_2", names(.)[4]))

  p_x.y_hat[[k]] <- list()

  # Given above defined subset of data pertaining to kth cell, define PMFs for X|Y = 0 and X|Y = 1

  for(m in 0:1){
    p_x.y_hat[[k]][[paste(m)]] <- train_data.subset[[k]] %>% filter(Y == m) %>% dplyr::select(-Y) %>% group_by(X_1, X_2, X_3) %>% summarize(p = n()/nrow(.), .groups = "drop") %>%
      arrange(X_1, X_2, X_3)
  }

  # Define subset of predictors pertaining to kth starting grid cell

  X.test[[k]] <- test_data %>% dplyr::select(c(delta, 2 + S_stop.N.indices[[k]]))
  X.test[[k]] <- X.test[[k]] %>% mutate(X_3 = X.test[[k]] %>% as.matrix() %>% .[, 2:10] %>% rowSums()) %>% dplyr::select(c(1, 6, 11)) %>% setNames(c("X_1", "X_2", names(.)[3]))

  # Define PMF values for Y|X = x for starting grid cells in test dataset

  p_y.x_hat.Bayes.test[[k]] <- p_y_hat[[k]]*ifelse(is.na(X.test[[k]] %>% left_join(p_x.y_hat[[k]][["1"]], by = c("X_1", "X_2", "X_3")) %>% .$p), 0,
    X.test[[k]] %>% left_join(p_x.y_hat[[k]][["1"]], by = c("X_1", "X_2", "X_3")) %>% .$p)/
    (p_y_hat[[k]]*ifelse(is.na(X.test[[k]] %>% left_join(p_x.y_hat[[k]][["1"]], by = c("X_1", "X_2", "X_3")) %>% .$p), 0,
    X.test[[k]] %>% left_join(p_x.y_hat[[k]][["1"]], by = c("X_1", "X_2", "X_3")) %>% .$p) +
    (1 - p_y_hat[[k]])*ifelse(is.na(X.test[[k]] %>% left_join(p_x.y_hat[[k]][["0"]], by = c("X_1", "X_2", "X_3")) %>% .$p), 0,
    X.test[[k]] %>% left_join(p_x.y_hat[[k]][["0"]], by = c("X_1", "X_2", "X_3")) %>% .$p))

  p_y.x_hat.Bayes.test[[k]] <- ifelse(is.na(p_y.x_hat.Bayes.test[[k]]), 0, p_y.x_hat.Bayes.test[[k]])
}

p_y.x_hat.Bayes.test <- suppressMessages(p_y.x_hat.Bayes.test %>% bind_cols() %>% as.matrix()); colnames(p_y.x_hat.Bayes.test) <- paste("p_y_", 1:625, ".x_hat.Bayes", sep = "")

Y_hat.Bayes.test <- ifelse(p_y.x_hat.Bayes.test >= 0.5, 1, 0); colnames(Y_hat.Bayes.test) <- paste("Y_", 1:625, "_hat.Bayes", sep = "")

# Evolve predicted starting grids using Naive Bayes model the appropriate number of time steps to compare resulting grids to provided stopping grids

Y_hat.Bayes.test.evolved <- suppressMessages(bind_cols(test_data[2], Y_hat.Bayes.test) %>% apply(1, function(x){generate_S_stop(x[-1], x[1])}) %>% t())

# Evaluate Bayes model on test data

Bayes_model_MAE <- mean(abs((test_data %>% dplyr::select(-c(1:2)) %>% as.matrix()) - Y_hat.Bayes.test.evolved))

# Train logistic regression models

p_y.x_hat.logistic.test <- list()

# Define logistic regression model for each unique starting grid cell and predictions for test data

for(k in 1:625){
  glm_temp <- train_data.subset[[k]] %>% glm(Y ~ X_1 + X_2 + X_3, data = ., family = "binomial")
  p_y.x_hat.logistic.test[[k]] <- predict(glm_temp, X.test[[k]], type = "response")
}

p_y.x_hat.logistic.test <- suppressMessages(p_y.x_hat.logistic.test %>% bind_cols() %>% as.matrix()); colnames(p_y.x_hat.logistic.test) <- paste("p_y_", 1:625, ".x_hat.logistic", sep = "")

Y_hat.logistic.test <- ifelse(p_y.x_hat.logistic.test >= 0.5, 1, 0); colnames(Y_hat.logistic.test) <- paste("Y_", 1:625, "_hat.logit", sep = "")

# Evolve predicted starting grids using collection of logistic regression models the appropriate number of time steps to compare resulting grids to provided stopping grids

Y_hat.logistic.test.evolved <- suppressMessages(bind_cols(test_data[2], Y_hat.logistic.test) %>% apply(1, function(x){generate_S_stop(x[-1], x[1])}) %>% t())

# Evaluate logistic regression models on test data

logistic_model_MAE <- mean(abs((test_data %>% dplyr::select(-c(1:2)) %>% as.matrix()) - Y_hat.logistic.test.evolved))

# Define and evaluate ensemble of Naive Bayes model and logistic model that uses the average of the conditional probabilities predicted by both models for its decision rule

p_y.x_hat.ensemble.test <- (p_y.x_hat.Bayes.test + p_y.x_hat.logistic.test)/2

Y_hat.ensemble.test <- ifelse(p_y.x_hat.ensemble.test >= 0.5, 1, 0); colnames(Y_hat.ensemble.test) <- paste("Y_", 1:625, "_hat.ensemble", sep = "")

Y_hat.ensemble.test.evolved <- suppressMessages(bind_cols(test_data[2], Y_hat.ensemble.test) %>% apply(1, function(x){generate_S_stop(x[-1], x[1])}) %>% t())

ensemble_model_MAE <- mean(abs((test_data %>% dplyr::select(-c(1:2)) %>% as.matrix()) - Y_hat.ensemble.test.evolved))

rm(i, I, j, J, k, m, glm_temp)
