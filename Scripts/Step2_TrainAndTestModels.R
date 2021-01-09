
# Train Naive Bayes model

p_y_hat <- train_data %>% dplyr::summarize(across(start_0:start_624, mean))

S_stop.N.indices <- list(); p_x.y_hat <- list(); data.stratified <- list()

for(k in 1:625){
  # As was done in function definition for define_aux_matrix(), determine indices (i,j) and (I,J) relating index for each cell and indices for cells in its neighborhood

  j <- ifelse(k %% 25 == 0, 25, k %% 25); i <- ((k - j) + 25)/25
  J <- ifelse((j - 1):(j + 1) %% 25 == 0, 25, (j - 1):(j + 1) %% 25); I <- ifelse((i - 1):(i + 1) %% 25 == 0, 25, (i - 1):(i + 1) %% 25)

  # Define vector of indices for neighboring cells for kth cell

  S_stop.N.indices[[k]] <- plyr::mdply(expand.grid(x = I, y = J), function(x, y){(x - 1)*25 + y}) %>% .[3] %>% pull()

  # Define subset of data corresponding to kth cell that includes the observations of the kth starting grid cell value, delta, the kth stopping grid cell value, and the corresponding
  # neighborhood sum of cells in the stopping grid

  data.stratified[[k]] <- train_data %>% dplyr::select(c(delta, 2 + k, 627 + S_stop.N.indices[[k]]))
  data.stratified[[k]] <- data.stratified[[k]] %>% mutate(X_3 = data.stratified[[k]] %>% as.matrix() %>% .[, 3:11] %>% rowSums()) %>% dplyr::select(c(2, 1, 7, 12)) %>%
    setNames(c("Y", "X_1", "X_2", names(.)[4]))

  p_x.y_hat[[k]] <- list()

  # Given above defined subset of data pertaining to kth cell, define PMFs for X|Y = 0 and X|Y = 1

  for(m in 0:1){
    p_x.y_hat[[k]][[paste(m)]] <- data.stratified[[k]] %>% filter(Y == m) %>% dplyr::select(-Y) %>% group_by(X_1, X_2, X_3) %>% summarize(p = n()/nrow(.), .groups = "drop") %>%
      arrange(X_1, X_2, X_3)
  }
}

# Define predictions for test data

X.test <- list(); p_y.x_hat.test <- list()

for(k in 1:625){

  # Define subset of predictors pertaining to kth starting grid cell

  X.test[[k]] <- test_data %>% dplyr::select(c(delta, 2 + S_stop.N.indices[[k]]))
  X.test[[k]] <- X.test[[k]] %>% mutate(X_3 = X.test[[k]] %>% as.matrix() %>% .[, 2:10] %>% rowSums()) %>% dplyr::select(c(1, 6, 11)) %>% setNames(c("X_1", "X_2", names(.)[3]))

  # Define PMF values for Y|X = x for starting grid cells in test dataset

  p_y.x_hat.test[[k]] <- p_y_hat[[k]]*(X.test[[k]] %>% left_join(p_x.y_hat[[k]][["1"]], by = c("X_1", "X_2", "X_3")) %>% .$p)/
    (p_y_hat[[k]]*(X.test[[k]] %>% left_join(p_x.y_hat[[k]][["1"]], by = c("X_1", "X_2", "X_3")) %>% .$p) +
       (1 - p_y_hat[[k]])*(X.test[[k]] %>% left_join(p_x.y_hat[[k]][["0"]], by = c("X_1", "X_2", "X_3")) %>% .$p))
}

p_y.x_hat.test <- bind_cols(p_y.x_hat.test) %>% as.matrix(); p_y.x_hat.test[is.na(p_y.x_hat.test)] <- 0; colnames(p_y.x_hat.test) <- paste("p_y_", 1:625, ".x_hat", sep = "")

Y_hat.Bayes.test <- ifelse(p_y.x_hat.test >= 0.5, 1, 0); colnames(Y_hat.Bayes.test) <- paste("Y", 1:625, "hat", sep = "_")

# Evolve predicted starting grids using Naive Bayes model the appropriate number of time steps to compare resulting grids to provided stopping grids

Y_hat.Bayes.test.evolved <- bind_cols(test_data[2], Y_hat.Bayes.test) %>% apply(1, function(x){
  generate_S_stop(x[-1], x[1])
})

# Evaluate Bayes model on test data

Bayes_model_MAE <- mean(abs((test_data %>% dplyr::select(-c(1:2)) %>% as.matrix()) - t(Y_hat.Bayes.test.evolved)))

# Train logistic regression models

Y_hat.logistic.test <- list()

# Define logistic regression model for each unique starting grid cell and predictions for test data

for(k in 1:625){
  glm_temp <- data.stratified[[k]] %>% glm(Y ~ X_1 + X_2 + X_3, data = ., family = "binomial")
  p_hat.logistic_temp <- predict(glm_temp, X.test[[k]], type = "response")
  Y_hat.logistic.test[[k]] <- ifelse(p_hat.logistic_temp >= 0.5, 1, 0)
}

Y_hat.logistic.test <- bind_cols(Y_hat.logistic.test) %>% as.matrix(); colnames(Y_hat.logistic.test) <- paste("Y_", 1:625, "_hat.logit", sep = "")

# Evolve predicted starting grids using collection of logistic regression models the appropriate number of time steps to compare resulting grids to provided stopping grids

Y_hat.logistic.test.evolved <- bind_cols(test_data[2], Y_hat.logistic.test) %>% apply(1, function(x){
  generate_S_stop(x[-1], x[1])
})

# Evaluate logistic regression models on test data

logistic_model_MAE <- mean(abs((test_data %>% dplyr::select(-c(1:2)) %>% as.matrix()) - t(Y_hat.logistic.test.evolved)))
