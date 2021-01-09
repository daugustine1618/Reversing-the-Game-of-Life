
# Example 1:

S_stop <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

S_stop.dim <- sqrt(length(S_stop))

S_stop.mtrx <- matrix(S_stop, c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_stop.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_stop.mtrx[S_stop.dim:1, ]))

delta <- 1

system.time(solve_S_start(S_stop, delta)); print(iter_ct)

S_start.mtrx <- matrix(game_log[[1]][["value"]], c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_start.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_start.mtrx[S_stop.dim:1, ]))

# Example 2:

S_stop <- c(1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0)

S_stop.dim <- sqrt(length(S_stop))

S_stop.mtrx <- matrix(S_stop, c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_stop.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_stop.mtrx[S_stop.dim:1, ]))

delta <- 5

system.time(solve_S_start(S_stop, delta)); print(iter_ct)

for(j in 1:delta){
  S_start.mtrx <- matrix(game_log[[j]][["value"]], c(S_stop.dim, S_stop.dim), byrow = TRUE)
  S_start.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_start.mtrx[S_stop.dim:1, ]))
}

rm(game_log, M, S_start.mtrx, S_stop.mtrx, delta, iter_ct, j, raw_data.dir.path, S_start.image, S_stop, S_stop.dim, S_stop.image, solution_found, test_data.filename,
   test_data.path, train_data.filename, train_data.path)
