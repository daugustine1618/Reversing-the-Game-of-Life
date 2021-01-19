
# Example 1 (5x5 grid with delta = 1):

S_stop <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

S_stop.dim <- sqrt(length(S_stop))

S_stop.mtrx <- matrix(S_stop, c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_stop.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_stop.mtrx[S_stop.dim:1, ]))

delta <- 1

system.time(S_start_seq <- aggregated_solve_S_start(S_stop, delta))

S_start.mtrx <- matrix(S_start_seq[[1]][["value"]], c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_start.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_start.mtrx[S_stop.dim:1, ]))

# Example 2 (4x4 grid with delta = 5):

S_stop <- c(1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0)

S_stop.dim <- sqrt(length(S_stop))

S_stop.mtrx <- matrix(S_stop, c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_stop.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_stop.mtrx[S_stop.dim:1, ]))

delta <- 5

system.time(S_start_seq <- aggregated_solve_S_start(S_stop, delta))

for(j in 1:delta){
  S_start.mtrx <- matrix(S_start_seq[[j]][["value"]], c(S_stop.dim, S_stop.dim), byrow = TRUE)
  S_start.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_start.mtrx[S_stop.dim:1, ]))
}

# Example 3 (25x25 grid with delta = 3):

S_stop <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

S_stop.dim <- sqrt(length(S_stop))

S_stop.mtrx <- matrix(S_stop, c(S_stop.dim, S_stop.dim), byrow = TRUE)

S_stop.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_stop.mtrx[S_stop.dim:1, ]))

delta <- 3

system.time(S_start_seq <- aggregated_solve_S_start(S_stop, delta))

for(j in 1:delta){
  S_start.mtrx <- matrix(S_start_seq[[j]][["value"]], c(S_stop.dim, S_stop.dim), byrow = TRUE)
  S_start.image <- image(1:S_stop.dim, 1:S_stop.dim, t(S_start.mtrx[S_stop.dim:1, ]))
}

stopCluster(reverse_GoL_cluster)

rm(S_start.mtrx, S_stop.mtrx, delta, j, S_start.image, S_stop, S_stop.dim, S_stop.image, reverse_GoL_cluster)
