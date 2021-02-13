
# Define names of packages to load for scripts. Not all are actually necessary

Packages <- c("bigmemory", "bigstatsr", "broom", "caret", "compiler", "caTools", "data.table", "doParallel", "doSNOW", "dplyr", "dslabs", "e1071", "fastAdaboost", "foreach", "formatR",
  "future", "gam", "genefilter", "ggplot2", "ggrepel", "gridExtra", "HistData", "ipc", "kernlab", "knitr", "Lahman", "lpSolve", "lubridate", "MASS", "matrixStats", "mvtnorm",
  "naivebayes", "parallel", "pdftools", "promises", "purrr", "randomForest", "ranger", "Rborist", "RColorBrewer", "recommenderlab", "recosystem", "reshape2", "rlist", "ROSE",
  "rpart", "rpart.plot", "rtweet", "rvest", "scales", "snow", "stringr", "svMisc", "svSocket", "textdata", "tibble", "tidyr", "tidytext", "tidyverse", "tree", "zoo")

# Download and install any packages not already installed and then load them

for(p in Packages){
  if(!require(p, character.only = TRUE)){install.packages(p, character.only = TRUE, repos = "http://cran.us.r-project.org")}
  library(p, character.only = TRUE)
}

# Set the number of logical processors to be used for parallel computing as well as the corresponding cluster of child processes, and register the back-end for the cluster

logical_CPUs <- detectCores(logical = TRUE)

reverse_GoL_cluster <- makeCluster(logical_CPUs); registerDoParallel(reverse_GoL_cluster)

# Define a function that takes a game grid in vector form and creates a matrix consisting of 0s and 1s, the locations of which being based on both the relationship between the
# indices of cells when game grids are represented as vectors and the indices of cells when game grids are represented as matrices as well as the relationship between the indices of
# cells when game grids are represented as matrices and the indices of their neighboring cells when game grids are represented as matrices

create_aux_matrix <- function(S.p1){
  S.length <- length(S.p1); S.mtrx.dim <- sqrt(S.length)
  if(S.mtrx.dim %% 1 != 0){
    print("Provided vector can't be expressed as a square matrix!")
  }else{
    M_aux <- matrix(nrow = S.length, ncol = S.length)
    for(k in 1:S.length){
      # Define index pair (i,j) of grid in matrix form corresponding to index k of grid in vector form
      j <- ifelse(k %% S.mtrx.dim == 0, S.mtrx.dim, k %% S.mtrx.dim)
      i <- ((k - j) + S.mtrx.dim)/S.mtrx.dim
      # Define indices for 9 cells in neighborhood of each cell, including the cell itself, in grid in matrix form
      for(J in ifelse((j - 1):(j + 1) %% S.mtrx.dim == 0, S.mtrx.dim, (j - 1):(j + 1) %% S.mtrx.dim)){
        for(I in ifelse((i - 1):(i + 1) %% S.mtrx.dim == 0, S.mtrx.dim, (i - 1):(i + 1) %% S.mtrx.dim)){
          # Define indices K of grid in vector form corresponding to index pairs (I,J) of grid in matrix form
          K <- (I - 1)*S.mtrx.dim + J
          M_aux[k, K] <- 1
        }
      }
    }
    # Define a value of 0 to elements not assigned a value of 1
    M_aux[is.na(M_aux)] <- 0
    return(M_aux)
  }
}

# Define a function that takes a game grid in vector form, evolves the grid forward one time step, and returns the resulting grid in vector form

evolve_S <- function(S.p2){
  N <- create_aux_matrix(S.p2) %*% S.p2
  ifelse(S.p2 == 1, ifelse(N %in% c(3, 4), 1, 0), ifelse(N == 3, 1, 0))
}

# Define a function that takes a game grid in vector form, evolves the game grid forward the specified number of time steps, and returns the resulting grid in vector form

generate_S_stop <- function(S_start.p, delta.p1){
  S_temp <- S_start.p
  for(j in 1:delta.p1){
    S_temp <- evolve_S(S_temp)
  }
  return(S_temp)
}

# Define a function that takes a game grid in vector form, that game grid's corresponding neighborhood sum matrix in vector form and a second game grid in vector form and checks the
# relationships between the three objects, returning TRUE if there are no inconsistencies between any of the objects' elements (assuming the objects follow the rules of Conway's Game
# of Life)

check_S_parent <- function(S_parent.p, S_child.p, N.p){
  all((S_child.p == 0 | (N.p %in% c(3, 4) & (N.p != 4 | S_parent.p == 1))) & (S_child.p == 1 | (N.p != 3 & (N.p != 4 | S_parent.p == 0))) &
        ((N.p %in% 1:9 | S_parent.p == 0) & (N.p %in% 0:8 | S_parent.p == 1) & S_parent.p %in% c(0, 1)))
}

# Define a function that takes a game grid in vector form, recursively solves it backwards the specified number of times, and then returns the first resulting grid ancestor that it
# generates as a solution in vector form

solve_S_ancestor <- function(i.p2, S_descendant.p, t.p){
  # Set seed for the sake of reproducibility
  set.seed(i.p2, sample.kind = "Rounding")
  # Define permutation of indices of ancestor grid's cells that are non-trivial and determine the living cells in S_stop
  indices <- M[which(S_descendant.p == 1),] %>% apply(1, function(x){which(x == 1)}) %>% as.vector() %>% unique() %>% sample()
  # Define the set of candidate values for each cell in parent grid, based on whether corresponding cells in S_stop are living or dead
  S.range <- lapply(as.list(1:length(S_descendant.p)), function(s){
    if(s %in% indices){
      c(0, 1)
    }else{
      0
    }
  })
  # Initialize index to 1 for all cells in parent grid, so that, to start, the first element of each set in S.range is assigned to the value of its corresponding cell in parent grid
  S.index <- rep(1, times = length(S_descendant.p))
  # Rather than using for() loop, use repeat() loop, as number of iterations of loop is unknown at execution
  repeat{
    # If complete_ind[1] has been reassigned a value of 1 by the start of each iteration of the repeat() loop, exit the loop, as this indicates a solution has already been found by
    # one of the other child processes in the cluster and there is therefore no need for this particular process to continue searching for a solution
    if(complete_ind[1] == 1){return(NULL)}
    # Assign values to parent grid for current iteration
    S_parent <- mapply(function(x, y){x[[y]]}, S.range, S.index)
    # Define vector of neighborhood sums corresponding to resulting parent grid
    N <- M %*% S_parent
    if(check_S_parent(S_parent, S_descendant.p, N) == TRUE){
      # Record ancestor grid in log of reversal process, then either assign indicator variable a value of "Y" so that function is exited or move on to solving previous ancestor grid
      game_log[[t.p]] <<- list(name = paste("S_", t.p - 1, sep = ""), value = S_parent)
      if(t.p == 1){
        solution_found <<- "Y"
      }else{
        solve_S_ancestor(i.p2, S_parent, t.p - 1)
      }
    }
    m <- match(TRUE, S.index[indices] < 2)
    # If either a solution is found or there are no more possible values to try for parent grid, exit the current repeat() loop, otherwise define the indices to be used to determine
    # the next candidate parent grid and continue
    if(solution_found == "Y" | is.na(m)){
      break
    }else{
      S.index[indices][index(S.index[indices]) < m] <- 1
      S.index[indices][m] <- S.index[indices][m] + 1
    }
  }
}

# Define a function that takes a game grid in vector form then first initializes a few objects before calling on solve_S_ancestor() to either return the first resulting grid ancestor
# sequence generated in vector form or return NULL if no such sequence exists

solve_S_start <- function(i.p1, S_stop.p1, delta.p2){
  solution_found <<- "N"
  game_log <<- list()
  M <<- create_aux_matrix(S_stop.p1)
  solve_S_ancestor(i.p1, S_stop.p1, delta.p2)
  if(solution_found == "Y"){complete_ind[1] <- 1}
  return(game_log)
}

# Define a function that simultaneously executes the solve_S_start() function in parallel on each logical processor

parallelize_solve_S_start  <- function(S_stop.p2, delta.p3){
  foreach(i = 1:logical_CPUs, .export = c("solve_S_start", "complete_ind", "create_aux_matrix", "solve_S_ancestor", "%>%", "check_S_parent", "index")) %dopar% {
    solve_S_start(i, S_stop.p2, delta.p3)
  }
}

# Define the final outside function that either returns the first solution found between the different child processes executed in parallel or returns a message indicating no solution
# exists

aggregated_solve_S_start <- function(S_stop.p3, delta.p4){
  # Initialize the 1 x 1 matrix, complete_ind, the value of which will be shared between all child processes executing solve_S_start() in parallel
  complete_ind <<- FBM(1, 1, init = 0)
  # Clean up the list of outputs returned by solve_S_start() by removing all NULL elements (including those that are nested)
  x <- parallelize_solve_S_start(S_stop.p3, delta.p4) %>% list.clean(., recursive = TRUE)
  if(complete_ind[1] == 1){
    # If complete_ind[1] == 1, then at least one of the remaining elements of x will consist of a complete sequence of ancestor grids, and the below statement will return the first such
    # element of x for which this is true
    return(x[[match(delta.p4, lapply(x, function(y){length(y)}) %>% cbind())]])
  }else{
    # If complete_ind[1] != 1, then none of the child processes executing solve_S_start() in parallel found a solution, so print the below message
    print("No solution exists!")
  }
}

# Unzip and import train and test datasets into R

data.dir.path <- paste(getwd(), "/Data", sep = "")

zipped_data.filename <- "conways-reverse-game-of-life-2020.zip"; zipped_data.path <- file.path(data.dir.path, zipped_data.filename)

unzip(zipfile = zipped_data.path, exdir = data.dir.path)

train_data.filename <- "train.csv"; train_data.path <- file.path(data.dir.path, train_data.filename)

test_data.filename <- "test.csv"; test_data.path <- file.path(data.dir.path, test_data.filename)

train_data <- read_csv(train_data.path); test_data <- read_csv(test_data.path)

rm(Packages, p, data.dir.path, test_data.filename, test_data.path, train_data.filename, train_data.path, zipped_data.filename, zipped_data.path)

