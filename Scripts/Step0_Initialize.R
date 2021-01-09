
# Define names of packages to load for scripts. Not all are actually necessary

Packages <- c("bigmemory", "broom", "caret", "compiler", "caTools", "data.table", "doParallel", "doSNOW", "dplyr", "dslabs", "e1071", "fastAdaboost", "foreach", "formatR", "gam",
              "genefilter", "ggplot2", "ggrepel", "gridExtra", "HistData", "kernlab", "knitr", "Lahman", "lpSolve", "lubridate", "MASS", "matrixStats", "mvtnorm", "naivebayes",
              "parallel", "pdftools", "purrr", "randomForest", "ranger", "Rborist", "RColorBrewer", "recommenderlab", "recosystem", "reshape2", "ROSE", "rpart", "rpart.plot", "rtweet",
              "rvest", "scales", "snow", "stringr", "svMisc", "textdata", "tibble", "tidyr", "tidytext", "tidyverse", "tree", "zoo")

# Download and install any packages not already installed and then load them

for(p in Packages){
  if(!require(p, character.only = TRUE)){install.packages(p, character.only = TRUE, repos = "http://cran.us.r-project.org")}
  library(p, character.only = TRUE)
}

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

solve_S_ancestor <- function(S_descendent.p, t.p){
  # Define row indices of auxiliary matrix M that correspond to those cells in S_stop that are living
  indices <- M[which(S_stop == 1),] %>% apply(1, function(x){which(x == 1)}) %>% as.vector() %>% unique() %>% sort()
  # Define the set of candidate values for each cell in parent grid, based on whether corresponding cells in S_stop are living or dead
  S.range <- lapply(as.list(1:length(S_descendent.p)), function(s){
    if(s %in% indices){
      c(0, 1)
    }else{
      0
    }
  })
  # Initialize index to 1 for all cells in parent grid, so that, to start, the first element of each set in S.range is assigned to the value of its corresponding cell in parent grid
  S.index <- rep(1, times = length(S_descendent.p))
  # Rather than using for() loop, use repeat() loop, as number of iterations of loop is unknown at execution
  repeat{
    iter_ct <<- iter_ct + 1
    # Assign values to parent grid for current iteration
    S_parent <- mapply(function(x, y){x[[y]]}, S.range, S.index)
    # Define vector of neighborhood sums corresponding to resulting parent grid
    N <- M %*% S_parent
    if(check_S_parent(S_parent, S_descendent.p, N) == TRUE){
      # Record ancestor grid in log of reversal process, then either assign indicator variable a value of "Y" so that function is exited or move on to solving previous ancestor grid
      game_log[[t.p]] <<- list(name = paste("S_", t.p - 1, sep = ""), value = S_parent)
      if(t.p == 1){
        solution_found <<- "Y"
      }else{
        solve_S_ancestor(S_parent, t.p - 1)
      }
    }
    m <- match(TRUE, S.index[indices] < 2)
    # If either a solution is found or there are no more possible values to try for parent grid, exit the current repeat() loop, otherwise define the indices to be used to determine
    # next candidate parent grid and continue
    if(solution_found == "Y" | is.na(m)){
      break
    }else{
      S.index[indices][index(S.index[indices]) < m] <- 1
      S.index[indices][m] <- S.index[indices][m] + 1
    }
  }
}

# Define a function that takes a game grid in vector form then first initializes a few objects before calling on solve_S_ancestor() to either ultimately return the first resulting grid
# ancestor generated in vector form or a message stating no grid ancestor exists for the input game grid and number of generations

solve_S_start <- function(S_stop.p, delta.p2){
  stopifnot(sqrt(length(S_stop.p)) %% 1 == 0)
  iter_ct <<- 1
  game_log <<- list(); solution_found <<- "N"
  M <<- create_aux_matrix(S_stop.p)
  solve_S_ancestor(S_stop.p, delta.p2)
  if(solution_found == "N"){print("No solution exists!")}
}

# Import train and test datasets into R

raw_data.dir.path <- paste(getwd(), "/Data", sep = "")

train_data.filename <- "train.csv"; train_data.path <- file.path(raw_data.dir.path, train_data.filename)

test_data.filename <- "test.csv"; test_data.path <- file.path(raw_data.dir.path, test_data.filename)

train_data <- read_csv(train_data.path); test_data <- read_csv(test_data.path)

rm(Packages, p)
