

get_tree_parameters <- function(contexts,
                                x,
                                prior = 0.5,
                                alphabet = c("0", "1", "2", "3")) {

  # 1. Convert sequence to a character vector

  if (length(x) == 1) {
    x <- strsplit(as.character(x), split = "")[[1]]
  } else {
    x <- as.character(x)
  }

  # 2. Validate alphabet

  if (!all(x %in% alphabet)) {
    stop("The sequence contains symbols not included in alphabet.")
  }

  T <- length(x)
  A <- length(alphabet)

  # 3. Convert symbols to integer indices

  x_idx <- match(x, alphabet)

  # 4. Compute transition probabilities for one context

  get_theta <- function(context) {

    context <- as.character(context)[1]

    if (identical(context, "")) {

      # Empty context: count all symbols
      counts <- tabulate(x_idx, nbins = A)

    } else {

      # Context stored most recent symbol first
      ctx <- strsplit(context, split = "")[[1]]
      L <- length(ctx)

      counts <- numeric(A)

      if (T > L) {

        # Reverse context into chronological order
        ctx_forward <- rev(ctx)

        # Identify matching positions
        matches <- rep(TRUE, T - L)

        for (j in seq_len(L)) {

          matches <- matches &
            (x[j:(T - L + j - 1)] == ctx_forward[j])

        }

        # Count symbols following matching contexts
        next_symbols <- x_idx[(L + 1):T][matches]

        counts <- tabulate(next_symbols, nbins = A)
      }
    }

    # 5. Dirichlet posterior mean

    probs <- (counts + prior) / sum(counts + prior)

    return(as.numeric(probs))
  }

  # 6. Compute probabilities for all contexts

  probability_list <- lapply(contexts, get_theta)

  # 7. Return results

  return(list(
    contexts = contexts,
    probabilities = probability_list
  ))
}



select_lambda_median <- function(trees_i, trees_j) {

  all_trees <- c(trees_i, trees_j)
  n <- length(all_trees)

  # Preallocate space for all unique pairs
  d_vals <- numeric(n * (n - 1) / 2)
  idx <- 0L

  if (n >= 2) {

    for (a in seq_len(n - 1)) {

      for (b in (a + 1):n) {

        d_ab <- tree_dissimilarity(
          all_trees[[a]],
          all_trees[[b]]
        )

        if (d_ab > 0) {
          idx <- idx + 1L
          d_vals[idx] <- d_ab
        }
      }
    }
  }

  # Keep only strictly positive distances
  d_vals <- d_vals[seq_len(idx)]

  1/median(d_vals)
}



individual_dissimilarity_mmd <- function(trees_i,
                                         posterior_i,
                                         trees_j,
                                         posterior_j,
                                         lambda = NULL) {

  if (is.null(lambda)) {
    lambda <- select_lambda_median(trees_i, trees_j)
  }

  w_i <- posterior_i / sum(posterior_i)
  w_j <- posterior_j / sum(posterior_j)

  n_i <- length(trees_i)
  n_j <- length(trees_j)

  # Initialize kernel matrices
  K_i <- matrix(1, n_i, n_i)
  K_j <- matrix(1, n_j, n_j)
  K_cross <- matrix(0, n_i, n_j)

  # First term: within-group kernel for i
  if (n_i > 1) {
    for (r in seq_len(n_i - 1)) {
      for (l in (r + 1):n_i) {

        value <- exp(
          -lambda * tree_dissimilarity(
            trees_i[[r]], trees_i[[l]]
          )
        )

        K_i[r, l] <- value
        K_i[l, r] <- value
      }
    }
  }

  # Second term: within-group kernel for j
  if (n_j > 1) {
    for (r in seq_len(n_j - 1)) {
      for (l in (r + 1):n_j) {

        value <- exp(
          -lambda * tree_dissimilarity(
            trees_j[[r]], trees_j[[l]]
          )
        )

        K_j[r, l] <- value
        K_j[l, r] <- value
      }
    }
  }

  # Cross term: between-group kernel
  for (r in seq_len(n_i)) {
    for (l in seq_len(n_j)) {

      K_cross[r, l] <- exp(
        -lambda * tree_dissimilarity(
          trees_i[[r]], trees_j[[l]]
        )
      )
    }
  }

  # Compute MMD terms
  term_i <- as.numeric(crossprod(w_i, K_i %*% w_i))
  term_j <- as.numeric(crossprod(w_j, K_j %*% w_j))
  cross <- as.numeric(crossprod(w_i, K_cross %*% w_j))

  # Squared MMD
  D <- term_i + term_j - 2 * cross

  return(D)
}



generate_data_fast2 <- function(ct_theta, N, alphabet = c("0", "1", "2", "3")) {

  contexts <- names(ct_theta)
  depth <- max(nchar(contexts))

  A <- length(alphabet)

  # Root-only model
  if (depth == 0L) {
    return(paste0(
      sample(alphabet, N, replace = TRUE,
             prob = ct_theta[["-"]]),
      collapse = ""
    ))
  }

  # Build lookup table for all possible histories
  n_states <- A^depth

  transition <- matrix(
    NA_real_,
    nrow = n_states,
    ncol = A
  )

  powers <- A^(0:(depth - 1L))

  for (state in 0:(n_states - 1L)) {

    digits <- (state %/% powers) %% A
    history <- paste0(alphabet[digits + 1L], collapse = "")

    # Find the first matching context
    for (j in seq_len(depth)) {

      ct <- substr(history, 1L, j)
      p <- ct_theta[[ct]]

      if (!is.null(p)) {
        transition[state + 1L, ] <- cumsum(p / sum(p))
        break
      }
    }
  }

  if (anyNA(transition)) {
    stop("Some histories have no matching context.")
  }

  # Initial observations
  x <- integer(N + depth)

  x[seq_len(depth)] <- sample.int(
    A, depth, replace = TRUE
  ) - 1L

  # Encode initial history
  state <- sum(
    x[depth:1L] * powers
  )

  # Generate sequence
  for (i in seq.int(depth + 1L, N + depth)) {

    u <- runif(1)

    k <- which(
      u <= transition[state + 1L, ]
    )[1L] - 1L

    x[i] <- k

    # Update history:
    # newest symbol becomes first digit
    state <- k + A * (state %% (n_states / A))
  }

  paste0(alphabet[x[(depth + 1L):(N + depth)] + 1L],
         collapse = "")
}



js_divergence <- function(p, q) {
  
  r <- 0.5 * (p + q)
  
  kl <- function(a, b) {

    idx <- a > 0
    sum(a[idx] * log(a[idx]/b[idx]))

  }
  
  0.5 * kl(p, r) + 0.5 * kl(q, r)

}



tree_dissimilarity <- function(tree1, tree2) {

  contexts1 <- tree1$contexts
  contexts2 <- tree2$contexts

  all_contexts <- union(contexts1, contexts2)
  common_contexts <- intersect(contexts1, contexts2)

  n_all <- length(all_contexts)
  n_common <- length(common_contexts)

  D_max <- log(2)

  # Contribution from non-shared contexts
  total <- (n_all - n_common) * D_max

  # Contributions from shared contexts
  if (n_common > 0) {

    idx1 <- match(common_contexts, contexts1)
    idx2 <- match(common_contexts, contexts2)

    for (j in seq_len(n_common)) {

      total <- total + js_divergence(
        tree1$probabilities[[idx1[j]]],
        tree2$probabilities[[idx2[j]]]
      )

    }
  }

  return(total / n_all)
}



