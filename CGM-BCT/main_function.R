
# Input parameters
# x_i: a categorical sequence (of class character)
# x_j: a categorical sequence (of class character)
# B: number of bootstrap replicates
# K: number of trees
# depth: maximum tree depth 
# alphabet: the alphabet (0, 1, 2...)

bootstrap_mmd <- function(x_i, x_j, B = 200, K = 5, depth = 10, alphabet = c("0", "1", "2", "3")) {

n_i <- nchar(x_i)
n_j <- nchar(x_j)

# 1. Compute the observed statistic

# Individual i

fit_i <- kBCT(x_i, k = K, depth = depth)
contexts_i <- fit_i$Contexts
posterior_i <- fit_i$Results$posterior

trees_i <- lapply(
  contexts_i,
  function(ct) get_tree_parameters(ct, x_i, alphabet = alphabet)
)

# Individual j

fit_j <- kBCT(x_j, k = K, depth = depth)
contexts_j <- fit_j$Contexts
posterior_j <- fit_j$Results$posterior

trees_j <- lapply(
  contexts_j,
  function(ct) get_tree_parameters(ct, x_j, alphabet = alphabet)
)
  
lambda <- select_lambda_median(trees_i, trees_j)

mmd_observed <- individual_dissimilarity_mmd(
  trees_i,
  posterior_i,
  trees_j,
  posterior_j,
  lambda = lambda
)
  
  
# 2. Fit the BCT pooled model

x_pooled <- paste0(x_i, x_j)

fit_pooled <- kBCT(x_pooled, k = K, depth = depth)
contexts_pooled <- fit_pooled$Contexts
posterior_pooled <- fit_pooled$Results$posterior

trees_pooled <- lapply(
  contexts_pooled,
  function(ct) get_tree_parameters(ct, x_pooled, alphabet = alphabet)
)
  
  
# Bootstrap replicates

mmd_bootstrap <- numeric(B)

w_pooled <- posterior_pooled/sum(posterior_pooled)

# Simulation trees

simulation_trees <- lapply(trees_pooled, function(tree) {

  if (length(tree$contexts) == 1 &&
      tree$contexts == "") {

    return(list(
      iid = TRUE,
      probabilities = tree$probabilities[[1]]
    ))
  }

  list(
    iid = FALSE,
    parameters = setNames(
      tree$probabilities,
      tree$contexts
    )
  )
})

for (b in seq_len(B)) {
  
  # 3. Sample a tree from the pooled posterior
  
  sampled_index <- sample.int(
  length(simulation_trees),
  size = 1,
  prob = w_pooled
)

tree_sampled <- simulation_trees[[sampled_index]]

if (tree_sampled$iid) {

  probs <- tree_sampled$probabilities

  x_i_star <- paste0(
    sample(alphabet, size = n_i,
           replace = TRUE, prob = probs),
    collapse = ""
  )

  x_j_star <- paste0(
    sample(alphabet, size = n_j,
           replace = TRUE, prob = probs),
    collapse = ""
  )

} else {

  x_i_star <- generate_data_fast2(
    tree_sampled$parameters, n_i, alphabet = alphabet
  )

  x_j_star <- generate_data_fast2(
    tree_sampled$parameters, n_j, alphabet = alphabet
  )
}
  
  # 5. Compute the bootstrap statistic
  
  fit_i_star <- kBCT(x_i_star, k = K, depth = depth)
  contexts_i_star <- fit_i_star$Contexts
  posterior_i_star <- fit_i_star$Results$posterior
  
  trees_i_star <- lapply(
    contexts_i_star,
    function(ct) get_tree_parameters(ct, x_i_star ,alphabet = alphabet)
  )
  
  fit_j_star <- kBCT(x_j_star, k = K, depth = depth)
  contexts_j_star <- fit_j_star$Contexts
  posterior_j_star <- fit_j_star$Results$posterior
  
  trees_j_star <- lapply(
    contexts_j_star,
    function(ct) get_tree_parameters(ct, x_j_star, alphabet = alphabet)
  )
  
  mmd_bootstrap[b] <- individual_dissimilarity_mmd(
  trees_i_star,
  posterior_i_star,
  trees_j_star,
  posterior_j_star,
  lambda = lambda
)

# print(b)

}
  

# p-value
  
p_value <- (1 + sum(mmd_bootstrap >= mmd_observed))/(B + 1)
  
return(p_value)
  
}
