library(igraph)

#load and show the text of exercises
load_and_show <- function(file) {
  code <- readLines(file)
  cat("```r\n", paste(code, collapse = "\n"), "\n```")
  eval(parse(text = code),envir = .GlobalEnv)
}

#plot foodweb using flow information
plot_as_flux <- function(g, ...) {
  lay <- layout_with_sugiyama(g)
  coords <- lay$layout
  coords[, 2] <- -coords[, 2]   # invert y-axis (basal bottom → predators top)
  
  # repartir los nodos de cada nivel a lo ancho, manteniendo su orden
  for (y in unique(coords[, 2])) {
    idx <- which(coords[, 2] == y)
    idx <- idx[order(coords[idx, 1])]
    n <- length(idx)
    coords[idx, 1] <- if (n == 1) 0.5 else seq(0, 1, length.out = n)
  }

  plot(g, layout = coords, ...)
}

# cumulative distribution function: P(K >= k)
cdf_get <- function(K) {
  k <- 0:max(K)
  p <- sapply(k, function(kk) mean(K >= kk))
  data.frame(k = k, p = p)
}

plot_foodweb_matrix <- function(A, log = FALSE,
                                show.labels = TRUE,
                                cex.lab = 0.7,
                                ncol = 100) {
  # checks
  stopifnot(is.matrix(A) || is.data.frame(A))
  A <- as.matrix(A)
  
  # replace NA with 0
  A[is.na(A)] <- 0
  
  # optional log-transform (very common for weights)
  if (log) {
    A <- log1p(A)
  }
  
  # make zeros invisible
  A_plot <- A
  A_plot[A_plot == 0] <- NA
  
  # greyscale palette (white = low, black = high)
  cols <- gray.colors(ncol, start = 1, end = 0)
  
  # plot matrix
  image(
    t(A_plot[nrow(A_plot):1, ]),
    col = cols,
    axes = FALSE
  )
  
  # labels
  if (show.labels) {
    axis(
      1,
      at = seq(0, 1, length.out = ncol(A)),
      labels = colnames(A),
      las = 2,
      cex.axis = cex.lab
    )
    axis(
      2,
      at = seq(0, 1, length.out = nrow(A)),
      labels = rev(rownames(A)),
      las = 2,
      cex.axis = cex.lab
    )
  }
  
  invisible(A)
}

get_MusRank <- function(web, mode = "ranking", niterations = 1000, seed = NULL) {
  # web: 2D incidence matrix/data.frame (plants x animals)
  # mode = "ranking" returns scores + sorted names; otherwise returns reordered web

  # Force 2D matrix
  web <- as.matrix(web)
  if (is.null(dim(web)) || length(dim(web)) != 2) {
    stop("web must be a 2D incidence matrix (plants x animals).")
  }

  # Drop empty plant rows (all zeros / NAs treated as 0)
  keep_rows <- rowSums(web != 0, na.rm = TRUE) > 0
  bweb <- web[keep_rows, , drop = FALSE]

  # If only one row survived, keep it as 2D
  if (is.null(dim(bweb)) || length(dim(bweb)) != 2) {
    bweb <- matrix(bweb, nrow = 1)
  }

  # Binarize
  bweb <- ifelse(is.na(bweb), 0, ifelse(bweb != 0, 1, 0))

  # Names (ensure they exist)
  if (is.null(rownames(bweb))) rownames(bweb) <- paste0("P", seq_len(nrow(bweb)))
  if (is.null(colnames(bweb))) colnames(bweb) <- paste0("A", seq_len(ncol(bweb)))

  plants  <- rownames(bweb)
  animals <- colnames(bweb)
  Np <- length(plants)
  Na <- length(animals)

  # Initialise fitness (animals) and complexity (plants)
  fitness    <- stats::runif(Na, 0, 5); names(fitness) <- animals
  complexity <- stats::runif(Np, 0, 5); names(complexity) <- plants

  # Iterations (reset accumulators each iteration)
  for (i in seq_len(niterations)) {
    f2 <- setNames(rep(0, Na), animals)
    c2 <- setNames(rep(0, Np), plants)

    # Animal fitness update: sum of complexities of partner plants
    for (a in animals) {
      p_pals <- plants[bweb[, a] != 0]
      f2[a] <- if (length(p_pals) > 0) sum(complexity[p_pals]) else 0
    }
    fmed <- mean(f2)

    # Plant complexity update: inverse of sum of inverse fitnesses of partner animals
    for (p in plants) {
      a_pals <- animals[bweb[p, ] != 0]
      if (length(a_pals) > 0) {
        denom <- sum(1 / fitness[a_pals])
        c2[p] <- if (is.finite(denom) && denom > 0) 1 / denom else 0
      } else {
        c2[p] <- 0
      }
    }
    cmed <- mean(c2)

    # Normalise (avoid division by zero)
    if (is.finite(fmed) && fmed > 0) fitness <- f2 / fmed
    if (is.finite(cmed) && cmed > 0) complexity <- c2 / cmed
  }

  # Rankings
  sorted_animals <- names(sort(fitness, decreasing = TRUE))
  sorted_plants  <- names(sort(complexity, decreasing = FALSE))

  if (mode == "ranking") {
    return(list(
      animal_fitness = fitness,
      plant_complexity = complexity,
      sorted_animals = sorted_animals,
      sorted_plants = sorted_plants
    ))
  } else {
    # Return reordered matrix (use original 'web' where possible)
    web2 <- web[intersect(sorted_plants, rownames(web)), intersect(sorted_animals, colnames(web)), drop = FALSE]
    return(web2)
  }
}

metric_to_color <- function(metric,
                              ncol = 100,
                              option = "viridis",
                              na_col = "grey90") {
  stopifnot(is.numeric(metric))

  # min–max scaling (ignore NA)
  mmin <- min(metric, na.rm = TRUE)
  mmax <- max(metric, na.rm = TRUE)

  scaled <- (metric - mmin) / (mmax - mmin + 1e-12)

  # viridis palette
  pal <- viridis::viridis(ncol, option = option)

  # map to colors
  cols <- pal[floor(scaled * (ncol - 1)) + 1]

  # handle NA / non-finite
  cols[!is.finite(metric)] <- na_col

  cols
}

membership_to_colors <- function(membership,
                                 palette = viridis::viridis,
                                 na_col = "grey90") {
  stopifnot(is.numeric(membership) || is.factor(membership))

  # handle NA
  cols <- rep(na_col, length(membership))

  ok <- !is.na(membership)
  if (!any(ok)) return(cols)

  # reindex memberships to 1..K (important!)
  mod_ids <- as.integer(factor(membership[ok]))
  n_mod   <- length(unique(mod_ids))

  pal <- palette(n_mod)

  cols[ok] <- pal[mod_ids]
  cols
}

# Samples m vectors randomly on a n-dimensional simplex
simplex_sampling <- function(m, n) {
  r <- list()
  for (j in 1:m) {
    dist <- c(sort(runif(n-1, 0, 1)), 1)
    r[[j]] <- c(dist[1], diff(dist))
  }
  return(r)
}

interaction_matrix <- function(web,gamma_avg,rho,delta){
  SA <- nrow(web)
  SP <- ncol(web)
  alphaA <- matrix(rho,SA,SA) + (1-rho) * diag(rep(1,SA))
  alphaP <- matrix(rho,SP,SP) + (1-rho) * diag(rep(1,SP))
  gammaA <- diag(rowSums(web)^-delta) %*% web
  gammaP <- diag(colSums(web)^-delta) %*% t(web)
  f <- sum(gammaA[web == 1] + gammaP[t(web) == 1] ) / (2 * sum(web==1))
  gammaA <- gamma_avg/f * diag(rowSums(web)^-delta) %*% web
  gammaP <- gamma_avg/f * diag(colSums(web)^-delta) %*% t(web)
  #a[is.nan(a)] = 0
  gammaA[is.nan(gammaA)]=0
  gammaP[is.nan(gammaP)]=0
  alpha <- rbind(cbind(alphaA,-gammaA),cbind(-gammaP,alphaP))
  out <- list(alpha = alpha, alphaA = alphaA, alphaP = alphaP, gammaA = gammaA, gammaP = gammaP)
  return(out)
}

#stability condition gamma_hat (average mutualistic strength at the stability threshold)
#inputs: web = mutualistic network (binary matrix),
#rho = mean field interspecific competition, delta = mutualistic trade-off
#output: gamma_hat = stability condition
gamma_hat <- function(web,rho,delta){
  f_eig <- function(gamma_avg,web,rho,delta){
    alpha <- interaction_matrix(web,gamma_avg,rho,delta)$alpha
    out <- (min(Re(eigen(alpha)$values)))^2
    out
  }
  out <- optimize(f_eig,c(0,1000),web = web, rho = rho, delta = delta)$minimum
  return(out)
}

dd_pair <- function(g1, g2, mode) {
  d1 <- degree_distribution(g1, mode = mode)
  d2 <- degree_distribution(g2, mode = mode)
  n <- max(length(d1), length(d2))
  m <- rbind(c(d1, rep(0, n - length(d1))),
             c(d2, rep(0, n - length(d2))))
  colnames(m) <- 0:(n - 1)
  m
}