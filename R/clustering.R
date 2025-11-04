#output directory
dir.create("project/outputs/figures", recursive = TRUE, showWarnings = FALSE)

library(mvtnorm)
library(cluster)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(scales)
library(plotly)
library(htmlwidgets)
library(tibble)
library(reticulate)

py_config()
main <- import_main()


#===============TASK1=================

#=============generate_hypercube_clusters function
generate_hypercube_clusters <- function(n,k,side_length,noise_sd=1.0){
  
  L <- side_length
  datapoints <- as.data.frame(matrix(0,nrow=n*k,ncol=n+1))
  colnames(datapoints) <- c(paste0("axis",seq(1:n)),"cluster")
  cov_mat <- (noise_sd^2)*diag(n)
  
  for (i in 1:n){
    rows <- ((i - 1) * k + 1):(i * k)
    mean_vec <- rep(0, n)
    mean_vec[i] <- L
    datapoints[rows, 1:n] <- mvtnorm::rmvnorm(k, mean = mean_vec, sigma = cov_mat)
    datapoints[rows, n + 1] <- i
  }
  datapoints$cluster <- factor(datapoints$cluster, levels = seq_len(n))
  datapoints
}
#n: The number of dimensions and clusters.
#k: The number of points per cluster.
#side-length: The distance of cluster centers from the origin along the axes.
#noise_sd: The standard deviation of the Guassian noise, controlling cluster spread.


#=============Simulation:Clustering data

#hyperparameters
grid_n <- c(6,5,4,3,2)
grid_L <- c(seq(1,3,0.5),4:10)
k <- 100

# run clusGap and return k-hat using the globalSEmax rule
estimate_k_clusgap <- function(x, K.max, nstart = 20, iter.max = 50, B = 100) {
  g <- clusGap(
    x,
    FUNcluster = kmeans,
    nstart = nstart,
    iter.max = iter.max,
    K.max = K.max,
    B = B
  )
  gap <- g$Tab[, "gap"]
  se  <- g$Tab[, "SE.sim"]
  k_hat <- as.integer(cluster::maxSE(gap, se, method = "globalSEmax"))
  # Safety: ensure k_hat within [1, K.max]
  k_hat <- max(1L, min(as.integer(k_hat), nrow(g$Tab)))
  return(k_hat)
}

#Simulation loop
set.seed(2025)
results <- map_dfr(grid_n, function(n) {
  map_dfr(grid_L, function(L) {
    # Generate data (n dims, n clusters, k points/cluster)
    dat <- generate_hypercube_clusters(n, k, side_length = L, noise_sd = 1.0)
    X <- as.matrix(dat[, 1:n, drop = FALSE])
    
    K.max <- max(n + 2, 8)
    
    k_hat <- estimate_k_clusgap(
      X,
      K.max = K.max,
      nstart = 20,
      iter.max = 50,
      B = 100  
    )
    
    tibble(
      n = n,
      side_length = L,
      k_hat = k_hat,
      true_k = n
    )
  })
})
    
#Visualization
label_breaks <- c(1, 1.5, 2, 2.5, 3, 4, 6, 8, 10)
p <- results %>%
  mutate(n = factor(n, levels = sort(unique(n), decreasing = TRUE))) %>%
  ggplot(aes(x = side_length, y = k_hat)) +
  geom_hline(aes(yintercept = true_k), linetype = "dashed") +
  geom_line(size = 0.9) +
  geom_point(size = 1) +
  scale_x_continuous(
    breaks = label_breaks,
    minor_breaks = grid_L,
    labels = scales::label_number(accuracy = 0.5)
  )+
  coord_cartesian(ylim = c(1, max(results$true_k) + 2)) +
  facet_wrap(~ n, nrow = 1, labeller = labeller(n = function(v) paste0("Dimension n = ", v))) +
  labs(
    title = "Estimated Number of Clusters vs. Side Length (Gap Statistic, globalSEmax)",
    subtitle = "Dashed line = true number of clusters (k = n). Lower side_length increases overlap and reduces estimated k.",
    x = "Side length (L)",
    y = "Estimated number of clusters (k̂)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p)
ggsave("project/outputs/figures/kmeans_cluster.png", width = 15, height = 5, dpi = 300)


#===============TASK2=================

# --- 1) Data generator ---
generate_shell_clusters <- function(n_shells, k_per_shell, max_radius, noise_sd = 0.1) {
  radii <- seq(from = max_radius / n_shells, to = max_radius, length.out = n_shells)
  pts <- vector("list", n_shells)
  for (i in seq_len(n_shells)) {
    r <- radii[i]
    theta <- runif(k_per_shell, 0, 2*pi)
    phi   <- acos(2*runif(k_per_shell) - 1)
    noisy_r <- pmax(r + rnorm(k_per_shell, 0, noise_sd), 0)
    x <- noisy_r * sin(phi) * cos(theta)
    y <- noisy_r * sin(phi) * sin(theta)
    z <- noisy_r * cos(phi)
    pts[[i]] <- data.frame(x, y, z, shell = factor(i))
  }
  do.call(rbind, pts)
}

# --- 2) Interactive 3D scatter ---
dat <- generate_shell_clusters(
  n_shells     = 4,
  k_per_shell  = 600,
  max_radius   = 3,
  noise_sd     = 0.1
)

dat$r <- sqrt(dat$x^2 + dat$y^2 + dat$z^2)

fig <- plot_ly(
  dat,
  x = ~x, y = ~y, z = ~z,
  type  = "scatter3d",
  mode  = "markers",
  color = ~shell,               
  marker = list(size = 2, opacity = 0.7),
  text   = ~paste0("shell = ", shell),
  hovertemplate = paste(
    "%{text}",
    "<br>x=%{x:.2f}",
    "<br>y=%{y:.2f}",
    "<br>z=%{z:.2f}",
    "<extra></extra>"
  )
) |>
  layout(
    title = "Concentric Shell Clusters (interactive)",
    scene = list(
      aspectmode = "cube",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    legend = list(title = list(text = "Shell"))
  )

fig
htmlwidgets::saveWidget(fig, "project/outputs/figures/shell_clusters_interative.html")

# --- 3) Spectral clustering ---
# clusGap-compatible wrapper: function(data, k) -> list(cluster = <int vector>)
spectral_gap_wrapper <- function(d_thresh = 1, random_state = 0L, n_init = 10L) {
  cache <- new.env(parent = emptyenv())
  function(data, k) {
    X <- as.matrix(data)
    
    if (is.null(cache$key) || !identical(cache$key, X)) {
      main$X_in <- X
      main$dth  <- as.numeric(d_thresh)
      
      py_run_string("
import numpy as np

def _pairwise_dists(X):
    G = X @ X.T
    sq = np.clip(np.diag(G), 0.0, None)
    D2 = sq[:, None] + sq[None, :] - 2.0*G
    np.maximum(D2, 0.0, out=D2)
    return np.sqrt(D2, dtype=X.dtype)

D  = _pairwise_dists(X_in)
A  = (D <= dth).astype(float)
np.fill_diagonal(A, 0.0)

deg = A.sum(axis=1)
with np.errstate(divide='ignore'):
    dinvsqrt = np.where(deg > 0, 1.0/np.sqrt(deg), 0.0)
Dinv = np.diag(dinvsqrt)
Lsym = np.eye(A.shape[0]) - (Dinv @ A @ Dinv)

w, U = np.linalg.eigh(Lsym)
idx = np.argsort(w)
U_sorted = U[:, idx]
")
      cache$U_sorted <- main$U_sorted
      cache$key <- X
    }
    
    # take first k eigenvectors, row-normalize
    V  <- cache$U_sorted[, 1:k, drop = FALSE]
    rn <- sqrt(rowSums(V^2)); rn[rn == 0] <- 1
    Y  <- V / rn
    
    # If there are fewer than k unique rows (degenerate bootstrap), add tiny jitter
    if (nrow(unique(Y)) < k) {
      Y <- Y + matrix(rnorm(length(Y), sd = 1e-8), nrow = nrow(Y))
    }
    
    set.seed(as.integer(random_state))
    km <- kmeans(Y, centers = k,
                 nstart = max(50, n_init), iter.max = 200,
                 algorithm = "Lloyd")   # still base kmeans, but avoids the empty-cluster error
    list(cluster = as.integer(km$cluster))
  }
}

# ---- 4) Simulation----
set.seed(44)
n_shells    <- 4
k_per_shell <- 100
noise_sd    <- 0.1
d_threshold <- 1

radii_grid <- seq(10, 0, by = -0.5)
Kmax <- 8
B    <- 30  

FUN <- spectral_gap_wrapper(d_thresh = d_threshold)

pick_k <- function(cg) {
  ks <- 1:nrow(cg$Tab)
  ks[maxSE(cg$Tab[,"gap"], cg$Tab[,"SE.sim"], method = "Tibs2001SEmax")]
}

res <- lapply(radii_grid, function(Rmax) {
  dat <- generate_shell_clusters(n_shells, k_per_shell, Rmax, noise_sd)
  X   <- as.matrix(dat[, c("x","y","z")])
  cg  <- clusGap(X, FUNcluster = FUN, K.max = Kmax, B = B, verbose = interactive())
  data.frame(max_radius = Rmax, k_hat = pick_k(cg))
})
df_res <- bind_rows(res)

# ---- PLOT ----
ggplot(df_res, aes(max_radius, k_hat)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_x_continuous(breaks = radii_grid) +
  labs(
    title = "Estimated clusters vs. max_radius (Spectral + clusGap)",
    subtitle = paste0("n_shells=", n_shells, ", k_per_shell=", k_per_shell,
                      ", noise_sd=", noise_sd, ", d_threshold=", d_threshold,
                      ", K.max=", Kmax, ", B=", B),
    x = "max_radius", y = "Estimated K"
  ) +
  theme_minimal(base_size = 13)

df_res

ggsave("project/outputs/figures/spectral_cluster_1.png", width = 10, height = 8, dpi = 300)


#d_threshold <- 0.8
d_threshold <- 0.8

FUN <- spectral_gap_wrapper(d_thresh = d_threshold)
res <- lapply(radii_grid, function(Rmax) {
  dat <- generate_shell_clusters(n_shells, k_per_shell, Rmax, noise_sd)
  X   <- as.matrix(dat[, c("x","y","z")])
  cg  <- clusGap(X, FUNcluster = FUN, K.max = Kmax, B = B, verbose = interactive())
  data.frame(max_radius = Rmax, k_hat = pick_k(cg))
})
df_res <- bind_rows(res)

# ---- PLOT ----
ggplot(df_res, aes(max_radius, k_hat)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_x_continuous(breaks = radii_grid) +
  labs(
    title = "Estimated clusters vs. max_radius (Spectral + clusGap)",
    subtitle = paste0("n_shells=", n_shells, ", k_per_shell=", k_per_shell,
                      ", noise_sd=", noise_sd, ", d_threshold=", d_threshold,
                      ", K.max=", Kmax, ", B=", B),
    x = "max_radius", y = "Estimated K"
  ) +
  theme_minimal(base_size = 13)

df_res

ggsave("project/outputs/figures/spectral_cluster_0.8.png", width = 10, height = 8, dpi = 300)

#d_threshold <- 1.2
d_threshold <- 1.2

FUN <- spectral_gap_wrapper(d_thresh = d_threshold)
res <- lapply(radii_grid, function(Rmax) {
  dat <- generate_shell_clusters(n_shells, k_per_shell, Rmax, noise_sd)
  X   <- as.matrix(dat[, c("x","y","z")])
  cg  <- clusGap(X, FUNcluster = FUN, K.max = Kmax, B = B, verbose = interactive())
  data.frame(max_radius = Rmax, k_hat = pick_k(cg))
})
df_res <- bind_rows(res)

# ---- PLOT ----
ggplot(df_res, aes(max_radius, k_hat)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_x_continuous(breaks = radii_grid) +
  labs(
    title = "Estimated clusters vs. max_radius (Spectral + clusGap)",
    subtitle = paste0("n_shells=", n_shells, ", k_per_shell=", k_per_shell,
                      ", noise_sd=", noise_sd, ", d_threshold=", d_threshold,
                      ", K.max=", Kmax, ", B=", B),
    x = "max_radius", y = "Estimated K"
  ) +
  theme_minimal(base_size = 13)

df_res

ggsave("project/outputs/figures/spectral_cluster_1.2.png", width = 10, height = 8, dpi = 300)
