KMEANS_K_DEFAULT = 20
ROBIN_DEFAULT_METHOD = "density"

init_kmeans = function(data, K) {
  ret = kmeans(data, K, nstart = KMEANS_K_DEFAULT)$centers
  return(ret)
}

init_random = function(data, K) {
  ret = data[sample(nrow(data), K), ]
  return(ret)
}

init_robin = function(data, K) {
  ret_rob <- robust_init_density(dist_matrix = dist(data), 
                                 data = data, 
                                 n_clusters = K, 
                                 method = ROBIN_DEFAULT_METHOD)
  ret = data[ret_rob$centers, ]
  return(ret)
}

init_custom = function(data, K) {
  ret = K
  return(ret)
}