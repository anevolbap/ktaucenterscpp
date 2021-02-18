KMEANS_K_DEFAULT = 20

initialize_kmeans = function(X, K) {
  ret = kmeans(X, K, nstart = KMEANS_K_DEFAULT)$centers
  return(ret)
}
  
