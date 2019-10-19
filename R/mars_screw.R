#' Intensity and saturation values of a picture from mars.
#'
#' A dataset containing the Intensity and Saturation values of a picture from Mars
#' taken from Rover Curiosity.
#' @format A list containing information about  pixels of a picture form mars
#' mainly containing red sand and  metal form Rover itself. List include
#'
#' \itemize{
#' \item SI_matrix: A matrix with 5063 rows and 128 columns.
#' Elements 1 to 64 of each row indicate the Saturation values of pixels in a square cell 8 x 8
#' whereas elements 65 to 128 of each row indicate the cell's Intensity  values.
#' \item geographic_matrix: An integer matrix of dimension 5063 x 2, each row indicates
#'  each square cell's locations (x-axis y-axis) at the picture.
#'  \item screw_index: the index corresponding to the screw observation (screw_index=4180)
#'}
#'
#'@examples
#'
#' # upload matrix 
#' A <- mars_screw$SI_matrix;
#' B <- mars_screw$geographic_matrix
#' screw_index <- mars_screw$screw_index
#' ## looking for the brightest cells
#' maximun_at_each_cell<- apply(A[ ,1:64], 1, max)
#' ten_brightest_cells <-order(maximun_at_each_cell, decreasing=TRUE)[1:10]
#'
#' ## plot locations where the ten brightest cells are.
#' plot(B, pch=19 )
#' points(B[ten_brightest_cells, ], pch=19,col="yellow" )
#'
#' ## plot locations where the screw are.
#' points(B[screw_index, ], pch=19,col="blue" )

#'
#' @source \url{https://www.nasa.gov/mission_pages/msl/multimedia/pia16225.html}
"mars_screw"

