#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//' rhoAux function
//'
//' Auxiliary function to compute rho and psi
// [[Rcpp::export]]
NumericVector rhoAux(NumericVector x, double cc) {
    return (1.792 -
            0.972 * pow(x, 2) / pow(cc, 2) +
            0.432 * pow(x, 4) / pow(cc, 4) -
            0.052 * pow(x, 6) / pow(cc, 6) +
            0.002 * pow(x, 8) / pow(cc, 8)) / 3.25;
}

//'rhoOpt function
//'
//' An implementation of quasi optimal rho functions following reference [1]
//'
//'@param x real number.
//'@param cc tunning constant.
//'@return rho(\code{x}/\code{c}\code{c})
//'@examples val = rhoOpt(x = 0.5, cc = 1)
//'@references [1] Salibian-Barrera, M., Willems, G., & Zamar,
//'   R. (2008). The fast-tau estimator for regression. Journal of
//'   Computational and Graphical Statistics, 17(3), 659-682.
//'
//'@export
// [[Rcpp::export]]
NumericVector rhoOpt(NumericVector x, double cc) {
    NumericVector out(x.size());

    out = pow(x, 2) / 2 / (3.25 * pow(cc, 2));
    out[abs(x) > 2*cc] = rhoAux(x[abs(x) > 2 * cc], cc);
    out[abs(x) > 3*cc] = 1;

    return out;
}

//'psiAux function
//'
//' An implementation of quasi optimal rho functions following reference [1]
//'
//'@param x real number.
//'@param cc tunning constant.
//'@return rho(\code{x}/\code{c}\code{c})
//'@examples val = rhoOpt(x = 0.5, cc = 1)
//'@references [1] Salibian-Barrera, M., Willems, G., & Zamar,
//'   R. (2008). The fast-tau estimator for regression. Journal of
//'   Computational and Graphical Statistics, 17(3), 659-682.
//'
//'@export
//[[Rcpp::export]]
NumericVector psiAux(NumericVector x, double cc) {
    return (-1.944 * x / pow(cc, 2) +
            1.728 * pow(x,3) / pow(cc, 4) -
            0.312 * pow(x,5) / pow(cc, 6) +
            0.016 * pow(x,7) / pow(cc, 8)) / 3.25;
}

//' psiOpt
//' psi function is the derivative of optimal rho function (see \code{\link{rhoOpt}})
//' @param x A real number.
//' @param cc a tunning constant.
//' @return rho''(\code{x}/ \code{c}\code{c})
//' @examples
//' psiOpt(x=0.5,cc=1)
//' @export
// [[Rcpp::export]]
NumericVector psiOpt(NumericVector x, double cc) {
    NumericVector out(x.size());
    out = x / (3.25 * pow(cc, 2));
    out[abs(x) > 2*cc] = psiAux(x[abs(x) > 2 * cc], cc);
    out[abs(x) > 3*cc] = 0;
    return out;
}

// [[Rcpp::export]]
NumericVector derpsiAux(NumericVector x, double cc) {
    return (-1.944 / pow(cc, 2) +
            5.184 * pow(x, 2) / pow(cc, 4) -
            1.56 * pow(x, 4) / pow(cc, 6) +
            0.112 * pow(x, 6) / pow(cc, 8)) / 3.25;
}

// [[Rcpp::export]]
NumericVector derpsiOpt(NumericVector x, double cc) {
    //' derpsiOpt
    //' the derivative of the psi function  (see \code{\link{psiOpt}})
    //'
    //' @param x A real number.
    //' @param cc a tunning constant.
    //' @return rho'(\code{x} / \code{c}\code{c})
    //' @examples
    //' psiOpt(x=0.5,cc=1)
    //'
    //' @export
    NumericVector out(x.size(), 1.0 / (3.25 * pow(cc, 2)));
    //derivative x^1
    //tmp2 <- (-1.944 * x / cc^2 + 1.728 * x^3 / cc^4 - 0.312 * x^5 / cc^6 + 0.016 * x^7 / cc^8) / 3.25
    // derivative of expresion above !!
    //tmp2 <- (-1.944/ cc^2 + (3*1.728) * x^2 / cc^4 - (5*0.312) * x^4 / cc^6 + (7*0.016) * x^6 / cc^8) / 3.25
    out[abs(x) > 2*cc] = derpsiAux(x[abs(x) > 2*cc], cc);
    out[abs(x) > 3*cc] = 0;
    return out;
}

//'  normal_consistency_constants
//'
//' @description constants previously computed in order to the M scale
//'   value be consistent with the standard normal distribution for
//'   the optimal rho function considered in \code{\link{rhoOpt}}.
//'   (Constant were computed from p=1 till p =400)
//'
//' @param p dimension where observation lives
//' @return cvalue
//'
//' @examples
//' p=5;
//' n=1000
//' X=matrix(rnorm(n*p), ncol=p)
//' dist=apply(X,1,function(t){sqrt(sum(t^2))})
//' s= Mscale(dist,b=0.5, c=normal_consistency_constants(p))
//'
//' ### variable s should be near from one for all p values between 1 and 400.
//'
//' @references [1] Maronna, R. A., Martin, R. D., Yohai, V. J., &
//'   Salibián-Barrera, M. (2018).
//' Robust statistics: theory and methods (with R). Wiley.
//' [2] Salibian-Barrera, M., Willems, G., & Zamar, R. (2008).
//' The fast-tau estimator for regression.
//' Journal of Computational and Graphical Statistics, 17(3), 659-682.
//'
//' @export
// [[Rcpp::export]]
double normal_consistency_constants(int p) {
    NumericVector vaux = {0.404629,0.6944748,0.8985921,1.063144,1.204321,1.329791,1.443817,1.548994,1.647149,1.739537,1.827075,1.910406,1.99017,2.066772,2.140529,2.211772,2.280742,2.347639,2.412622, 2.475882,2.537545,2.597723,2.656494,2.714016,2.770276,2.825434,2.879547,2.932612,2.984741,3.035955,3.08632,3.135869,3.184648, 3.232684,3.279986,3.326633,3.372634,3.418005,3.462781,3.506981, 3.550627,3.593741,3.636342,3.678449,3.720075,3.761236,3.80195, 3.842231,3.88208,3.921557,3.960609,3.999286,4.037581,4.075532, 4.113132,4.15039,4.187306,4.223932,4.260211,4.296178,4.331845, 4.367245,4.402347,4.437173,4.471721,4.506001,4.540023,4.573796, 4.60732,4.6406,4.673648,4.706454,4.739037,4.771396,4.803549, 4.835481,4.867181,4.898696,4.930009,4.961124,4.992046,5.022776, 5.053319,5.083678,5.113838,5.143877,5.173678,5.203352,5.232814, 5.262137,5.291296,5.320304,5.34915,5.377833,5.406382,5.434777, 5.463004,5.491069,5.519042,5.546857,5.574513,5.602029,5.629455, 5.656713,5.683843,5.710848,5.737726,5.764478,5.791107,5.817614, 5.844001,5.870269,5.89642,5.922455,5.948373,5.974178,5.999901, 6.025461,6.05095,6.076328,6.101599,6.126763,6.151823,6.176781, 6.201638,6.226395,6.251055,6.275624,6.300084,6.324457,6.348727,6.372907,6.397029,6.421037,6.444946,6.468764,6.492492,6.516168,6.539737,6.563222,6.586623,6.609943,6.63318,6.656336,6.679371,6.702367,6.725285,6.748124,6.770886,6.793581,6.816188,6.83873,6.861204,6.883556,6.905896,6.928129,6.950294,6.972394,6.994428,7.016405,7.038284,7.060111,7.081862,7.103564,7.12519,7.146749,7.168252,7.189681,7.211047,7.232349,7.253587,7.274746,7.295884,7.316926,7.337911,7.358844,7.379737,7.400556,7.421316,7.442017,7.462661,7.483248,7.503778,7.524252,7.54467,7.565039,7.585339,7.605592,7.625797,7.645945,7.666039,7.68608,7.706069,7.726007,7.74591,7.765748,7.785535,7.805273,7.82496,7.844558,7.864147,7.883689,7.903179,7.922623,7.942019,7.961367,7.980669,7.999924,8.019132,8.038295,8.057414,8.076486,8.095513,8.114495,8.133433,8.152326,8.171176,8.189983,8.208746,8.227467,8.246149,8.264782,8.283357,8.301944,8.320429,8.338918,8.35733,8.375742,8.394071,8.412383,8.430632,8.448894,8.46707,8.485221,8.503334,8.521408,8.539443,8.557449,8.575411,8.593294,8.61118,8.629047,8.646855,8.664627,8.682362,8.700059,8.717721,8.735375,8.752975,8.7705,8.788031,8.805526,8.822969,8.840417,8.857788,8.875119,8.892458,8.909744,8.926998,8.94422,8.961409,8.978539,8.995658,9.012743,9.029797,9.046818,9.063807,9.080764,9.097691,9.114586,9.131449,9.148281,9.165083,9.181853,9.198593,9.215303,9.231982,9.248631,9.265251,9.28184,9.2984,9.314931,9.331432,9.347904,9.364348,9.380762,9.397148,9.413505,9.429834,9.446134,9.462407,9.478651,9.494868,9.511057,9.527219,9.543353,9.55946,9.57554,9.591593,9.607619,9.623618,9.639591,9.655538,9.671458,9.687352,9.70322,9.719062,9.734879,9.75067,9.766435,9.782175,9.797895,9.813593,9.829224,9.84487,9.86049,9.876084,9.891652,9.907196,9.922714,9.938208,9.953677,9.969121,9.984542,9.999938,10.01531,10.03066,10.04598,10.06129,10.07656,10.09182,10.10705,10.12226,10.13745,10.15261,10.16776,10.18288,10.19798,10.21304,10.22809,10.24312,10.25813,10.27312,10.2881,10.30303,10.31794,10.33285,10.34773,10.3626,10.37743,10.39222,10.40702,10.42181,10.43656,10.45129,10.46599,10.48068,10.49537,10.51,10.52463,10.53923,10.55382,10.56839,10.58294,10.59746,10.61197,10.62646,10.64092,10.65537,10.6698,10.68421,10.6986,10.71297,10.72732,10.74165,10.75597,10.77026,10.78454,10.79879,10.81303,10.82725,10.84145,10.85564,10.8698,10.88395,10.89807,10.91218,10.92627,10.9403,10.95436,10.96841,10.98243,10.99643,11.01041,11.02437,11.03832,11.05225,11.06617,11.08006,11.09393};
    return vaux(p-1);
}

// [[Rcpp::export]]
NumericVector constC1(NumericVector p) {
    NumericVector out(p.size(), 1);
    return out;
}

// [[Rcpp::export]]
NumericVector constC2(NumericVector p) {
    double aMaronna = 2.9987;
    double bMaronna = -0.4647;
    return aMaronna * pow(p, bMaronna);
}

// [[Rcpp::export]]
double my_median(NumericVector x) {
    std::size_t size = x.size();
    std::sort(x.begin(), x.end());
    if (size  % 2 == 0) return (x[size / 2 - 1] + x[size / 2]) / 2.0;
    return x[size / 2];
}

//' Mscale
//' the M scale of an univariate sample (see reference below)
//'
//' @param u an univariate sample of size n.
//' @param b the desired break down point
//' @param cc a tuning constant, if consistency to standard normal
//'   distribution is desired use
//' \code{\link{normal_consistency_constants}}
//' @return the Mscale value
//' @examples
//' Mscale(u=rnorm(100),c=1)
//'
//' @importFrom stats median
//' @references Maronna, R. A., Martin, R. D., Yohai, V. J., & Salibian-Barrera, M. (2018).
//' Robust statistics: theory and methods (with R). Wiley.
//' @export
// [[Rcpp::export]]
double Mscale(NumericVector u, double b, double cc) {
    double sn = my_median(abs(u)) / 0.6745;
    if (sn == 0) return sn;

    NumericVector var;
    double quantity = mean(rhoOpt(var, cc)) - b; // indicates if the root is on the right
    if (quantity == 0) return(sn);

    while (quantity > 0) {
        sn = 1.5 * sn;
        quantity = mean(rhoOpt(u/sn, cc)) - b;
    }

    int i = 0;
    double err = 1;
    double AA, BB;
    //  double B;
    double factorAB;
    while ((i < 1000) & (err > 1e-10)) {
        var = u / sn;
        AA = mean(rhoOpt(var, cc));
        BB = mean(psiOpt(var, cc) * var);
        factorAB = (AA - BB - b) / (2.0 * AA - BB - 2.0 * b);
        err = fabs(factorAB - 1);
        sn = sn * fabs(factorAB);
        i = i + 1;
    }
    return sn;
}
