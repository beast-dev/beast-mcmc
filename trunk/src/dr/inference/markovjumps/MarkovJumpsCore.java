package dr.inference.markovjumps;

import java.util.Arrays;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *
 * A base class for implementing Markov chain-induced counting processes (markovjumps) in BEAST using BEAGLE
 * This work is supported by NSF grant 0856099
 *
 * Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 * Journal of Mathematical Biology, 56, 391-412.
 *
 */

public class MarkovJumpsCore {

    public MarkovJumpsCore(int stateCount) {
        this.stateCount = stateCount;
        this.stateCount2 = stateCount * stateCount;
        auxInt = new double[stateCount2];
        tmp1 = new double[stateCount2];
        tmp2 = new double[stateCount2];
    }

//SEXP aux_mat1(SEXP x, SEXP y){
//  int i, j, nx;
//  double *vec, scalar;
//  SEXP ans;
//
//  nx = length(x);
//  x = coerceVector(x, REALSXP);
//  y = coerceVector(y, REALSXP);
//
//  vec = REAL(x);
//  scalar = REAL(y)[0];
//
//
//  PROTECT(ans = allocMatrix(REALSXP, nx, nx));
//
//  for(i = 0; i < nx; i++) {
//    for(j = 0; j < nx; j++){
//      if ((vec[i]-vec[j]<10e-7) && (vec[i]-vec[j]>-10e-7)){
//	REAL(ans)[i + nx*j] = exp(vec[j]*scalar)*scalar;
//      }else{
//	REAL(ans)[i + nx*j] = (exp(vec[i]*scalar) - exp(vec[j]*scalar))/(vec[i] - vec[j]);
//      }
//    }
//  }
//
//  UNPROTECT(1);
//  return(ans);
//}

    private void populateAuxInt(double[] eval, double scalar, double[] auxInt) {

        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (Math.abs(eval[i] - eval[j]) < 1E-7) {
                    auxInt[index] = Math.exp(eval[j] * scalar) * scalar;
                } else {
                    auxInt[index] = (Math.exp(eval[i] * scalar) -
                            Math.exp(eval[j] * scalar)) /
                            (eval[i] - eval[j]);
                }
                index++;
            }
        }
    }

    public void computeCondMeanMarkovJumps(double[] evec,
                                           double[] ievc,
                                           double[] eval,
                                           double[] rateReg,
                                           double   time,
                                           double[] transitionProbs,
                                           double[] countMatrix) {
        computeJointMeanMarkovJumps(evec, ievc, eval, rateReg, time, countMatrix);
        for(int i=0; i<stateCount2; i++) {
            countMatrix[i] /= transitionProbs[i];
        }
    }

//joint.mean.markov.jumps = function(rate.eigen, regist.matrix, interval.len){
//
//  if (!("eigen" %in% class(rate.eigen)))
//    stop("Error: object \"rate.eigen\" is not of class \"eigen\"")
//
//  if (!("revmc" %in% class(rate.eigen)))
//    stop("Error: object \"rate.eigen\" is not of class \"revmc\"")
//
//  if (!is.matrix(regist.matrix) || (nrow(regist.matrix) != ncol(regist.matrix)))
//    stop("\"regist.matrix\" must be a square matrix")
//
//  if (nrow(rate.eigen$rate.matrix) != nrow(regist.matrix))
//    stop("dimensions of rate and register matrices do not match")
//
//  if (prod((regist.matrix == 1) + (regist.matrix == 0)) == 0)
//    stop("all entries of \"regist.matrix\" must be either 0 or 1")
//
//  space.size = dim(regist.matrix)[1]
//
//  zero.matrix = matrix(0,space.size,space.size)
//  factorial.moments = zero.matrix
//  rate.reg = rate.eigen$rate.matrix*regist.matrix
//
//  if (-sum(rate.eigen$stationary*diag(rate.eigen$rate.matrix))*interval.len < 0.001){## if time interval is very small do first order Taylor expansion
//    factorial.moments = rate.reg*interval.len +
//      (rate.eigen$rate.matrix%*%rate.reg + rate.reg%*%rate.eigen$rate.matrix)*(interval.len)^2/2 +
//        rate.eigen$rate.matrix%*%rate.reg%*%rate.eigen$rate.matrix*(interval.len)^3/6
//  }else{
//    int.matrix = .Call("aux_mat1", rate.eigen$values, interval.len)
//
//    factorial.moments = rate.eigen$vectors%*%
//      (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
//        rate.eigen$invvectors
//  }
//
//  return(factorial.moments)
//}

    public void computeJointMeanMarkovJumps(double[] evec,
                                            double[] ievc,
                                            double[] eval,
                                            double[] rateReg,
                                            double   time,
                                            double[] countMatrix) {

        // TODO Implement Taylor expansion for small time

        int index;

        // Equation (37) from Minin and Suchard
        populateAuxInt(eval,time,auxInt);

        // Equation (36) from Minin and Suchard
        // Take rate.reg%*%rate.eigen$vectors
        index = 0;
        for(int i=0; i<stateCount; i++) {
            for(int j=0; j<stateCount; j++) {
                tmp1[index] = 0;
                for(int k=0; k<stateCount; k++) {
                    tmp1[index] += rateReg[i * stateCount + k] * evec[k * stateCount + j];
                }
                index++;
            }
        }

        // Take int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors)
        index = 0;
        for(int i=0; i<stateCount; i++) {
            for(int j=0; j<stateCount; j++) {
                tmp2[index] = 0;
                for(int k=0; k<stateCount; k++) {
                    tmp2[index] += ievc[i * stateCount  + k] * tmp1[k * stateCount + j];
                }
                tmp2[index] *= auxInt[i * stateCount + j];
                index++;
            }
        }

        // Take (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
        //        rate.eigen$invvectors
        index = 0;
        for(int i=0; i<stateCount; i++) {
            for(int j=0; j<stateCount; j++) {
                tmp1[index] = 0;
                for(int k=0; k<stateCount; k++) {
                    tmp1[index] += tmp2[i * stateCount + k] * ievc[k * stateCount + j];
                }
                index++;
            }
        }

        //Take factorial.moments = rate.eigen$vectors%*%
        //      (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
        //        rate.eigen$invvectors
        index = 0;
        for(int i=0; i<stateCount; i++) {
            for(int j=0; j<stateCount; j++) {
                countMatrix[index] = 0;
                for(int k=0; k<stateCount; k++) {
                    countMatrix[index] += evec[i * stateCount + k] * tmp1[k * stateCount + j];
                }
                index++;
            }
        }
    }

    public static void fillRegistrationMatrix(double[] matrix, int dim) {
        Arrays.fill(matrix,1.0);
        for(int i=0; i<dim; i++) {
            matrix[i*dim + i] = 0;
        }
    }

    public static void fillRegistrationMatrix(double[] matrix, int from, int to, int dim) {
        Arrays.fill(matrix,0.0);
        matrix[from*dim + to] = 1.0;
    }

    public static void swapRows(double[] matrix, int swap1, int swap2, int dim) {
        for(int i=0; i<dim; i++) {
            double tmp = matrix[swap1 * dim + i];
            matrix[swap1 * dim + i] = matrix[swap2 * dim + i];
            matrix[swap2 * dim + i] = tmp;
        }
    }

    public static void swapCols(double[] matrix, int swap1, int swap2, int dim) {
        for(int i=0; i<dim; i++) {
            double tmp = matrix[i * dim + swap1];
            matrix[i * dim + swap1] = matrix[i * dim + swap2];
            matrix[i * dim + swap2] = tmp;
        }
    }

    public static void makeComparableToRPackage(double[] matrix) {
        if (matrix.length != 16) {
            throw new RuntimeException("Function constucted for nucleotides");
        }
        swapRows(matrix,1,2,4);
        swapCols(matrix,1,2,4);
    }

    private int stateCount;
    private int stateCount2;
    private double[] auxInt;
    private double[] tmp1;
    private double[] tmp2;
}
