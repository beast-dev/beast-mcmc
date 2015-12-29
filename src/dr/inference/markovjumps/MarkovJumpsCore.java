/*
 * MarkovJumpsCore.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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
        expEvalScalar = new double[stateCount];
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

        for (int i = 0; i < stateCount; i++) {
            expEvalScalar[i] = Math.exp(eval[i] * scalar);
        }
       
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (Math.abs(eval[i] - eval[j]) < 1E-7) {
                    auxInt[index] = expEvalScalar[i] * scalar;
                } else {
                    auxInt[index] = (expEvalScalar[i] -
                            expEvalScalar[j]) /
                            (eval[i] - eval[j]);
                }
                index++;
            }
        }
    }

    public void computeCondStatMarkovJumps(double[] evec,
                                           double[] ievc,
                                           double[] eval,
                                           double[] rateReg,
                                           double   time,
                                           double[] transitionProbs,
                                           double[] countMatrix) {
        computeJointStatMarkovJumps(evec, ievc, eval, rateReg, time, countMatrix);
        for(int i=0; i<stateCount2; i++) {
            countMatrix[i] /= transitionProbs[i];
        }
    }

    public void computeCondStatMarkovJumpsPrecompute(double[] evec,
                                           double[] ievc,
                                           double[] eval,
                                           double[] ievcRateRegEvc,
                                           double   time,
                                           double[] transitionProbs,
                                           double[] countMatrix) {
        computeJointStatMarkovJumpsPrecompute(evec, ievc, eval, ievcRateRegEvc, time, countMatrix);
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

    public void computeJointStatMarkovJumps(double[] evec,
                                            double[] ievc,
                                            double[] eval,
                                            double[] rateReg,
                                            double   time,
                                            double[] countMatrix) {
        // Equation (37) from Minin and Suchard
        populateAuxInt(eval,time,auxInt);

        // Equation (36) from Minin and Suchard
        // Take rate.reg%*%rate.eigen$vectors
        matrixMultiply(rateReg, evec, stateCount, tmp1);

        // Take int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors)
        matrixMultiply(ievc, tmp1, stateCount, tmp2);
        for (int i = 0; i < stateCount2; i++) {
            tmp2[i] *= auxInt[i];
        }

        // Take (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
        //        rate.eigen$invvectors
        matrixMultiply(tmp2, ievc, stateCount, tmp1);

        // Take factorial.moments = rate.eigen$vectors%*%
        //      (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
        //        rate.eigen$invvectors
        matrixMultiply(evec, tmp1, stateCount, countMatrix);
    }

    public void computeJointStatMarkovJumpsPrecompute(double[] evec,
                                            double[] ievc,
                                            double[] eval,
                                            double[] ievcRateRegEvc,
                                            double   time,
                                            double[] countMatrix) {
        // Equation (37) from Minin and Suchard
        populateAuxInt(eval,time,auxInt);

        // Equation (36) from Minin and Suchard                      
        // Take int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors)
        for (int i = 0; i < stateCount2; i++) {
            tmp2[i] = auxInt[i] * ievcRateRegEvc[i];
        }

        // Take (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
        //        rate.eigen$invvectors
        matrixMultiply(tmp2, ievc, stateCount, tmp1);

        // Take factorial.moments = rate.eigen$vectors%*%
        //      (int.matrix*(rate.eigen$invvectors%*%rate.reg%*%rate.eigen$vectors))%*%
        //        rate.eigen$invvectors
        matrixMultiply(evec, tmp1, stateCount, countMatrix);
    }

    // Computes C = A %*% B for square matrices A and B
    public static void matrixMultiply(final double[] A,
                                      final double[] B,
                                      final int dim,
                                      final double[] C) {
        int index = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                C[index] = 0;
                for (int k = 0; k < dim; k++) {
                    C[index] += A[i * dim + k] * B[k * dim + j];
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
       fillRegistrationMatrix(matrix,from,to,dim,1.0);
    }

    public static void fillRegistrationMatrix(double[] matrix, int from, int to, int dim, double value) {
        Arrays.fill(matrix,0.0);
        matrix[from*dim + to] = value;
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
        if (matrix.length == 16) {
            swapRows(matrix,1,2,4);
            swapCols(matrix,1,2,4);
        } else if (matrix.length == 4) {
            swapCols(matrix,1,2,1);
        } else {
             throw new RuntimeException("Function constructed for nucleotides");
        }
    }

    private int stateCount;
    private int stateCount2;
    private double[] auxInt;
    private double[] tmp1;
    private double[] tmp2;
    private double[] expEvalScalar;
}
