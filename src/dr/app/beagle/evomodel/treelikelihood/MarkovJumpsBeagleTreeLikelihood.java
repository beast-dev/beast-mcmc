package dr.app.beagle.evomodel.treelikelihood;

//import dr.evolution.alignment.PatternList;
//import dr.evolution.datatype.DataType;
//import dr.evomodel.tree.TreeModel;
//import dr.evomodel.branchratemodel.BranchRateModel;
//import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
//import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
//import dr.inference.model.MatrixParameter;

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
public class MarkovJumpsBeagleTreeLikelihood  { //extends AncestralStateBeagleTreeLikelihood {

    public MarkovJumpsBeagleTreeLikelihood(int stateCount) {
        this.stateCount = stateCount;
    }

//    public MarkovJumpsBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
//                                           BranchSiteModel branchSiteModel, SiteRateModel siteRateModel,
//                                           BranchRateModel branchRateModel, boolean useAmbiguities,
//                                           PartialsRescalingScheme scalingScheme, DataType dataType, String tag,
//                                           SubstitutionModel substModel, MatrixParameter registerMatrixParameter) {
//        super(patternList, treeModel, branchSiteModel, siteRateModel, branchRateModel, useAmbiguities,
//                scalingScheme, dataType, tag, substModel);
//
//        this.registerMatrixParameter = registerMatrixParameter;
//        addVariable(registerMatrixParameter);
//        fillRegisterMatrix(registerMatrixParameter, registerMatrix);
//    }

//    private void fillRegisterMatrix(MatrixParameter matrixParameter, double[] matrix) {
//        final int colDim = matrixParameter.getColumnDimension();
//        final int rowDim = matrixParameter.getRowDimension();
//        if (matrix == null ||
//                matrix.length != colDim * rowDim) {
//            matrix = new double[colDim * rowDim];
//        }
//        int index = 0;
//        for(int i=0; i<rowDim; i++) {
//            for(int j=0; j<colDim; j++) {
//                if (i == j) {
//                    matrix[index] = 0;
//                } else {
//                    matrix[index] = matrixParameter.getParameterValue(i,j);
//                }
//                index++;
//            }
//        }
//    }

    protected void fillRegisterMatrix(double[] matrix) {
        if (matrix == null) {
            matrix = new double[stateCount*stateCount];
        }
        Arrays.fill(matrix,1.0);
        for(int i=0; i<stateCount; i++) {
            matrix[i*stateCount + i] = 0;
        }
    }

    protected void computeMeanMarkovJumps(SubstitutionModel substModel,
                                        double[] Lr,
                                        double time,
                                        double[] countMatrix) {

        if (transitionProbs == null) {
            transitionProbs = new double[stateCount * stateCount];
        }

        if (spectralMatrices == null) {
            spectralMatrices = new double[stateCount][];
            for(int i=0; i<stateCount; i++) {
                spectralMatrices[i] = new double[stateCount * stateCount];
            }
        }

        // Get transition probability matrix
        substModel.getTransitionProbabilities(time,transitionProbs);
        Arrays.fill(countMatrix,0); // Zero all counts

        final EigenDecomposition eigenDecomposition = substModel.getEigenDecomposition();
        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        // Get spectral matrices
        // Equation (34) from Minin and Suchard
        // TODO Only recompute when eigenDecomposition changes
        int index;
        for(int i=0; i<stateCount; i++) {
            // For each i in Equation (34), compute element (jk) of B_i
            double[] Bi = spectralMatrices[i];
            index = 0;
            for(int j=0; j<stateCount; j++) {
                for(int k=0; k<stateCount; k++) {
                    Bi[index] = evec[j * stateCount + i] * ievc[i * stateCount + k];
                    index++;
                }
            }
        }

//        double[] Lr = registerMatrix;

        for(int i=0; i<stateCount; i++) {
            double[] Bi = spectralMatrices[i];
            for(int j=0; j<stateCount; j++) {
                double[] Bj = spectralMatrices[j];

                // Equation (37) from Minin and Suchard
                // Compute I_{ij}
                double auxInc;
                if(Math.abs(eval[i] - eval[j]) < 1E-7) {
                    auxInc = Math.exp(eval[j] * time) * time;
                } else {
                    auxInc = (Math.exp(eval[i] * time) -
                              Math.exp(eval[j] * time)) /
                             (eval[i] - eval[j]);
                }

                // Equation (36) from Minin and Suchard
                index = 0;
                // For each element (kl) in M matrix
                for(int k=0; k<stateCount; k++) {
                    for(int l=0; l<stateCount; l++) {

                        // Do matrix * matrix * matrix
                        for(int m=0; m<stateCount; m++) {
                            for(int n=0; n<stateCount; n++) {
                                countMatrix[index] += Bi[k * stateCount + m] *
                                                      Lr[m * stateCount + n] *
                                                      Bj[n * stateCount + l];
                            }
                        }
                        countMatrix[index] *= auxInc;

                        // Next element
                        index++;
                    }
                }             
            }
        }
    }

    public static void main(String[] args) {
        MarkovJumpsBeagleTreeLikelihood markovjumps = new MarkovJumpsBeagleTreeLikelihood(4);
        markovjumps.computeMeanMarkovJumps(null,null,0,null);
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

//mean.markov.jumps = function(rate.matrix, register.matrix, my.time){
//
//  space.size = nrow(rate.matrix)
//
//  # obtain an eigen decompostition of the rate matrix
//  eigen.struct = eigen(rate.matrix)
//
//  # invert eigen vectors
//
//  inv.eigen = solve(eigen.struct$vectors)
//
//  # get transition probability matrix
//  trans.probs = MatrixExp(rate.matrix, my.time)
//
//  zero.matrix = diag(rep(0,space.size))
//  factorial.moments = zero.matrix
//
//  # get spectral matrices
//  spectral.mat = list(space.size)
//
//  for (i in c(1:space.size)){
//    D.matrix = zero.matrix
//    D.matrix[i,i] = 1
//
//    spectral.mat[[i]] = eigen.struct$vectors%*%D.matrix%*%inv.eigen
//  }
//
//  aux.int = 0
//
//  for (i in c(1:space.size)){
//    for (j in c(1:space.size)){
//      if (abs(eigen.struct$values[i] - eigen.struct$values[j]) < 10^-7){
//        aux.int = exp(eigen.struct$values[j] * my.time)*my.time
//      }else{
//        aux.int = (exp(eigen.struct$values[i] * my.time) -
//                   exp(eigen.struct$values[j] * my.time))/
//                     (eigen.struct$values[i] - eigen.struct$values[j])
//      }
//
//      factorial.moments = factorial.moments +
//        spectral.mat[[i]]%*%register.matrix%*%spectral.mat[[j]]*aux.int
//    }
//  }
//
//  return(factorial.moments/trans.probs)
//}
        
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

    private int stateCount;
//    private double[] registerMatrix;
    private double[] transitionProbs;
    private double[][] spectralMatrices;
//    private MatrixParameter registerMatrixParameter;
}
