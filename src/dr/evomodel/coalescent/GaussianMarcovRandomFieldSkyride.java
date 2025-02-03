package dr.evomodel.coalescent;

import dr.inference.distribution.RandomField;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.GaussianMarkovRandomField;
import no.uib.cipr.matrix.DenseVector;

import java.util.ArrayList;
import java.util.List;

// this code adapts GMRFMultilocusSkyrideLikelihood to the new framework
// TODO still to manage the gradients
public class GaussianMarcovRandomFieldSkyride extends GaussianMarkovRandomField {

    private List<MatrixParameter> covariates;
    private List<Parameter> betaList;
    private List<Parameter> deltaList;

    public static final double LOG_TWO_TIMES_PI = 1.837877;

    private SymmetricTriDiagonalMatrix Q;
    private int[] firstObservedIndex;
    private int[] lastObservedIndex;
    private int[] recIndices;
    private int[] distIndices;
    private final List<Parameter> covPrecParametersRecent;
    private final List<Parameter> covPrecParametersDistant;
    private boolean missingValues;
    private boolean modelWithCovariates;
//    private List<RadomField.WeightProvider> weightProviderList;

    public GaussianMarcovRandomFieldSkyride(int dim,
                                            Parameter precisionParameter,
                                            Parameter logPopSizes,
                                            Parameter lambdaParameter,
                                            RandomField.WeightProvider weightProvider,
                                            List<Parameter> firstObservedIndexParameter,
                                            List<Parameter> lastObservedIndexParameter,
                                            List<Parameter> covPrecParametersRecent,
                                            List<Parameter> covPrecParametersDistant,
                                            Parameter recentIndices,
                                            Parameter distantIndices) {
        super("a", dim, precisionParameter, logPopSizes, lambdaParameter, weightProvider, false);
        // the following is needed for the case with missing values
        // TODO all these indexes should be handled by the weight provider and not here
        if (firstObservedIndexParameter != null) initializeFirstObservedIndex(firstObservedIndexParameter, recentIndices);
        if (lastObservedIndexParameter != null) initializeLastObservedIndex(lastObservedIndexParameter, distantIndices);
        this.covPrecParametersRecent = covPrecParametersRecent;
        this.covPrecParametersDistant = covPrecParametersDistant;
        setupCovs();
        if (lastObservedIndexParameter != null || firstObservedIndexParameter != null) {
            this.missingValues = true;
            throw new RuntimeException("Missing values case not yet implemented");
//            setupGMRFWeightsForMissingCov();
            // TODO transfer to weight provider; it builds on:
            //  firstObservedIndex, covPrecParametersRecent, lastObservedIndex, covPrecParametersDistant
        } // TODO: missing Values ONLY for the case for covariates
    }


    // TODO: this case should be dealt with in the parser
    //   MatrixParameter dMatrix should just be written as a list with one element
    // private MatrixParameter dMatrix;
//    private Parameter betaParameter;
//    public GaussianMarcovRandomFieldwithCovariates(int dim,
//                                                   Parameter precisionParameter,
//                                                   Parameter logPopSizes,
//                                                   Parameter lambdaParameter,
//                                                   RandomField.WeightProvider weightProvider,
//
//                                                   MatrixParameter dMatrix,
//                                                   Parameter betaParameter,
//
//                                                   List<Parameter> firstObservedIndexParameter,
//                                                   List<Parameter> lastObservedIndexParameter,
//                                                   List<Parameter> covPrecParametersRecent,
//                                                   List<Parameter> covPrecParametersDistant,
//                                                   Parameter recentIndices,
//                                                   Parameter distantIndices) {
//        this(dim, precisionParameter, logPopSizes, lambdaParameter, weightProvider,
//                firstObservedIndexParameter, lastObservedIndexParameter, covPrecParametersRecent, covPrecParametersDistant,
//                recentIndices, distantIndices);
//        this.modelWithCovariates = true;
//        this.betaParameter = betaParameter; // emptyness check should be in the parser
//        addVariable(betaParameter);
//        addVariable(dMatrix);
//        this.dMatrix = dMatrix; // emptyness check should be in the parser
//        if (dMatrix.getRowDimension() != dim) {
//            throw new RuntimeException("Covariate dimension (" + dMatrix.getRowDimension() + ") does not match the dimension of the field (" + dim + ")");
//        } else if (dMatrix.getColumnDimension() != dMatrix.getRowDimension()) {
//            throw new RuntimeException("Incorrect covariate dimensions (" + dMatrix.getColumnDimension() + " != "
//                    + dMatrix.getRowDimension() + ")");
//        }
//    }

    public GaussianMarcovRandomFieldSkyride(int dim,
                                            Parameter precisionParameter,
                                            Parameter logPopSizes,
                                            Parameter lambdaParameter,
                                            RandomField.WeightProvider weightProvider,

                                            List<MatrixParameter> covariates,
                                            List<Parameter> betaList,
                                            List<Parameter> deltaList,

                                            List<Parameter> firstObservedIndexParameter,
                                            List<Parameter> lastObservedIndexParameter,
                                            List<Parameter> covPrecParametersRecent,
                                            List<Parameter> covPrecParametersDistant,
                                            Parameter recentIndices,
                                            Parameter distantIndices
                                                   ) {
        this(dim, precisionParameter, logPopSizes, lambdaParameter, weightProvider,
                firstObservedIndexParameter, lastObservedIndexParameter, covPrecParametersRecent, covPrecParametersDistant,
                recentIndices, distantIndices);
        this.modelWithCovariates = true;
        this.covariates = covariates; // TODO emptyness check should be in the parser
        for (MatrixParameter cov : covariates) {
            if (cov.getRowDimension() != dim) {
                throw new RuntimeException("Covariate dimension (" + cov.getRowDimension() + ") does not match the dimension of the field (" + dim + ")");
            }
            addVariable(cov);
        }
        this.betaList = betaList; // TODO emptyness check should be in the parser
        if (betaList.size() != covariates.size()) {
            throw new RuntimeException("beta.size(" + betaList.size() + ") != covariates.size(" + covariates.size() + ")");
        }
        for (Parameter betaParam : betaList) {
            addVariable(betaParam);
        }
        this.deltaList = initializeDelta(deltaList, betaList);
        if (deltaList != null) {
            for (Parameter dParam : deltaList) {
                addVariable(dParam);
            }
        }
    }

    @Override
    public double logPdf(double[] x) {

        double currentLikelihood = 0.0;
        if (missingValues) {
            throw new RuntimeException("Missing values case not yet implemented");
//            currentLikelihood += handleMissingValues(); // TODO this function should act on the weight provider, not here
        }

        DenseVector currentGamma = new DenseVector(x); // TODO: maybe write everything in terms of dense vector? or change just this
        if (modelWithCovariates) {
            updateGammaWithCovariates(currentGamma); // this is a mean adjustment; TODO alternatively update getMean()
        }
        currentLikelihood = logPdf(currentGamma.getData(), getMean(), precisionParameter.getParameterValue(0), getQ(), isImproper(), getLogDeterminant());
        //  TODO manage getMean()

        if (lambdaParameter.getParameterValue(0) == 1) {
            currentLikelihood -= (dim - 1) / 2.0 * LOG_TWO_TIMES_PI;
        } else {
            currentLikelihood -= dim / 2.0 * LOG_TWO_TIMES_PI;
        }
        return currentLikelihood;
    }

    // TODO: these two must go in the weight provider for the missing cov case
    // this class takes the weight matrix for the missing covariates and multiplies the lower triangular
    // part (diag included) of the matrix by the precision parameter
    // it does the same for the firstObs - 2, firstObs - 2 element ???? // TODO ask marc
//    private SymmTridiagMatrix getScaledWeightMatrixForMissingCovRecent(double precision, int covIndex, int firstObs) {
//        SymmTridiagMatrix a = weightMatricesForMissingCovRecent.get(covIndex).copy();
//        for (int i = 0; i < a.numRows() - 1; i++) {
//            a.set(i, i, a.get(i, i) * precision);
//            a.set(i + 1, i, a.get(i + 1, i) * precision);
//        }
//        a.set(firstObs - 2, firstObs - 2,
//                a.get(firstObs - 2, firstObs - 2) * precision);
//        return a;
//    }

// TODO ask marc same as above but for the fieldLength - lastObs - 1, fieldLength - lastObs - 1

//    private SymmTridiagMatrix getScaledWeightMatrixForMissingCovDistant(double precision, int covIndex, int lastObs) {
//        SymmTridiagMatrix a = weightMatricesForMissingCovDistant.get(covIndex).copy();
//        for (int i = 0; i < a.numRows() - 1; i++) {
//            a.set(i, i, a.get(i, i) * precision);
//            a.set(i + 1, i, a.get(i + 1, i) * precision);
//        }
//        a.set(fieldLength - lastObs - 1, fieldLength - lastObs - 1,
//                a.get(fieldLength - lastObs - 1, fieldLength - lastObs - 1) * precision);
//        return a;
//    }

    protected void updateGammaWithCovariates(DenseVector currentGamma) {

        final int N = currentGamma.size();
        double[] update = new double[N];
//        if (dMatrix != null) { // TODO dMatrix can be seen as a list with only one element; this should be handled in the parser
//            final int K = dMatrix.getColumnDimension();
//            for (int i = 0; i < N; ++i) {
//                for (int j = 0; j < K; ++j) {
//                    update[i] += dMatrix.getParameterValue(i, j) * betaParameter.getParameterValue(j);
//                }
//            }
//        }
        for (int k = 0; k < betaList.size(); ++k) {
            Parameter b = betaList.get(k);
            Parameter d = deltaList.get(k);
            final int J = b.getDimension();
            MatrixParameter covariate = covariates.get(k);
            boolean transposed = isTransposed(N, J, covariate);

            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < J; ++j) {
                    update[i] += getCovariateValue(covariate, i, j, transposed) * b.getParameterValue(j) *
                            d.getParameterValue(j);
                }
            }
        }

        for (int i = 0; i < N; ++i) {
            currentGamma.set(i, currentGamma.get(i) - update[i]);
        }
    }

    // TODO: implement
//    private double handleMissingValues() {
//
//        assert (covPrecParametersRecent != null);
//        assert (covariates != null);
//        assert (covPrecParametersDistant != null);
//
//        int numMissing;
//        DenseVector tempVectMissingCov;
//        SymmTridiagMatrix missingCovQ;
//        DenseVector tempVectMissingCov2;
//        int numMissingRecent;
//
//        double currentLike = 0.0;
//
//        if (lastObservedIndex != null) {
//            for (int i = 0; i < covPrecParametersDistant.size(); i++) {
//
//                numMissing = dim - lastObservedIndex[i];
//                tempVectMissingCov = new DenseVector(numMissing);
//                tempVectMissingCov2 = new DenseVector(numMissing);
//
//                missingCovQ = getScaledWeightMatrixForMissingCovDistant(covPrecParametersDistant.get(i).getParameterValue(0), i,
//                        lastObservedIndex[i]);
//
//                for (int j = 0; j < numMissing; j++) {
//                    tempVectMissingCov.set(j, covariates.get(distIndices[i] - 1).getParameterValue(0, lastObservedIndex[i] + j) -
//                            covariates.get(distIndices[i] - 1).getParameterValue(0, lastObservedIndex[i] - 1));
//                }
//
//                missingCovQ.mult(tempVectMissingCov, tempVectMissingCov2);
//                currentLike += 0.5 * (numMissing) * Math.log(covPrecParametersDistant.get(i).getParameterValue(0))
//                        - 0.5 * tempVectMissingCov.dot(tempVectMissingCov2);
//            }
//        }
//
//        if (firstObservedIndex != null) {
//            for (int i = 0; i < covPrecParametersRecent.size(); i++) {
//
//                numMissingRecent = firstObservedIndex[i] - 1;
//                tempVectMissingCov = new DenseVector(numMissingRecent);
//                tempVectMissingCov2 = new DenseVector(numMissingRecent);
//
//                missingCovQ = getScaledWeightMatrixForMissingCovRecent(covPrecParametersRecent.get(i).getParameterValue(0), i,
//                        firstObservedIndex[i]);
//
//                for (int j = 0; j < numMissingRecent; j++) {
//                    tempVectMissingCov.set(j, covariates.get(recIndices[i] - 1).getParameterValue(0, j) -
//                            covariates.get(recIndices[i] - 1).getParameterValue(0, firstObservedIndex[i] - 1));
//                }
//
//                missingCovQ.mult(tempVectMissingCov, tempVectMissingCov2);
//                currentLike += 0.5 * (numMissingRecent) * Math.log(covPrecParametersRecent.get(i).getParameterValue(0))
//                        - 0.5 * tempVectMissingCov.dot(tempVectMissingCov2);
//            }
//        }
//        return currentLike;
//    }


    private void initializeFirstObservedIndex(List<Parameter> firstObservedIndexParameter, Parameter recentIndices) {
        int size = firstObservedIndexParameter.size();
        this.firstObservedIndex = new int[size];
        for (int i = 0; i < size; i++) {
            this.firstObservedIndex[i] = (int) firstObservedIndexParameter.get(i).getParameterValue(0);
        }
        this.recIndices = new int[size];
        if (recentIndices != null) {
            // indices specify which covariates require default unobserved covariate data prior
            for (int i = 0; i < size; i++) {
                this.recIndices[i] = (int) recentIndices.getParameterValue(i);
            }
        } else {
            // If specific covariates not specified by indices, need default unobserved covariate data prior for all covariates
            for (int i = 0; i < size; i++) {
                this.recIndices[i] = i + 1;
            }
        }
    }

    private void initializeLastObservedIndex(List<Parameter> lastObservedIndexParameter, Parameter distantIndices) {
        int size = lastObservedIndexParameter.size();
        this.lastObservedIndex = new int[size];
        for (int i = 0; i < size; i++) {
            this.lastObservedIndex[i] = (int) lastObservedIndexParameter.get(i).getParameterValue(0);
        }
        this.distIndices = new int[size];
        if (distantIndices != null) {
            // indices specify which covariates require default unobserved covariate data prior
            for (int i = 0; i < size; i++) {
                this.distIndices[i] = (int) distantIndices.getParameterValue(i);
            }
        } else {
            // If specific covariates not specified by indices, need default unobserved covariate data prior for all covariates
            for (int i = 0; i < size; i++) {
                this.distIndices[i] = i + 1;
            }
        }
    }

    private List<Parameter> initializeDelta(List<Parameter> deltaList, List<Parameter> betaList) {
        if (deltaList != null) {
            return deltaList;
        }
        List<Parameter> newDelta = new ArrayList<>();
        for (int i = 0; i < betaList.size(); i++) {
            Parameter deltaParam = new Parameter.Default(1.0);
//            deltaParam.setParameterValue(0, 1.0);
            newDelta.add(deltaParam);
        }
        return newDelta;
    }

    private void setupCovs() {
        if (covPrecParametersRecent != null) {
            for (Parameter covPrecRecent : covPrecParametersRecent) {
                addVariable(covPrecRecent);
            }
        }
        if (covPrecParametersDistant != null) {
            for (Parameter covPrecDistant : covPrecParametersDistant) {
                addVariable(covPrecDistant);
            }
        }
    }


    private double getCovariateValue(MatrixParameter matrix, int i, int j, boolean transposed) {
        if (transposed) {
            return matrix.getParameterValue(j, i);
        } else {
            return matrix.getParameterValue(i, j);
        }
    }

    private boolean isTransposed(int N, int J, MatrixParameter matrix) {
        if (J == matrix.getRowDimension() && N == matrix.getColumnDimension()) {
            return true;
        } else if (J == matrix.getColumnDimension() && N == matrix.getRowDimension()) {
            return false;
        } else {
            throw new RuntimeException("Incorrect dimensions in " + matrix.getId() + " (r=" + matrix.getRowDimension() +
                    ",c=" + matrix.getColumnDimension() + ")");
        }
    }


}