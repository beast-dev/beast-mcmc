package dr.evomodel.continuous;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;

import java.util.List;

/**
 * Integrated multivariate trait likelihood that assumes a fully-conjugate prior on the root.
 * The fully-conjugate prior is a multivariate normal distribution with a precision scaled by
 * diffusion process
 *
 * @author Marc A. Suchard
 */
public class FullyConjugateMultivariateTraitLikelihood extends IntegratedMultivariateTraitLikelihood {

    public FullyConjugateMultivariateTraitLikelihood(String traitName,
                                                     TreeModel treeModel,
                                                     MultivariateDiffusionModel diffusionModel,
                                                     CompoundParameter traitParameter,
                                                     List<Integer> missingIndices,
                                                     boolean cacheBranches,
                                                     boolean scaleByTime,
                                                     boolean useTreeLength,
                                                     BranchRateModel rateModel,
                                                     Model samplingDensity,
                                                     boolean reportAsMultivariate,
                                                     double[] rootPriorMean,
                                                     double rootPriorSampleSize,
                                                     boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);

        // fully-conjugate multivariate normal with own mean and prior sample size
        this.rootPriorMean = rootPriorMean;
        this.rootPriorSampleSize = rootPriorSampleSize;

        priorInformationKnown = false;
    }


//    private double getLogPrecisionDetermination() {
//        return Math.log(diffusionModel.getDeterminantPrecisionMatrix()) + dimTrait * Math.log(rootPriorSampleSize);
//    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == diffusionModel) {
            priorInformationKnown = false;
        }
    }

    public void restoreState() {
        super.restoreState();
        priorInformationKnown = false;
    }

    public void makeDirty() {
        super.makeDirty();
        priorInformationKnown = false;
    }

    protected double integrateLogLikelihoodAtRoot(double[] conditionalRootMean,
                                                  double[] marginalRootMean,
                                                  double[][] notUsed,
                                                  double[][] treePrecisionMatrix, double conditionalRootPrecision) {
        final double square;
        final double marginalPrecision = conditionalRootPrecision + rootPriorSampleSize;
        final double marginalVariance = 1.0 / marginalPrecision;

        // square : (Ay + Bz)' (A+B)^{-1} (Ay + Bz)

        // A = conditionalRootPrecision * treePrecisionMatrix
        // B = rootPriorSampleSize * treePrecisionMatrix

        if (dimTrait > 1) {

            computeWeightedAverage(conditionalRootMean, conditionalRootPrecision, rootPriorMean, rootPriorSampleSize,
                    dimTrait, marginalRootMean);
            
            square = computeQuadraticProduct(marginalRootMean, treePrecisionMatrix, marginalRootMean, dimTrait) * marginalPrecision;

        } else {
            // 1D is very simple
            final double x = conditionalRootMean[0] * conditionalRootPrecision + rootPriorMean[0] * rootPriorSampleSize;
            square = x * x * treePrecisionMatrix[0][0] * marginalVariance;
        }

        if (!priorInformationKnown) {
            setRootPriorSumOfSquares(treePrecisionMatrix);
        }

        final double retValue = 0.5 * (dimTrait * Math.log(rootPriorSampleSize * marginalVariance) - zBz + square);

        if (DEBUG) {
            System.err.println("(Ay+Bz)(A+B)^{-1}(Ay+Bz) = " + square);
            System.err.println("density = " + retValue);
            System.err.println("zBz = " + zBz);
        }

        return retValue;
    }

    protected double integrateLogLikelihoodAtRootFromFullTreeMatrix(double[][] treeTraitPrecisionMatrix,
                                                                    double[] tipTraits) {

        throw new RuntimeException("Not yet implemented (for debugging purposes only)");
    }

    private void setRootPriorSumOfSquares(double[][] treePrecisionMatrix) {

        zBz = computeQuadraticProduct(rootPriorMean, treePrecisionMatrix, rootPriorMean, dimTrait) * rootPriorSampleSize;
        priorInformationKnown = true;
    }

    protected double[][] computeMarginalRootMeanAndVariance(double[] conditionalRootMean, 
                                                            double[][] notUsed,
                                                            double[][] treeVarianceMatrix, 
                                                            double conditionalRootPrecision) {

        final double[][] outVariance = tmpM; // Use a temporary buffer, will stay valid for only a short while

        computeWeightedAverage(conditionalRootMean, conditionalRootPrecision, rootPriorMean, rootPriorSampleSize,
                dimTrait, conditionalRootMean);

        final double totalVariance = 1.0 / (conditionalRootPrecision + rootPriorSampleSize);
        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                outVariance[i][j] = treeVarianceMatrix[i][j] * totalVariance;
            }
        }

        return outVariance;
    }

    private double[] rootPriorMean;
    private double rootPriorSampleSize;

    private boolean priorInformationKnown = false;
    private double zBz; // Prior sum-of-squares contribution
}
