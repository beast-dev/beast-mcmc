package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.WrappedNormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.evomodelxml.continuous.hmc.IntegratedLoadingsGradientParser;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.matrixAlgebra.*;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.util.StopWatch;
import dr.util.TaskPool;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.continuous.hmc.LinearOrderTreePrecisionTraitProductProvider.castTreeTrait;
import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;
import static dr.math.matrixAlgebra.missingData.MissingOps.weightedAverage;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrook
 */
public class IntegratedLoadingsGradient implements GradientWrtParameterProvider, VariableListener, Reportable {

    private final TreeTrait<List<WrappedNormalSufficientStatistics>> fullConditionalDensity;
    private final IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood;
    private final ContinuousTraitPartialsProvider partialsProvider;
    protected final int dimTrait;
    protected final int dimFactors;
    protected final int dimPartials;
    private final Tree tree;
    private final Likelihood likelihood;
    private final double[] data;
    private final boolean[] missing;
    private final ThreadUseProvider threadUseProvider;
    private final RemainderCompProvider remainderCompProvider;
    private final TaskPool taskPool;

    public IntegratedLoadingsGradient(TreeDataLikelihood treeDataLikelihood,
                                      ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                      IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                      ContinuousTraitPartialsProvider partialsProvider,
                                      TaskPool taskPool,
                                      ThreadUseProvider threadUseProvider,
                                      RemainderCompProvider remainderCompProvider) {


        this.factorAnalysisLikelihood = factorAnalysisLikelihood;
        this.partialsProvider = partialsProvider;

        String traitName = factorAnalysisLikelihood.getModelName();

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);

        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
        this.tree = treeDataLikelihood.getTree();

        this.dimTrait = factorAnalysisLikelihood.getDataDimension();
        this.dimFactors = factorAnalysisLikelihood.getNumberOfFactors();
        this.dimPartials = partialsProvider.getTraitDimension();

        Parameter dataParameter = factorAnalysisLikelihood.getParameter();
        this.data = dataParameter.getParameterValues();
        dataParameter.addVariableListener(this);

        this.missing = getMissing(factorAnalysisLikelihood.getMissingDataIndices(), dataParameter.getDimension());

        List<Likelihood> likelihoodList = new ArrayList<>();
        likelihoodList.add(treeDataLikelihood);
        likelihoodList.add(factorAnalysisLikelihood);
        this.likelihood = new CompoundLikelihood(likelihoodList);

        this.taskPool = (taskPool != null) ? taskPool :
                new TaskPool(tree.getExternalNodeCount(), 1);

        if (this.taskPool.getNumTaxon() != tree.getExternalNodeCount()) {
            throw new IllegalArgumentException("Incorrectly specified TaskPool");
        }

        this.threadUseProvider = threadUseProvider;
        this.remainderCompProvider = remainderCompProvider;

        if (TIMING) {
            int length = 5;
            stopWatches = new StopWatch[length];
            for (int i = 0; i < length; ++i) {
                stopWatches[i] = new StopWatch();
            }

            System.out.println("WARNING: " + IntegratedLoadingsGradientParser.PARSER_NAME +
                    " is running serially (not in parallel).");

        }
    }

    private boolean[] getMissing(List<Integer> missingIndices, int length) {
        boolean[] missing = new boolean[length];

        for (int i : missingIndices) {
            missing[i] = true;
        }

        return missing;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return factorAnalysisLikelihood.getLoadings();
    }

    @Override
    public int getDimension() {
        return dimFactors * dimTrait;
    }

    private int getGradientDimension() {
        return dimFactors * dimTrait;
    }

    private ReadableMatrix shiftToSecondMoment(WrappedMatrix variance, ReadableVector mean) {

        assert (variance.getMajorDim() == variance.getMinorDim());
        assert (variance.getMajorDim() == mean.getDim());

        final int dim = variance.getMajorDim();

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                variance.set(i, j, variance.get(i, j) + mean.get(i) * mean.get(j));
            }
        }

        return variance;
    }

    private static WrappedNormalSufficientStatistics getWeightedAverage(ReadableVector m1, ReadableMatrix p1,
                                                                        ReadableVector m2, ReadableMatrix p2) {

        assert (m1.getDim() == m2.getDim());
        assert (p1.getDim() == p2.getDim());

        assert (m1.getDim() == p1.getMinorDim());
        assert (m1.getDim() == p1.getMajorDim());

        int dim = m1.getDim();

        final WrappedVector m12 = new WrappedVector.Raw(new double[m1.getDim()], 0, dim);
        final DenseMatrix64F p12 = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F v12 = new DenseMatrix64F(dim, dim);

        final WrappedMatrix wP12 = new WrappedMatrix.WrappedDenseMatrix(p12);
        final WrappedMatrix wV12 = new WrappedMatrix.WrappedDenseMatrix(v12);

        MissingOps.add(p1, p2, wP12);
        safeInvert2(p12, v12, false);

        weightedAverage(m1, p1, m2, p2, m12, wV12, dim);

        return new WrappedNormalSufficientStatistics(m12, wP12, wV12);
    }

    @Override
    public double[] getGradientLogDensity() {

        if (TIMING) {
            stopWatches[2].start();
        }

        final double[][] gradients = new double[this.taskPool.getNumThreads()][getGradientDimension()];

        final ReadableVector gamma = new WrappedVector.Parameter(factorAnalysisLikelihood.getPrecision());
        final ReadableMatrix loadings = ReadableMatrix.Utils.transposeProxy(
                new WrappedMatrix.MatrixParameter(factorAnalysisLikelihood.getLoadings()));

        final double[] rawGamma = factorAnalysisLikelihood.getPrecision().getParameterValues();
        final double[] transposedLoadings = ReadableMatrix.Utils.toArray(
                new WrappedMatrix.MatrixParameter(factorAnalysisLikelihood.getLoadings()));

        if (DEBUG) {
            System.err.println("G : " + gamma);
            System.err.println("L : " + loadings);
        }

        assert (gamma.getDim() == dimTrait);
        assert (loadings.getMajorDim() == dimFactors);
        assert (loadings.getMinorDim() == dimTrait);

        // [E(F) Y^t - E(FF^t)L]\Gamma
        // E(FF^t) = V(F) + E(F)E(F)^t

        // Y: N x P
        // F: N x K
        // L: K x P

        // (K x N)(N x P)(P x P) - (K x N)(N x K)(K x P)(P x P)
        // sum_{N} (K x 1)(1 x P)(P x P) - (K x K)(K x P)(P x P)

        if (TIMING) {
            stopWatches[2].stop();
            stopWatches[3].start();
        }

        if (remainderCompProvider.computeRemainder()) {
            likelihood.getLogLikelihood(); // forces full likelihood evaluation (shouldn't add extra work beyond computing remainders)
        }

        final List<WrappedNormalSufficientStatistics> allStatistics =
                fullConditionalDensity.getTrait(tree, null); // TODO Need to test if faster to load inside loop


        if (TIMING) {
            stopWatches[3].stop();
        }

        assert (allStatistics.size() == tree.getExternalNodeCount());

        if (!threadUseProvider.usePool() || TIMING) {
            for (int taxon = 0, end = tree.getExternalNodeCount(); taxon < end; ++taxon) {
                computeGradientForOneTaxon(0, taxon,
                        loadings, transposedLoadings,
                        gamma, rawGamma,
                        allStatistics.get(taxon), gradients);
            }
        } else {

            this.taskPool.fork(
                    (taxon, thread) -> computeGradientForOneTaxon(thread, taxon,
                            loadings, transposedLoadings,
                            gamma, rawGamma,
                            allStatistics.get(taxon), gradients)
            );
        }

        return join(gradients);
    }


    private void computeGradientForOneTaxon(final int index,
                                            final int taxon,
                                            final ReadableMatrix loadings,
                                            final double[] transposedLoadings,
                                            final ReadableVector gamma,
                                            final double[] rawGamma,
                                            final WrappedNormalSufficientStatistics statistic,
                                            final double[][] gradArray) {

        if (TIMING) {
            stopWatches[0].start();
        }

//        final WrappedVector y = getTipData(taxon);
        final WrappedNormalSufficientStatistics dataKernel = getTipKernel(taxon);

        final ReadableVector meanKernel = dataKernel.getMean();
        final ReadableMatrix precisionKernel = dataKernel.getPrecision();

        if (DEBUG) {
//            System.err.println("Y " + taxon + " : " + y);
            System.err.println("YM" + taxon + " : " + meanKernel);
            System.err.println("YP" + taxon + " : " + precisionKernel);
        }

//        for (WrappedNormalSufficientStatistics statistic : statistics) {  // TODO Maybe need to re-enable

        final ReadableVector meanFactor = statistic.getMean();
        final WrappedMatrix precisionFactor = statistic.getPrecision();
        final WrappedMatrix varianceFactor = statistic.getVariance();

        if (DEBUG) {
            System.err.println("FM" + taxon + " : " + meanFactor);
            System.err.println("FP" + taxon + " : " + precisionFactor);
            System.err.println("FV" + taxon + " : " + varianceFactor);
        }

        WrappedNormalSufficientStatistics convolution = getWeightedAverage(
                meanFactor, precisionFactor,
                meanKernel, precisionKernel);

        convolution = partialsProvider.partitionNormalStatistics(convolution, factorAnalysisLikelihood);

        final ReadableVector mean = convolution.getMean();
//            final ReadableMatrix precision = convolution.getPrecision();
        final WrappedMatrix variance = convolution.getVariance();

        if (DEBUG) {
            System.err.println("CM" + taxon + " : " + mean);
//                System.err.println("CP" + taxon + " : " + precision);
            System.err.println("CV" + taxon + " : " + variance);
        }

        final ReadableMatrix secondMoment = shiftToSecondMoment(variance, mean);
//            final ReadableMatrix product = ReadableMatrix.Utils.productProxy(
//                    secondMoment, loadings
//            );

        if (DEBUG) {
            System.err.println("S" + taxon + " : " + secondMoment);
//                System.err.println("P" + taxon + " : " + product);
        }

        double[] moment = ReadableMatrix.Utils.toArray(secondMoment);

        if (TIMING) {
            stopWatches[0].stop();
            stopWatches[1].start();
        }

        for (int factor = 0; factor < dimFactors; ++factor) {
            double factorMean = mean.get(factor);

            for (int trait = 0; trait < dimTrait; ++trait) {
                if (!missing[taxon * dimTrait + trait]) {

                    double product = 0.0;
                    for (int k = 0; k < dimFactors; ++k) {
                        product += moment[factor * dimFactors + k] // secondMoment.get(factor, k)
                                * transposedLoadings[trait * dimFactors + k]; // loadings.get(k, trait);
                    }

                    gradArray[index][factor * dimTrait + trait] +=
                            (factorMean // mean.get(factor)
                                    * data[taxon * dimTrait + trait] //y.get(trait)
                                    - product)
//                                         - product.get(factor, trait))
                                    * rawGamma[trait]; // gamma.get(trait);

                }
            }
        }

        if (TIMING) {
            stopWatches[1].stop();
        }
//        }
    }

    private static double[] join(double[][] array) {

        int nRows = array.length;
        int nCols = array[0].length;
        double[] result = array[0];

        for (int row = 1; row < nRows; ++row) {
            double[] temp = array[row];
            for (int col = 0; col < nCols; ++col) {
                result[col] += temp[col];
            }
        }

        return result;
    }

//    private WrappedVector getTipData(int taxonIndex) {
//        return new WrappedVector.Parameter(dataParameter, taxonIndex * dimTrait, dimTrait);
//    }

    private WrappedNormalSufficientStatistics getTipKernel(int taxonIndex) {
        double[] buffer = partialsProvider.getTipPartial(taxonIndex, false);
        return new WrappedNormalSufficientStatistics(buffer, 0, dimPartials, null, PrecisionType.FULL);
    }

    public enum ThreadUseProvider {
        PARALLEL {
            @Override
            boolean usePool() {
                return true;
            }
        },

        SERIAL {
            @Override
            boolean usePool() {
                return false;
            }
        };

        abstract boolean usePool();


    }

    @Override
    public String getReport() {

        String report = "";

        if (TIMING) {
            report += timingInfo();
        }

        report += GradientWrtParameterProvider.getReportAndCheckForError(
                this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                null);

        if (TIMING) {
            report += timingInfo();
        }

        return report;
    }

    private String timingInfo() {
        StringBuilder sb = new StringBuilder("\nTiming in IntegratedLoadingsGradient\n");
        for (StopWatch stopWatch : stopWatches) {
            sb.append("\t").append(stopWatch.toString()).append("\n");
            stopWatch.reset();
        }
        return sb.toString();
    }

    public enum RemainderCompProvider {

        FULL {
            @Override
            boolean computeRemainder() {
                return true;
            }
        },

        SKIP {
            @Override
            boolean computeRemainder() {
                return false;
            }
        };

        abstract boolean computeRemainder();


    }

    private StopWatch[] stopWatches;
    private static final boolean TIMING = false;

    private static final boolean DEBUG = false;


    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        throw new RuntimeException("Trait data is not cached");
    }
}
