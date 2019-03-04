package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.WrappedNormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.*;
import dr.math.matrixAlgebra.missingData.MissingOps;
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
public class IntegratedLoadingsGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeTrait<List<WrappedNormalSufficientStatistics>> fullConditionalDensity;
    private final TreeDataLikelihood treeDataLikelihood;
    private final IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood;
    private final int dimTrait;
    private final int dimFactors;
    private final Tree tree;
    private final Likelihood likelihood;
    private final Parameter data;
    private final boolean[] missing;
    private final TaxonTaskPool taxonTaskPool;

    private IntegratedLoadingsGradient(TreeDataLikelihood treeDataLikelihood,
                                       ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                       IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                       TaxonTaskPool taxonTaskPool) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.factorAnalysisLikelihood = factorAnalysisLikelihood;

        String traitName = factorAnalysisLikelihood.getModelName();

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
        this.tree = treeDataLikelihood.getTree();

        this.dimTrait = factorAnalysisLikelihood.getDataDimension();
        this.dimFactors = factorAnalysisLikelihood.getNumberOfFactors();

        this.data = factorAnalysisLikelihood.getParameter();
        this.missing = getMissing(factorAnalysisLikelihood.getMissingDataIndices(), data.getDimension());

        List<Likelihood> likelihoodList = new ArrayList<Likelihood>();
        likelihoodList.add(treeDataLikelihood);
        likelihoodList.add(factorAnalysisLikelihood);
        this.likelihood = new CompoundLikelihood(likelihoodList);

        this.taxonTaskPool = (taxonTaskPool != null) ? taxonTaskPool :
                new TaxonTaskPool(tree.getExternalNodeCount(), 1);

        if (this.taxonTaskPool.getNumTaxon() != tree.getExternalNodeCount()) {
            throw new IllegalArgumentException("Incorrectly specified TaxonTaskPool");
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
    public Parameter getParameter() { return factorAnalysisLikelihood.getLoadings(); }

    @Override
    public int getDimension() {
        return dimFactors * dimTrait;
    }

    private ReadableMatrix shiftToSecondMoment(WrappedMatrix variance, ReadableVector mean) {

        assert(variance.getMajorDim() == variance.getMinorDim());
        assert(variance.getMajorDim()== mean.getDim());

        final int dim = variance.getMajorDim();

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                variance.set(i,j, variance.get(i,j) + mean.get(i) * mean.get(j));
            }
        }

        return variance;
    }

    private WrappedNormalSufficientStatistics getWeightedAverage(ReadableVector m1, ReadableMatrix p1,
                                                                 ReadableVector m2, ReadableMatrix p2) {

        assert (m1.getDim() == m2.getDim());
        assert (p1.getDim() == p2.getDim());

        assert (m1.getDim() == p1.getMinorDim());
        assert (m1.getDim() == p1.getMajorDim());

        final WrappedVector m12 = new WrappedVector.Raw(new double[m1.getDim()], 0, dimFactors);
        final DenseMatrix64F p12 = new DenseMatrix64F(dimFactors, dimFactors);
        final DenseMatrix64F v12 = new DenseMatrix64F(dimFactors, dimFactors);

        final WrappedMatrix wP12 = new WrappedMatrix.WrappedDenseMatrix(p12);
        final WrappedMatrix wV12 = new WrappedMatrix.WrappedDenseMatrix(v12);

        MissingOps.add(p1, p2, wP12);
        safeInvert2(p12, v12, false);

        weightedAverage(m1, p1, m2, p2, m12, wV12, dimFactors);

        return new WrappedNormalSufficientStatistics(m12, wP12, wV12);
    }

    @Override
    public double[] getGradientLogDensity() {

        final double[][] gradients = new double[taxonTaskPool.getNumThreads()][getDimension()];

        final ReadableVector gamma = new WrappedVector.Parameter(factorAnalysisLikelihood.getPrecision());
        final ReadableMatrix loadings = ReadableMatrix.Utils.transposeProxy(
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
        
        final List<WrappedNormalSufficientStatistics> allStatistics =
                fullConditionalDensity.getTrait(tree, null); // TODO Need to test if faster to load inside loop

        assert (allStatistics.size() == tree.getExternalNodeCount());

        taxonTaskPool.fork(new TaxonTaskPool.TaxonCallable() {
            @Override
            public void execute(int taxon, int thread) {
                computeGradientForOneTaxon(thread, taxon, loadings, gamma, allStatistics.get(taxon), gradients);
            }
        });


        return join(gradients);
    }


    private void computeGradientForOneTaxon(final int index,
                                            final int taxon,
                                            final ReadableMatrix loadings,
                                            final ReadableVector gamma,
                                            final WrappedNormalSufficientStatistics statistic,
                                            final double[][] gradArray) {

        final WrappedVector y = getTipData(taxon);
        final WrappedNormalSufficientStatistics dataKernel = getTipKernel(taxon);

        final ReadableVector meanKernel = dataKernel.getMean();
        final ReadableMatrix precisionKernel = dataKernel.getPrecision();

        if (DEBUG) {
            System.err.println("Y " + taxon + " : " + y);
            System.err.println("YM" + taxon + " : " + meanKernel);
            System.err.println("YP" + taxon + " : " + precisionKernel);
        }

//        final List<WrappedNormalSufficientStatistics> statistics =
//                fullConditionalDensity.getTrait(tree, tree.getExternalNode(taxon)); // TODO Suspect faster here

//        for (WrappedNormalSufficientStatistics statistic : statistics) {  // TODO Maybe need to re-enable

            final ReadableVector meanFactor = statistic.getMean();
            final WrappedMatrix precisionFactor = statistic.getPrecision();
            final WrappedMatrix varianceFactor = statistic.getVariance();

            if (DEBUG) {
                System.err.println("FM" + taxon + " : " + meanFactor);
                System.err.println("FP" + taxon + " : " + precisionFactor);
                System.err.println("FV" + taxon + " : " + varianceFactor);
            }

            final WrappedNormalSufficientStatistics convolution = getWeightedAverage(
                    meanFactor, precisionFactor,
                    meanKernel, precisionKernel);

            final ReadableVector mean = convolution.getMean();
            final ReadableMatrix precision = convolution.getPrecision();
            final WrappedMatrix variance = convolution.getVariance();

            if (DEBUG) {
                System.err.println("CM" + taxon + " : " + mean);
                System.err.println("CP" + taxon + " : " + precision);
                System.err.println("CV" + taxon + " : " + variance);
            }

            final ReadableMatrix secondMoment = shiftToSecondMoment(variance, mean);
            final ReadableMatrix product = ReadableMatrix.Utils.productProxy(
                    secondMoment, loadings
            );

            if (DEBUG) {
                System.err.println("S" + taxon + " : " + secondMoment);
                System.err.println("P" + taxon + " : " + product);
            }

            for (int factor = 0; factor < dimFactors; ++factor) {
                for (int trait = 0; trait < dimTrait; ++trait) {
                    if (!missing[taxon * dimTrait + trait]) {
                        gradArray[index][factor * dimTrait + trait] +=
                                (mean.get(factor) * y.get(trait) - product.get(factor, trait))
                                        * gamma.get(trait);

                    }
                }
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

    @Override
    public String getReport() {
        return getReport(getGradientLogDensity());
    }

    private WrappedVector getTipData(int taxonIndex) {
        return new WrappedVector.Parameter(data, taxonIndex * dimTrait, dimTrait);
    }

    private WrappedNormalSufficientStatistics getTipKernel(int taxonIndex) {
        double[] buffer = factorAnalysisLikelihood.getTipPartial(taxonIndex, false);
        return new WrappedNormalSufficientStatistics(buffer, 0, dimFactors, null, PrecisionType.FULL);
    }

    private MultivariateFunction numeric = new MultivariateFunction() {

        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                factorAnalysisLikelihood.getLoadings().setParameterValue(i, argument[i]);
            }

            treeDataLikelihood.makeDirty();
            factorAnalysisLikelihood.makeDirty();

            return treeDataLikelihood.getLogLikelihood() + factorAnalysisLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return factorAnalysisLikelihood.getLoadings().getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return Double.NEGATIVE_INFINITY;
        }

        @Override
        public double getUpperBound(int n) {
            return Double.POSITIVE_INFINITY;
        }
    };

    private String getReport(double[] gradient) {

        String result = new WrappedVector.Raw(gradient).toString();

        if (NUMERICAL_CHECK) {

            Parameter loadings = factorAnalysisLikelihood.getLoadings();
            double[] savedValues = loadings.getParameterValues();
            double[] testGradient = NumericalDerivative.gradient(numeric, loadings.getParameterValues());
            for (int i = 0; i < savedValues.length; ++i) {
                loadings.setParameterValue(i, savedValues[i]);
            }

            result += "\nNumerical estimate: \n" +
                    new WrappedVector.Raw(testGradient) +
                    " @ " + new WrappedVector.Raw(loadings.getParameterValues()) + "\n";
        }

        return result;
    }

    private static final boolean DEBUG = false;
    private static final boolean NUMERICAL_CHECK = true;

    private static final String PARSER_NAME = "integratedFactorAnalysisLoadingsGradient";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood)
                    xo.getChild(TreeDataLikelihood.class);

            IntegratedFactorAnalysisLikelihood factorAnalysis = (IntegratedFactorAnalysisLikelihood)
                    xo.getChild(IntegratedFactorAnalysisLikelihood.class);

            DataLikelihoodDelegate likelihoodDelegate = treeDataLikelihood.getDataLikelihoodDelegate();

            if (!(likelihoodDelegate instanceof ContinuousDataLikelihoodDelegate)) {
                throw new XMLParseException("TODO");
            }

            ContinuousDataLikelihoodDelegate continuousDataLikelihoodDelegate =
                    (ContinuousDataLikelihoodDelegate) likelihoodDelegate;

            TaxonTaskPool taxonTaskPool = (TaxonTaskPool) xo.getChild(TaxonTaskPool.class);

            // TODO Check dimensions, parameters, etc.

            return new IntegratedLoadingsGradient(
                    treeDataLikelihood,
                    continuousDataLikelihoodDelegate,
                    factorAnalysis,
                    taxonTaskPool);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Generates a gradient provider for the loadings matrix when factors are integrated out";
        }

        @Override
        public Class getReturnType() {
            return IntegratedLoadingsGradient.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(IntegratedFactorAnalysisLikelihood.class),
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(TaxonTaskPool.class, true),
        };
    };
}
