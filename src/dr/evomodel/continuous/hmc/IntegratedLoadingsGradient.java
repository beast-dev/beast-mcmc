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
import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert;
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
    private final List<Integer> missingIndices;

    private IntegratedLoadingsGradient(TreeDataLikelihood treeDataLikelihood,
                               ContinuousDataLikelihoodDelegate likelihoodDelegate,
                               IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.factorAnalysisLikelihood = factorAnalysisLikelihood;

        String traitName = factorAnalysisLikelihood.getModelName(); // TODO Is this correct?

        String fcdName = WrappedTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addWrappedFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
        this.tree = treeDataLikelihood.getTree();

        this.dimTrait = factorAnalysisLikelihood.getTraitDimension();
        this.dimFactors = factorAnalysisLikelihood.getNumberOfFactors();

        this.data = factorAnalysisLikelihood.getParameter();
        this.missingIndices = factorAnalysisLikelihood.getMissingIndices();

        List<Likelihood> likelihoodList = new ArrayList<Likelihood>();
        likelihoodList.add(treeDataLikelihood);
        likelihoodList.add(factorAnalysisLikelihood);
        this.likelihood = new CompoundLikelihood(likelihoodList);

    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return factorAnalysisLikelihood.getLoadings(); // TODO May need to work with vech(L)
    }

    @Override
    public int getDimension() {
        return dimFactors * dimTrait; // TODO May need to work with vech(L)
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

        final WrappedVector m12 = new WrappedVector.Raw(new double[m1.getDim()], 0, dimTrait);
        final DenseMatrix64F p12 = new DenseMatrix64F(dimFactors, dimFactors);
        final DenseMatrix64F v12 = new DenseMatrix64F(dimFactors, dimFactors);

        final WrappedMatrix wP12 = new WrappedMatrix.WrappedDenseMatrix(p12);
        final WrappedMatrix wV12 = new WrappedMatrix.WrappedDenseMatrix(v12);

        MissingOps.add(p1, p2, wP12);
        safeInvert(p12, v12, false);

        weightedAverage(m1, p1, m2, p2, m12, wV12, dimFactors);

        return new WrappedNormalSufficientStatistics(m12, wP12, wV12);
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradient = new double[getDimension()];

        ReadableVector gamma = new WrappedVector.Parameter(factorAnalysisLikelihood.getPrecision());
        ReadableMatrix loadings = new WrappedMatrix.MatrixParameter(factorAnalysisLikelihood.getLoadings());

        if (DEBUG) {
            System.err.println("G : " + gamma);
            System.err.println("L : " + loadings);
        }

        assert (gamma.getDim() == dimTrait);
        assert (loadings.getMajorDim() == dimFactors);
        assert (loadings.getMinorDim() == dimTrait);

        // [E(F) Y^t - E(FF^t)L]\Gamma
        // E(FF^t) = V(F) - E(F)E(F)^t
        // error in comment ^ could indicate error in implementation

        // Y: N x P
        // F: N x K
        // L: K x P

        // (K x N)(N x P)(P x P) - (K x N)(N x K)(K x P)(P x P)
        // sum_{N} (K x 1)(1 x P)(P x P) - (K x K)(K x P)(P x P)
        
        for (int taxon = 0; taxon < tree.getExternalNodeCount(); ++taxon) {

            WrappedVector y = getTipData(taxon);
            WrappedNormalSufficientStatistics dataKernel = getTipKernel(taxon);

            ReadableVector meanKernel = dataKernel.getMean();
            ReadableMatrix precisionKernel = dataKernel.getPrecision();

            if (DEBUG) {
                System.err.println("Y " + taxon + " : " + y);
                System.err.println("YM" + taxon + " : " + meanKernel);
                System.err.println("YP" + taxon + " : " + precisionKernel);
            }

            // TODO Work with fullConditionalDensity
            List<WrappedNormalSufficientStatistics> statistics =
                    fullConditionalDensity.getTrait(tree, tree.getExternalNode(taxon));

            for (WrappedNormalSufficientStatistics statistic : statistics) {

                ReadableVector meanFactor = statistic.getMean();
                WrappedMatrix precisionFactor = statistic.getPrecision();
                WrappedMatrix varianceFactor = statistic.getVariance();

                if (DEBUG) {
                    System.err.println("FM" + taxon + " : " + meanFactor);
                    System.err.println("FP" + taxon + " : " + precisionFactor);
                    System.err.println("FV" + taxon + " : " + varianceFactor);
                }

                WrappedNormalSufficientStatistics convolution = getWeightedAverage(
                        meanFactor, precisionFactor,
                        meanKernel, precisionKernel);

                ReadableVector mean = convolution.getMean();
                ReadableMatrix precision = convolution.getPrecision();
                WrappedMatrix variance = convolution.getVariance();

                if (DEBUG) {
                    System.err.println("CM" + taxon + " : " + mean);
                    System.err.println("CP" + taxon + " : " + precision);
                    System.err.println("CV" + taxon + " : " + variance);
                }

                ReadableMatrix secondMoment = shiftToSecondMoment(variance, mean);
                ReadableMatrix product = ReadableMatrix.Utils.productProxy(secondMoment,loadings);

                if (DEBUG) {
                    System.err.println("S" + taxon + " : " + secondMoment);
                    System.err.println("P" + taxon + " : " + product);
                }

                double[] contribution;

                if (DEBUG) {
                        contribution = new double[dimTrait * dimFactors];
                }

                int offset = 0;
                for (int factor = 0; factor < dimFactors; ++factor) {
                    for (int trait = 0; trait < dimTrait; ++trait) {

                        // TODO Handle missing values with ...
                        missingIndices.contains(offset);

                        if (DEBUG) {
                            contribution[factor * dimTrait + trait] =
                                    (mean.get(factor) * y.get(trait) - product.get(factor, trait))
                                            * gamma.get(trait);
                        }


                        gradient[factor * dimTrait + trait] +=
                                (mean.get(factor) * y.get(trait) - product.get(factor, trait))
                                        * gamma.get(trait);

                        ++offset;
                    }
                }

                if (DEBUG) {
                    System.err.println("C" + taxon + " : " + new WrappedVector.Raw(contribution));
                    System.err.println();
                }

            }
        }

        if (DEBUG) {
            System.err.println(getReport(gradient));
        }

        return gradient;
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
        return new WrappedNormalSufficientStatistics(buffer, 0, dimTrait, null, PrecisionType.FULL);
    }

    private MultivariateFunction numeric = new MultivariateFunction() {

        // TODO Handle vech(loadings)
        // TODO Transform each parameter into (-\infty,\infty)

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

         if (DEBUG) {

             Parameter loadings = factorAnalysisLikelihood.getLoadings();
             double[] savedValues = loadings.getParameterValues();
             double[] testGradient = NumericalDerivative.gradient(numeric, loadings.getParameterValues());
             for (int i = 0; i < savedValues.length; ++i) {
                 loadings.setParameterValue(i, savedValues[i]);
             }

             result += "\nDebug info: \n" +
                     new WrappedVector.Raw(testGradient) +
                     " @ " + new WrappedVector.Raw(loadings.getParameterValues()) + "\n";
         }

         return result;
    }

    private static final boolean DEBUG = true;

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

            // TODO Check dimensions, parameters, etc.

            return new IntegratedLoadingsGradient(
                    treeDataLikelihood,
                    continuousDataLikelihoodDelegate,
                    factorAnalysis);
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
        };
    };
}
