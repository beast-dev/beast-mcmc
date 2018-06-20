package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.evomodel.treedatalikelihood.preorder.WrappedMeanPrecision;
import dr.evomodel.treedatalikelihood.preorder.WrappedTipFullConditionalDistributionDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.continuous.hmc.LinearOrderTreePrecisionTraitProductProvider.castTreeTrait;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrook
 */
public class IntegratedLoadingsGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeTrait<List<WrappedMeanPrecision>> fullConditionalDensity;

    private final TreeDataLikelihood treeDataLikelihood;
    private final IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood;
    private final int dimTrait;
    private final int dimFactors;
    private final Tree tree;
    private final Likelihood likelihood;

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
        return dimTrait * dimTrait; // TODO May need to work with vech(L)
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradient = new double[getDimension()];

        Parameter gamma = factorAnalysisLikelihood.getPrecision();
        assert (gamma.getDimension() == dimTrait);

        // [E(F) Y^t + E(FF^t)L]\Gamma

        // Y: N x P
        // F: N x K
        // L: K x P

        // (K x N)(N x P)(P x P) - (K x N)(N x K)(K x P)(P x P)
        // sum_{N} (K x 1)(1 x P)(P x P) - (K x 1)(1 x K)(K x P)(P x P)
        
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {

            WrappedVector y = getTipData(i);

            // TODO Work with fullConditionalDensity
            List<WrappedMeanPrecision> statistics =
                    fullConditionalDensity.getTrait(tree, tree.getExternalNode(i));

            for (WrappedMeanPrecision meanPrecision : statistics) {
                System.err.println(i + " : " + meanPrecision.getMean());
                System.err.println(i + " : " + meanPrecision.getPrecision());
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
        return new WrappedVector.Raw(factorAnalysisLikelihood.getTipPartial(taxonIndex, false),
                0, dimTrait);
    }

//    private WrappedMeanPrecision getTipData(int tipIndex) {
//        double[] buffer = factorAnalysisLikelihood.getTipPartial(tipIndex, false);
//
//        return new WrappedMeanPrecision(buffer, 0, dimTrait, null, PrecisionType.FULL);
//    }

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
                     " @ " + new WrappedVector.Raw(loadings.getParameterValues());
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
