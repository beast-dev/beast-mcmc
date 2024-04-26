package dr.inference.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.geo.distributions.MatrixVonMisesFisherDistribution;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

public class MatrixVonMisesFisherGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String MATRIX_VONMISES_FISHER_GIBBS = "matrixVonMisesFisherGibbsOperator";
    private final FactorAnalysisOperatorAdaptor adaptor;
    private final DenseMatrix64F C;
    private final MatrixVonMisesFisherDistribution vonMisesFisherDistribution;

    MatrixVonMisesFisherGibbsOperator(double weight, FactorAnalysisOperatorAdaptor adaptor) {
        setWeight(weight);
        this.adaptor = adaptor;

        this.C = new DenseMatrix64F(adaptor.getNumberOfTraits(), adaptor.getNumberOfFactors());
        this.vonMisesFisherDistribution =
                new MatrixVonMisesFisherDistribution(adaptor);

    }

    @Override
    public String getOperatorName() {
        return MATRIX_VONMISES_FISHER_GIBBS;
    }

    @Override
    public double doOperation() {

        double[] draw = vonMisesFisherDistribution.nextRandom();

        double[] kBuffer = new double[adaptor.getNumberOfFactors()];

        for (int i = 0; i < adaptor.getNumberOfTraits(); i++) {
            int offset = i;
            for (int j = 0; j < adaptor.getNumberOfFactors(); j++) {
                kBuffer[j] = draw[offset];
                offset += adaptor.getNumberOfTraits();
            }
            adaptor.setLoadingsForTraitQuietly(i, kBuffer);
        }
        adaptor.fireLoadingsChanged();
        return 0;
    }


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            IntegratedFactorAnalysisLikelihood factorLikelihood =
                    (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
            FactorAnalysisOperatorAdaptor.IntegratedFactors adaptor =
                    new FactorAnalysisOperatorAdaptor.IntegratedFactors(factorLikelihood, treeLikelihood);
            double weight = xo.getDoubleAttribute(WEIGHT);
            return new MatrixVonMisesFisherGibbsOperator(weight, adaptor);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(IntegratedFactorAnalysisLikelihood.class),
                    AttributeRule.newDoubleRule(WEIGHT)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return MatrixVonMisesFisherGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return MATRIX_VONMISES_FISHER_GIBBS;
        }
    };
}
