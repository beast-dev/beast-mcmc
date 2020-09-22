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
    private final DenseMatrix64F D;
    private final DenseMatrix64F V;
    private final DenseMatrix64F F;
    private final DenseMatrix64F Y;
    private final DenseMatrix64F C;
    private final MatrixVonMisesFisherDistribution vonMisesFisherDistribution;

    MatrixVonMisesFisherGibbsOperator(double weight, FactorAnalysisOperatorAdaptor adaptor) {
        setWeight(weight);
        this.adaptor = adaptor;
        this.V = new DenseMatrix64F(adaptor.getNumberOfTraits(), adaptor.getNumberOfFactors());
        this.D = new DenseMatrix64F(adaptor.getNumberOfFactors(), adaptor.getNumberOfFactors());
        this.F = new DenseMatrix64F(adaptor.getNumberOfTaxa(), adaptor.getNumberOfFactors());
        this.Y = new DenseMatrix64F(adaptor.getNumberOfTaxa(), adaptor.getNumberOfTraits());
        this.C = new DenseMatrix64F(adaptor.getNumberOfTraits(), adaptor.getNumberOfFactors());
        this.vonMisesFisherDistribution =
                new MatrixVonMisesFisherDistribution(adaptor.getNumberOfTraits(), adaptor.getNumberOfFactors());

    }

    @Override
    public String getOperatorName() {
        return MATRIX_VONMISES_FISHER_GIBBS;
    }

    @Override
    public double doOperation() {
        splitLoadings();
        fillFactors();
        fillTraits();
        double maxPrecision = getMaximumPrecision();


        CommonOps.multTransA(F, Y, C);
        for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
            double scaledNorm = D.get(i, i) * maxPrecision;
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                C.set(j, i, C.get(j, i) * scaledNorm);
            }
        }

        vonMisesFisherDistribution.setC(C.getData());
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


    private void splitLoadings() {
        int offset = 0;
        for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
            double sumSquares = 0;
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                sumSquares += adaptor.getLoadingsValue(offset);
                offset++;
            }
            double norm = Math.sqrt(sumSquares);
            D.set(i, i, norm);
            double invNorm = 1.0 / norm;
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                V.set(offset, adaptor.getLoadingsValue(offset) * invNorm);
            }
        }
    }

    private void fillFactors() {
        adaptor.drawFactors();
        for (int i = 0; i < adaptor.getNumberOfTaxa(); i++) {
            for (int j = 0; j < adaptor.getNumberOfFactors(); j++) {
                F.set(i, j, adaptor.getFactorValue(j, i));
            }
        }
    }

    private void fillTraits() {
        //TODO: fill in missing traits
        for (int i = 0; i < adaptor.getNumberOfTaxa(); i++) {
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                Y.set(i, j, adaptor.getDataValue(j, i));
            }
        }
    }

    private double getMaximumPrecision() {
        double maxPrec = 0;
        for (int i = 0; i < adaptor.getNumberOfTraits(); i++) {
            if (adaptor.getColumnPrecision(i) > maxPrec) {
                maxPrec = adaptor.getColumnPrecision(i);
            }
        }
        return maxPrec;
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
