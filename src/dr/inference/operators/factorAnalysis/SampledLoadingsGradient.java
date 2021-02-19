package dr.inference.operators.factorAnalysis;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;


public class SampledLoadingsGradient implements GradientWrtParameterProvider, ModelListener {

    private final MatrixParameterInterface loadings;
    private final LatentFactorModel factorModel;
    private final NewLoadingsGibbsOperator gibbsOperator;
    private final double[][] means;
    private final double[][][] precisions;
    private final int nFactors;
    private final int nTraits;
    private boolean statisticsKnown = false;

    SampledLoadingsGradient(LatentFactorModel factorModel, NewLoadingsGibbsOperator gibbsOperator) {

        this.factorModel = factorModel;
        this.loadings = factorModel.getLoadings();
        this.gibbsOperator = gibbsOperator;
        gibbsOperator.setStatisticsOnly(true);
        this.nFactors = loadings.getColumnDimension();
        this.nTraits = loadings.getRowDimension();
        this.means = new double[nTraits][nFactors];
        this.precisions = new double[nTraits][nFactors][nFactors];

        factorModel.addModelListener(this);
    }

    @Override
    public Likelihood getLikelihood() {
        return factorModel;
    }

    @Override
    public Parameter getParameter() {
        return loadings;
    }

    @Override
    public int getDimension() {
        return loadings.getDimension();
    }

    private void updateStatistics() {
        gibbsOperator.getAdaptor().drawFactors();
        for (int i = 0; i < nTraits; i++) {
            gibbsOperator.drawI(i);
            double[] mean = gibbsOperator.getCurrentMean(i);
            double[][] prec = gibbsOperator.getCurrentPrecision(i);
            System.arraycopy(mean, 0, means[i], 0, nFactors);
            for (int j = 0; j < nFactors; j++) {
                System.arraycopy(prec[j], 0, precisions[i][j], 0, nFactors);
            }
        }

        statisticsKnown = true;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        if (object != loadings) {
            statisticsKnown = false;
        }
    }

    @Override
    public void modelRestored(Model model) {
        statisticsKnown = false;
    }


    @Override
    public double[] getGradientLogDensity() {
        if (!statisticsKnown) {
            updateStatistics();
        }
        double[] gradient = new double[getDimension()];
        double error[] = new double[loadings.getColumnDimension()];

        for (int i = 0; i < nTraits; i++) {
            double[] mean = means[i];
            double[][] prec = precisions[i];
            for (int j = 0; j < nFactors; j++) {
                error[j] = loadings.getParameterValue(i, j) - mean[j];
            }
            Matrix pMat = new Matrix(prec);
            Vector v = null;
            try {
                v = pMat.product(new Vector(error));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

            for (int j = 0; j < nFactors; j++) {
                gradient[j * nTraits + i] = -v.component(j);
            }

        }
        return gradient;
    }


    private static final String SAMPLED_LOADINGS_GRADIENT = "sampledLoadingsGradient";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            LatentFactorModel factorModel = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
            NewLoadingsGibbsOperator gibbsOperator = (NewLoadingsGibbsOperator) xo.getChild(NewLoadingsGibbsOperator.class);
            return new SampledLoadingsGradient(factorModel, gibbsOperator);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(NewLoadingsGibbsOperator.class),
                    new ElementRule(LatentFactorModel.class)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return SampledLoadingsGradient.class;
        }

        @Override
        public String getParserName() {
            return SAMPLED_LOADINGS_GRADIENT;
        }
    };


}
