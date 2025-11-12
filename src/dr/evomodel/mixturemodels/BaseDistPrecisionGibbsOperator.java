package dr.evomodel.mixturemodels;

import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.WishartGammalDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.WishartDistribution;
import dr.math.distributions.WishartStatistics;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;


public class BaseDistPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String BASE_DIST_PREC_OPERATOR = "baseDistPrecisionGibbsOperator";
    public static final String WEIGHT = "weight";
    public static final String BASE_DIST_NUM = "baseDistNum";


    private double pathWeight = 1.0;
    private GenPolyaUrnProcessPrior gpuProcess;
    private final Parameter mean;
    private final MatrixParameterInterface precision;
    private final int baseDistNum;
    private final int dim;
    private double priorDf;
    private SymmetricMatrix priorInverseScaleMatrix;
    private WishartGammalDistributionModel priorModel = null;
    private static final boolean DEBUG = false;
    private double numberObservations;

    public BaseDistPrecisionGibbsOperator(
            GenPolyaUrnProcessPrior gpuProcess,
            WishartStatistics priorDistribution,
            int baseDistNum,
            double weight) {
        super();

        this.gpuProcess = gpuProcess;
        MultivariateNormalDistributionModel density = (MultivariateNormalDistributionModel) gpuProcess.getParametricBaseDist().get(baseDistNum);
        this.mean = density.getMeanParameter();
        this.precision = density.getPrecisionMatrixParameter();
        this.dim = mean.getDimension();
        this.baseDistNum = baseDistNum;

        setupWishartStatistics(priorDistribution);

        setWeight(weight);
    }


    private void setupWishartStatistics(WishartStatistics priorDistribution) {

        this.priorDf = priorDistribution.getDF();
        this.priorInverseScaleMatrix = null;
        double[][] scale = priorDistribution.getScaleMatrix();
        if (scale != null) {
            this.priorInverseScaleMatrix = (SymmetricMatrix) (new SymmetricMatrix(scale)).inverse();
        }
    }

    private void incrementOuterProduct(double[][] S,
                                       GenPolyaUrnProcessPrior gpu) {

        double[] mean = gpu.getParametricBaseDist().get(baseDistNum).getMean();

        CompoundParameter urp = gpu.getUniquelyRealizedParameters();
        int[] isCatActive = gpu.getIsCatActive();

        numberObservations = 0;

        for (int k = 0; k < gpu.maxCategoryCount; k++) {
            if(isCatActive[k] == 1) {
                double[] data = urp.getParameter(k).getParameterValues();
                for (int i = 0; i < dim; i++) {
                    // if(i == 1) {
                    //     System.out.println("category " + k + " dimension " + i + " is: " + data[i]);
                    // }
                    data[i] -= mean[i];
                }

                for (int i = 0; i < dim; i++) {  // symmetric matrix,
                    for (int j = i; j < dim; j++) {
                        S[j][i] = S[i][j] += data[i] * data[j];
                    }
                }
                numberObservations += 1;
            }
        }
    }


    private double[][] getOperationScaleMatrixAndSetObservationCount() {

        double[][] S = new double[dim][dim];
        SymmetricMatrix S2;
        SymmetricMatrix inverseS2 = null;
        numberObservations = 0;

        incrementOuterProduct(S, gpuProcess);

        try {
            S2 = new SymmetricMatrix(S);
            if (pathWeight != 1.0) {
                S2 = (SymmetricMatrix) S2.product(pathWeight);
            }
            if (priorInverseScaleMatrix != null) {
                S2 = priorInverseScaleMatrix.add(S2);
            }
            inverseS2 = (SymmetricMatrix) S2.inverse();

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        assert inverseS2 != null;

        return inverseS2.toComponents();
    }

    public double doOperation() {

        //setupWishartStatistics(priorModel);

        final double[][] scaleMatrix = getOperationScaleMatrixAndSetObservationCount();
        final double treeDf = numberObservations;

        final double df = priorDf + treeDf * pathWeight;

        double[][] draw = WishartDistribution.nextWishart(df, scaleMatrix);

        if (DEBUG) {
            System.err.println("draw = " + new Matrix(draw));
        }

        for (int i = 0; i < dim; i++) {
            Parameter column = precision.getParameter(i);
            for (int j = 0; j < dim; j++)
                column.setParameterValueQuietly(j, draw[j][i]);
        }
        precision.fireParameterChangedEvent();

        return 0;
    }

    public void setPathParameter(double beta) {
        if (beta < 0 || beta > 1) {
            throw new IllegalArgumentException("Illegal path weight of " + beta);
        }
        pathWeight = beta;
    }

    public String getOperatorName() {
        return BASE_DIST_PREC_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return BASE_DIST_PREC_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            int bdNum = xo.getIntegerAttribute(BASE_DIST_NUM);

            //WishartStatistics ws = (WishartStatistics) xo.getChild(WishartStatistics.class);
            MultivariateDistributionLikelihood prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);

            GenPolyaUrnProcessPrior gpu = (GenPolyaUrnProcessPrior) xo.getChild(GenPolyaUrnProcessPrior.class);

            return new BaseDistPrecisionGibbsOperator(gpu, (WishartStatistics) prior.getDistribution(), bdNum, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a Gibbs sampler for the precision matrix of a GPU process multivariate normal base distribution.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newDoubleRule(BASE_DIST_NUM, false),
                //new ElementRule(WishartStatistics.class, false),
                new ElementRule(MultivariateDistributionLikelihood.class, false),
                new ElementRule(GenPolyaUrnProcessPrior.class, false),
        };
    };


}