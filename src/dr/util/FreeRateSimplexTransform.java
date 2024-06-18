package dr.util;

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

public class FreeRateSimplexTransform extends Transform.MultivariateTransform {

    private Parameter weights;

    public FreeRateSimplexTransform(int dim, Parameter weights) {
        super(dim);
        this.weights = weights;
        this.outputDimension = dim + 1;
    }

    public FreeRateSimplexTransform(int dim) {
        super(dim);
        weights = new Parameter.Default(dim, (double) 1 /dim);
        this.outputDimension = dim + 1;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public String getTransformName() {
        return "freeRateSimplex";
    }


    @Override
    protected double[] transform(double[] values) {

        // This is a transformation of an n-dimensional vector to another n-dimensional vector but what comes _out_ is a
        // simplex of one greater dimension. The weights parameter has dimension n+1.

        double[] out = new double[values.length + 1];

        double denominator = 0;

        for(int i=0; i<dim; i++){
            denominator +=  values[i]*weights.getParameterValue(i);
        }

        denominator += weights.getParameterValue(dim);

        for(int i=0; i<dim; i++){

            out[i] = values[i]/denominator;

        }

        out[dim] = 1/denominator;
        return(out);

    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        Matrix partialsMatrix = new Matrix(dim, dim);

        double sqrtDenominator;

        // reminder: dim is one less than the dimension of the simplex

        double tempSum = 0;

        for(int i=0; i<dim; i++){
            tempSum += weights.getParameterValue(i)*values[i];
        }

        sqrtDenominator = 1-tempSum;

        assert sqrtDenominator > 0;

        double denominator = sqrtDenominator * sqrtDenominator;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (j == i) {
                    partialsMatrix.set(i, j, weights.getParameterValue(dim)*(sqrtDenominator + values[i]*weights.getParameterValue(i)) / denominator);
                } else {
                    partialsMatrix.set(i, j, (values[i]*weights.getParameterValue(j)*weights.getParameterValue(dim)) / denominator);
                }
            }
        }

        double logJacobian = 0;
        try {
            logJacobian = Math.log(Math.abs(partialsMatrix.determinant()));
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        return logJacobian;
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("not implemented");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        throw new RuntimeException("not implemented");
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String WEIGHTS = "weights";
        public static final String DIMENSION = "dimension";
        public static final String FREERATE_SIMPLEX_TRANSFORM = "freeRateSimplexTransform";
        @Override
        public Class getReturnType() {
            return FreeRateSimplexTransform.class;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            if(xo.hasChildNamed(WEIGHTS)){
                Parameter weights = (Parameter) xo.getElementFirstChild(WEIGHTS);
                return new FreeRateSimplexTransform(weights.getDimension()-1, weights);
            } else if(xo.hasAttribute(DIMENSION)) {
                int dimension = xo.getIntegerAttribute(DIMENSION);
                return new FreeRateSimplexTransform(dimension);
            } else {
                throw new XMLParseException("FreeRateSimplexTransform must have either a dimension attribute or" +
                        "a set of weights");
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new OrRule(
                    AttributeRule.newDoubleRule(DIMENSION),
                    new ElementRule(WEIGHTS, Parameter.class))
            };
        }

        @Override
        public String getParserDescription() {
            return "Transform from a (n-1)-dimensional vector of positive real numbers to a n-dimensional weighted " +
                    "simplex";
        }

        @Override
        public String getParserName() {
            return FREERATE_SIMPLEX_TRANSFORM;
        }


    };

}
