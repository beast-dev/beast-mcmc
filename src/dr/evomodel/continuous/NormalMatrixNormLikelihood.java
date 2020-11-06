package dr.evomodel.continuous;

import dr.inference.distribution.NormalStatisticsProvider;
import dr.inference.model.*;
import dr.inference.operators.repeatedMeasures.MultiplicativeGammaGibbsHelper;
import dr.math.distributions.MultivariateGammaLikelihood;
import dr.util.Transform;
import dr.xml.*;

public class NormalMatrixNormLikelihood extends MultivariateGammaLikelihood implements MultiplicativeGammaGibbsHelper, NormalStatisticsProvider {

    private final int rowDimension;
    private final Parameter columnNorms;

    public NormalMatrixNormLikelihood(String name, Parameter normalPrecision, Parameter columnNorms, int rowDimension) {
        super(name, setupShape(rowDimension, columnNorms.getDimension()), setupScale(normalPrecision), setupSquaredNorms(columnNorms));
        this.columnNorms = columnNorms;

        this.rowDimension = rowDimension;
    }

    private static final Parameter setupShape(int rowDimension, int colDimension) {
        Parameter.Default shape = new Parameter.Default(colDimension, ((double) rowDimension) / 2.0);
        return shape;
    }

    private static final Parameter setupScale(Parameter precision) {
        TransformedParameter variance = new TransformedParameter(precision, new Transform.ReciprocalTransform());
        ScaledParameter scale = new ScaledParameter(2.0, variance);
        return scale;
    }

    private static final Parameter setupSquaredNorms(Parameter norms) {
        TransformedParameter squaredNorms = new TransformedParameter(norms, new Transform.PowerTransform(2.0));
        return squaredNorms;
    }

    @Override
    public double computeSumSquaredErrors(int column) {
        return data.getParameterValue(column);
    }

    @Override
    public Parameter getParameter() {
        return columnNorms;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] grad = super.getGradientLogDensity();
        for (int i = 0; i < dim; i++) {
            grad[i] *= 2 * columnNorms.getParameterValue(i);
        }
        return grad;
    }

    @Override
    public int getRowDimension() {
        return rowDimension;
    }

    @Override
    public int getColumnDimension() {
        return dim;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {


        private static final String NORMAL_NORM_LIKELIHOOD = "normalMatrixNormLikelihood";
        private static final String GLOBAL_PRECISION = "globalPrecision";
        private static final String MATRIX = "matrix";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = xo.hasAttribute(XMLObject.ID) ? xo.getId() : NORMAL_NORM_LIKELIHOOD;
            XMLObject pxo = xo.getChild(GLOBAL_PRECISION);
            Parameter globalPrecision = (Parameter) pxo.getChild(Parameter.class);

            XMLObject mxo = xo.getChild(MATRIX);
            ScaledMatrixParameter matrix = (ScaledMatrixParameter) mxo.getChild(ScaledMatrixParameter.class);

            Parameter scaleParameter = matrix.getScaleParameter();

            if (scaleParameter.getDimension() != globalPrecision.getDimension()) {
                throw new XMLParseException("incompatible dimensions: the " + GLOBAL_PRECISION + " parameter in with id `" +
                        globalPrecision.getId() + "` had dimension " + globalPrecision.getDimension() +
                        ", while the " + MATRIX + " parameter with id `" + matrix.getId() + "` has " +
                        scaleParameter.getDimension() + " columns."
                );
            }

            return new NormalMatrixNormLikelihood(id, globalPrecision, scaleParameter,
                    matrix.getRowDimension());
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(GLOBAL_PRECISION, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(MATRIX, new XMLSyntaxRule[]{
                            new ElementRule(ScaledMatrixParameter.class)
                    })

            };
        }

        @Override
        public String getParserDescription() {
            return "gamma likelihood on the squared norm of the columns of a matrix";
        }

        @Override
        public Class getReturnType() {
            return NormalMatrixNormLikelihood.class;
        }

        @Override
        public String getParserName() {
            return NORMAL_NORM_LIKELIHOOD;
        }
    };

    @Override
    public double getNormalMean(int dim) {
        return 0;
    }

    @Override
    public double getNormalSD(int dim) {
        int row = dim / rowDimension;
        return Math.sqrt(scale.getParameterValue(row) / 2);
    }
}
