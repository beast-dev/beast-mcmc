package dr.inference.distribution;

import dr.inference.model.TransformedMatrixParameter;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.util.SVDTransform;
import dr.util.Transform;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class OrthogonalRotationProvider implements NormalStatisticsHelpers.MatrixNormalStatisticsHelper {

    private final TransformedMatrixParameter matrixParameter;
    private final NormalStatisticsHelpers.IndependentNormalStatisticsProvider prior;
    private final SVDTransform transform;
    private final NormalStatisticsHelpers.MatrixNormalStatisticsHelper priorHelper;

    OrthogonalRotationProvider(TransformedMatrixParameter matrixParameter,
                               NormalStatisticsHelpers.IndependentNormalStatisticsProvider prior) {
        this.matrixParameter = matrixParameter;
        this.prior = prior;
        Transform transform = matrixParameter.getTransform();
        if (!(transform instanceof SVDTransform)) {
            throw new RuntimeException("Must be an svdTransform.");
        }
        this.transform = (SVDTransform) transform;
        this.priorHelper = prior.matrixNormalHelper(matrixParameter.getRowDimension(), matrixParameter.getColumnDimension());
    }


    @Override
    public double getScalarPrecision() {
        return priorHelper.getScalarPrecision();
    }

    @Override
    public double[] precisionMeanProduct(int col) {
        return new double[matrixParameter.getColumnDimension()];
    }

    @Override
    public double[][] getColumnPrecision(int col) {
        matrixParameter.getParameterValue(0, 0); //force update if necessary
        double[][] U = transform.getU();
        Matrix uMat = new Matrix(U);
        Matrix originalPrec = new Matrix(priorHelper.getColumnPrecision(col));
        Matrix precMat;
        try {
            Matrix intermediate = uMat.product(originalPrec);
            precMat = intermediate.productWithTransposed(uMat);

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
            throw new RuntimeException("Dimensions don't match.");
//            precMat = null;
        }
        return precMat.toComponents();
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String ORTHO_PROVIDER = "orthogonalRotationProvider";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TransformedMatrixParameter parameter =
                    (TransformedMatrixParameter) xo.getChild(TransformedMatrixParameter.class);
            NormalStatisticsHelpers.IndependentNormalStatisticsProvider prior =
                    (NormalStatisticsHelpers.IndependentNormalStatisticsProvider) xo.getChild(NormalStatisticsHelpers.IndependentNormalStatisticsProvider.class);
            return new OrthogonalRotationProvider(parameter, prior);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[0]; //TODO
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return OrthogonalRotationProvider.class;
        }

        @Override
        public String getParserName() {
            return ORTHO_PROVIDER;
        }
    };
}
