package dr.inference.model;

import dr.util.Transform;
import dr.xml.*;

public class TransformedMatrixParameter extends TransformedMultivariateParameter implements MatrixParameterInterface {

    private final MatrixParameterInterface matrix;
    private final MatrixParameterInterface transformedMatrix;
    private final int nRows;
    private final int nCols;
    private boolean needToUpdate = true;


    public TransformedMatrixParameter(MatrixParameterInterface matrix, Transform.MultivariateTransform transform) {
        super(matrix, transform);
        this.matrix = matrix;
        this.nRows = matrix.getRowDimension();
        this.nCols = matrix.getColumnDimension();

        this.transformedMatrix = new FastMatrixParameter(matrix.getId() + ".transformed", nRows, nCols, 0.0);
        addVariableListener((VariableListener) matrix);
    }


    @Override
    public double getParameterValue(int row, int col) {
        update();
        return transformedMatrix.getParameterValue(row, col);
    }

    @Override
    public Parameter getParameter(int column) {
        update();
        return transformedMatrix.getParameter(column);
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        matrix.setParameterValue(row, col, inverse(value));
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        matrix.setParameterValueQuietly(row, col, inverse(value));
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        matrix.setParameterValueNotifyChangedAll(row, col, inverse(value));
    }

    @Override
    public double[] getColumnValues(int col) {
        update();
        return transformedMatrix.getColumnValues(col);
    }

    @Override
    public double[][] getParameterAsMatrix() {
        update();
        return transformedMatrix.getParameterAsMatrix();
    }

    @Override
    public int getColumnDimension() {
        return matrix.getColumnDimension();
    }

    @Override
    public int getRowDimension() {
        return matrix.getRowDimension();
    }

    @Override
    public int getUniqueParameterCount() {
        return transformedMatrix.getUniqueParameterCount();
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return transformedMatrix.getUniqueParameter(index);
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        update();
        transformedMatrix.copyParameterValues(destination, offset);
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        update();
        matrix.setAllParameterValuesQuietly(inverse(values), offset);
    }

    @Override
    public String toSymmetricString() {
        return transformedMatrix.toSymmetricString();
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return transformedMatrix.isConstrainedSymmetric(); //TODO: get from transform
    }

    @Override
    protected void update() {
        if (needToUpdate) {
            unTransformedValues = matrix.getParameterValues();
            transformedValues = transform(unTransformedValues);
            transformedMatrix.setAllParameterValuesQuietly(transformedValues, 0);
            transformedMatrix.fireParameterChangedEvent();
            needToUpdate = false;
        }
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        needToUpdate = true;
        update(); //TODO: try to do this lazily
        super.variableChangedEvent(variable, index, type);

    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private final String TRANSFORMED_MATRIX = "transformedMatrixParameter";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameterInterface param = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            Transform.MultivariateTransform transform =
                    (Transform.MultivariateTransform) xo.getChild(Transform.MultivariateTransform.class);
            return new TransformedMatrixParameter(param, transform);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
                    new ElementRule(Transform.MultivariateTransform.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Transformed matrix";
        }

        @Override
        public Class getReturnType() {
            return TransformedMultivariateParameter.class;
        }

        @Override
        public String getParserName() {
            return TRANSFORMED_MATRIX;
        }
    };
}
