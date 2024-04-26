package dr.inference.model;

import dr.util.Transform;


public class TransformedMatrixParameter extends TransformedMultivariateParameter implements MatrixParameterInterface {

    private final int rowDim;
    private final int colDim;

    public TransformedMatrixParameter(Parameter parameter, Transform.MultivariableTransform transform,
                                      Boolean inverse, int rowDim, int colDim) {
        super(parameter, transform, inverse);
        this.rowDim = rowDim;
        this.colDim = colDim;
    }


    public TransformedMatrixParameter(MatrixParameterInterface parameter, Transform.MultivariableTransform transform,
                                      Boolean inverse) {
        this(parameter, transform, inverse, parameter.getRowDimension(), parameter.getColumnDimension());
    }


    @Override
    public double getParameterValue(int row, int col) {
        return getParameterValue(index(row, col));
    }

    @Override
    public Parameter getParameter(int column) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        setParameterValue(index(row, col), value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        setParameterValueQuietly(index(row, col), value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        setParameterValueNotifyChangedAll(index(row, col), value);
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double[][] getParameterAsMatrix() {
        return MatrixParameterInterface.getParameterAsMatrix(this);
    }

    @Override
    public int getColumnDimension() {
        return colDim;
    }

    @Override
    public int getRowDimension() {
        return rowDim;
    }

    @Override
    public int getUniqueParameterCount() {
        return 1;
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return parameter;
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

}


