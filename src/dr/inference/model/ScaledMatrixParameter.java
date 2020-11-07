package dr.inference.model;

import dr.xml.*;

public class ScaledMatrixParameter extends Parameter.Abstract implements MatrixParameterInterface, VariableListener {

    private final MatrixParameterInterface matrixParameter;
    private final Parameter scaleParameter;
    private final double[][] columnBuffers;
    private final double[][] storedColumnBuffers;
    private Boolean valuesKnown = false;
    private Boolean storedValuesKnown = false;
    private final String name;
    private final boolean[] renormalize;
    private boolean doNotPropogateChangesUp = false;

    public ScaledMatrixParameter(String name, MatrixParameterInterface matrixParameter, Parameter scaleParameter) {
        this.name = name;

        this.matrixParameter = matrixParameter;
        this.scaleParameter = scaleParameter;
        this.columnBuffers = new double[matrixParameter.getColumnDimension()][matrixParameter.getRowDimension()];
        this.storedColumnBuffers = new double[matrixParameter.getColumnDimension()][matrixParameter.getRowDimension()];

        this.renormalize = new boolean[scaleParameter.getDimension()];

        scaleParameter.addParameterListener(this);
        matrixParameter.addParameterListener(this);
    }

    private void updateBuffer() {
        if (!valuesKnown) {
            for (int col = 0; col < matrixParameter.getColumnDimension(); col++) {
                for (int row = 0; row < matrixParameter.getRowDimension(); row++) {

                    columnBuffers[col][row] = matrixParameter.getParameterValue(row, col) *
                            scaleParameter.getParameterValue(col);
                }
            }

            valuesKnown = true;
        }
    }

    public MatrixParameterInterface getMatrixParameter() {
        return matrixParameter;
    }

    public Parameter getScaleParameter() {
        return scaleParameter;
    }

    @Override
    public double getParameterValue(int row, int col) {
        updateBuffer();
        return columnBuffers[col][row];
    }

    @Override
    public Parameter getParameter(int column) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double getParameterValue(int dim) {
        updateBuffer();
        int col = dim / matrixParameter.getRowDimension();
        int row = dim - col * matrixParameter.getRowDimension();

        return columnBuffers[col][row];

    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public String getParameterName() {
        return name;
    }

    @Override
    public void addBounds(Bounds<Double> bounds) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public Bounds<Double> getBounds() {
        throw new RuntimeException("not implemented");
    }

    @Override
    public void addDimension(int index, double value) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double removeDimension(int index) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        matrixParameter.setParameterValueQuietly(row, col, value);
        renormalize[col] = true;
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double[] getColumnValues(int col) {
        updateBuffer();
        return columnBuffers[col].clone();
    }

    @Override
    public double[][] getParameterAsMatrix() {
        updateBuffer();
        return columnBuffers.clone();
    }

    @Override
    public int getDimension() {
        return matrixParameter.getDimension();
    }

    @Override
    public int getColumnDimension() {
        return matrixParameter.getColumnDimension();
    }

    @Override
    public int getRowDimension() {
        return matrixParameter.getRowDimension();
    }

    @Override
    public int getUniqueParameterCount() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        updateBuffer();
        int newOffset = offset;
        for (int i = 0; i < matrixParameter.getColumnDimension(); i++) {
            System.arraycopy(columnBuffers[i], 0, destination, newOffset, matrixParameter.getRowDimension());
            newOffset += matrixParameter.getRowDimension();
        }
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        if (!doNotPropogateChangesUp) fireParameterChangedEvent();

        valuesKnown = false; //TODO: do only update necessary indices
    }

    @Override
    public void fireParameterChangedEvent() {
        for (int col = 0; col < renormalize.length; col++) {
            if (renormalize[col]) {
                doNotPropogateChangesUp = true;
                double norm = 0.0;
                for (int row = 0; row < matrixParameter.getRowDimension(); row++) {
                    double val = matrixParameter.getParameterValue(row, col);
                    norm += val * val;
                }

                norm = Math.sqrt(norm);

                scaleParameter.setParameterValueQuietly(col, norm);
                for (int row = 0; row < matrixParameter.getRowDimension(); row++) {
                    matrixParameter.setParameterValue(row, col, matrixParameter.getParameterValue(row, col) / norm);
                }
                renormalize[col] = false;
            }
        }
        if (doNotPropogateChangesUp) {
            scaleParameter.fireParameterChangedEvent();
            matrixParameter.fireParameterChangedEvent();
            doNotPropogateChangesUp = false;
        }

        fireParameterChangedEvent(-1, Parameter.ChangeType.VALUE_CHANGED);
    }

    private void transferBuffer(double[][] src, double[][] dest) {
        for (int i = 0; i < src.length; i++) {
            System.arraycopy(src[i], 0, dest[i], 0, src[i].length);
        }
    }

    @Override
    protected void storeValues() {
        scaleParameter.storeParameterValues();
        matrixParameter.storeParameterValues();

        if (valuesKnown) {
            transferBuffer(columnBuffers, storedColumnBuffers);
        }
        storedValuesKnown = valuesKnown;
    }

    @Override
    protected void restoreValues() {
        scaleParameter.restoreParameterValues();
        matrixParameter.restoreParameterValues();

        if (storedValuesKnown) {
            transferBuffer(storedColumnBuffers, columnBuffers);
        }
        valuesKnown = storedValuesKnown;
    }

    @Override
    protected void acceptValues() {
        scaleParameter.acceptParameterValues();
        matrixParameter.acceptParameterValues();
    }

    @Override
    public String getDimensionName(int dim) {
        int column = dim / matrixParameter.getRowDimension();
        int row = dim - (column * matrixParameter.getRowDimension());
        return getStatisticName() + "." + (column + 1) + "." + (row + 1);
    }

    @Override
    protected void adoptValues(Parameter source) {
        throw new RuntimeException("not implemented");
    }

    public static final String SCALED_MATRIX = "scaledMatrixParameter";
    public static final String SCALE = "scale";
    public static final String MATRIX = "matrix";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter scale = (Parameter) xo.getChild(SCALE).getChild(Parameter.class);
            MatrixParameterInterface matrix =
                    (MatrixParameterInterface) xo.getChild(MATRIX).getChild(MatrixParameterInterface.class);

            final String name = xo.hasId() ? xo.getId() : null;


            return new ScaledMatrixParameter(name, matrix, scale);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(SCALE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(MATRIX, new XMLSyntaxRule[]{
                            new ElementRule(MatrixParameterInterface.class)
                    })
            };
        }

        @Override
        public String getParserDescription() {
            return "parameter that scales each column of a matrix";
        }

        @Override
        public Class getReturnType() {
            return ScaledMatrixParameter.class;
        }

        @Override
        public String getParserName() {
            return SCALED_MATRIX;
        }
    };
}
