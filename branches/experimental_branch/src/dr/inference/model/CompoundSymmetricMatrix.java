package dr.inference.model;

/**
 * @author Marc Suchard
 */
public class CompoundSymmetricMatrix extends MatrixParameter {

    private Parameter diagonalParameter;
    private Parameter offDiagonalParameter;

    private boolean asCorrelation = false;

    private int dim;

    public CompoundSymmetricMatrix(Parameter diagonals, Parameter offDiagonal, boolean asCorrelation) {
        super(MATRIX_PARAMETER);
        diagonalParameter = diagonals;
        offDiagonalParameter = offDiagonal;
        addParameter(diagonalParameter);
        addParameter(offDiagonal);
        dim = diagonalParameter.getDimension();
        this.asCorrelation = asCorrelation;
    }

    public double[] getAttributeValue() {
        double[] stats = new double[dim * dim];
        int index = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                stats[index] = getParameterValue(i, j);
                index++;
            }
        }        
        return stats;
    }

    public double[] getDiagonals() {
        return diagonalParameter.getParameterValues();
    }

    public double getOffDiagonal() {
        return offDiagonalParameter.getParameterValue(0);
    }

    public double getParameterValue(int row, int col) {
        if (row != col) {
            if (asCorrelation) {
                return offDiagonalParameter.getParameterValue(0) *
                        Math.sqrt(diagonalParameter.getParameterValue(row) * diagonalParameter.getParameterValue(col));
            }
            return offDiagonalParameter.getParameterValue(0);
        }
        return diagonalParameter.getParameterValue(row);
    }

    public double[][] getParameterAsMatrix() {
        final int I = dim;
        double[][] parameterAsMatrix = new double[I][I];
        for (int i = 0; i < I; i++) {
            parameterAsMatrix[i][i] = getParameterValue(i, i);
            for (int j = i + 1; j < I; j++) {
                parameterAsMatrix[j][i] = parameterAsMatrix[i][j] = getParameterValue(i, j);
            }
        }
        return parameterAsMatrix;
    }

    public int getColumnDimension() {
        return diagonalParameter.getDimension();
    }

    public int getRowDimension() {
        return diagonalParameter.getDimension();
    }

}
