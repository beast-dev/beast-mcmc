package dr.inference.model;

/**
 * @author Marc Suchard
 */
public class DiagonalMatrix extends MatrixParameter {

	private Parameter diagonalParameter;

	public DiagonalMatrix(Parameter param) {
		super(MATRIX_PARAMETER);
		addParameter(param);
		diagonalParameter = param;
	}

//	public DiagonalMatrix(String name, Parameter parameter) {
//		Parameter.Default(name, parameters);
//	}

	public double getParameterValue(int row, int col) {
		if (row != col)
			return 0.0;
		return diagonalParameter.getParameterValue(row);
	}

	public double[][] getParameterAsMatrix() {
		final int I = getDimension();
		double[][] parameterAsMatrix = new double[I][I];
		for (int i = 0; i < I; i++) {
			parameterAsMatrix[i][i] = diagonalParameter.getParameterValue(i);
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
