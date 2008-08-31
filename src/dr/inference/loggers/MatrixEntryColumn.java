package dr.inference.loggers;

/**
 * @author Marc A. Suchard
 */
public class MatrixEntryColumn extends NumberColumn {

	private int indexI;
	private int indexJ;
	private double[][] mat;

	public MatrixEntryColumn(String name, int indexI, int indexJ, double[][] mat) {
		super(name + "_" + indexI + "/" + indexJ);

		this.mat = mat;

		if (indexI < 0 || indexJ < 0 || indexI >= mat.length || indexJ >= mat[0].length)
			throw new RuntimeException("Out of bounds");

		this.indexI = indexI;
		this.indexJ = indexJ;
	}

	public double getDoubleValue() {
		return mat[indexI][indexJ];
	}
}
