package dr.evomodel.continuous;

import dr.math.SparseMatrixExponential;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 14, 2007
 * Time: 7:33:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class TopographicalMap {

	public TopographicalMap(double[][] map) {
		xDim = map.length;
		yDim = map[0].length;
		this.map = map;
		order = setUpIndices();
		setUpRandomWalkMatrix();
	}

	public int getXDim() {
		return xDim;
	}

	public int getYDim() {
		return yDim;
	}

	public int getOrder() {
		return order;
	}

	public int getNonZeroSize() {
		return nonZeroElements;
	}


	private int setUpIndices() {
		int index = 0;
		indexByXY = new int[xDim * yDim];
		int[] tmpXYByIndex = new int[xDim * yDim];
		for (int x = 0; x < xDim; x++) {
			for (int y = 0; y < yDim; y++) {
				final int xy = x * yDim + y;
				if (!Double.isNaN(map[x][y])) {
					tmpXYByIndex[index] = xy;
					indexByXY[xy] = index;
					index++;
				} else
					indexByXY[xy] = -1;
			}
		}
		xyByIndex = new int[index];
		xyByIndex = tmpXYByIndex;
		return index;
	}

	public double getCTMCProbability(double[] start, double[] stop, double time) {
		int startIndex = getIndex((int) start[0], (int) start[1]);
		int stopIndex = getIndex((int) stop[0], (int) stop[1]);

//		startIndex = stopIndex = 0;

		if (startIndex == -1 || stopIndex == -1)
			return 0;
//		System.err.println("Indices: "+startIndex+" -> "+stopIndex);
//		System.err.println("Time = "+time);
		double prob = matrixExp.getExponentialEntry(startIndex, stopIndex, time);
//		System.err.println("prob = "+prob);
		return prob;
	}


	public SparseMatrixExponential getMatrix() {
		return matrixExp;
	}

	public boolean isValidPoint(int x, int y) {
		return (getIndex(x, y) != -1);
	}

	public int getIndex(int x, int y) {
		if (x < 0 || x >= xDim || y < 0 || y >= yDim)
			return -1;
		return indexByXY[x * yDim + y];
	}

	public int getX(int index) {
		return xyByIndex[index] / yDim;
	}

	public int getY(int index) {
		final int xy = xyByIndex[index];
		final int dim = yDim;
		final int mod = xy % dim;
		return mod;
	}

	private void setUpRandomWalkMatrix() {

		List<GraphEdge> edgeList = new ArrayList<GraphEdge>();
		for (int index = 0; index < order; index++) {
			int x = getX(index);
			int y = getY(index);
			int dest;

			dest = getIndex(x - 1, y - 1);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x - 1, y - 1)
				));

			dest = getIndex(x - 1, y);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x - 1, y)
				));

			dest = getIndex(x - 1, y + 1);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x - 1, y + 1)
				));

			dest = getIndex(x, y - 1);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x, y - 1)
				));

			dest = getIndex(x, y + 1);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x, y + 1)
				));

			dest = getIndex(x + 1, y - 1);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x + 1, y - 1)
				));

			dest = getIndex(x + 1, y);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x + 1, y)
				));

			dest = getIndex(x + 1, y + 1);
			if (dest != -1)
				edgeList.add(new GraphEdge(
						index, dest, getWeight(x, y, x + 1, y + 1)
				));


		}


		double[] rowTotal = new double[order];
		nonZeroElements = edgeList.size() + order;  // don't forget the diagonals

		matrixExp = new SparseMatrixExponential(order, nonZeroElements);

//		int index = 0;
		for (GraphEdge edge : edgeList) {
			int start = edge.getStart();
			int stop = edge.getStop();
			double weight = edge.getWeight();
			rowTotal[start] -= weight;
//			bandMatrix.set(start, stop, weight);
			matrixExp.addEntry(start, stop, weight);
		}

		double norm = 0;
		for (int i = 0; i < order; i++) {
			if (-2 * rowTotal[i] > norm)
				norm = -2 * rowTotal[i];
//			bandMatrix.set(i, i, rowTotal[i]);
			matrixExp.addEntry(i, i, rowTotal[i]);

		}

		matrixExp.setNorm(norm);
	}


	public void doStuff() {
		System.err.println(matrixExp.getExponentialEntry(0, 1, 1000));
		System.err.println(matrixExp.getExponentialEntry(0, 0, 1000));


	}

	public static void writeEigenvectors(String file, double[][] mat) throws IOException {

		PrintWriter writer;
		if (file.endsWith("gz"))
			writer = new PrintWriter(
					new OutputStreamWriter
							(new GZIPOutputStream
									(new FileOutputStream(file))));
		else
			writer = new PrintWriter(new FileWriter(file));

		final int dim = mat.length;
		writer.println(dim);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++)
				writer.println(mat[i][j]);
		}
		writer.close();

	}

	public static void writeEigenvalues(String file, double[] vec) throws IOException {

		PrintWriter writer;
		if (file.endsWith("gz"))
			writer = new PrintWriter(
					new OutputStreamWriter
							(new GZIPOutputStream
									(new FileOutputStream(file))));
		else
			writer = new PrintWriter(new FileWriter(file));

		writer.println(vec.length);
		for (double d : vec)
			writer.println(d);
		writer.close();

	}

	public static int extractInt(String line) {
		StringTokenizer st = new StringTokenizer(line);
		st.nextToken();
		return Integer.parseInt(st.nextToken());
	}

	public static final String defaultInvalidString = "*";

	public static double[][] readGRASSAscii(String file) throws IOException {
		return readGRASSAscii(file, defaultInvalidString);
	}

	public static double[][] readGRASSAscii(String file, String invalidString) throws IOException {
		double[][] map;
		BufferedReader reader = getReader(file);
		String northLine = reader.readLine();
		String southLine = reader.readLine();
		String eastLine = reader.readLine();
		String westLine = reader.readLine();
		int numberRows = extractInt(reader.readLine());
		int numberCols = extractInt(reader.readLine());
		map = new double[numberRows][numberCols];
		for (int i = 0; i < numberRows; i++) {
			String line = reader.readLine();
			StringTokenizer st = new StringTokenizer(line);
			for (int j = 0; j < numberCols; j++) {
				String item = st.nextToken();
				if (item.compareTo(invalidString) == 0) {
					map[i][j] = Double.NaN;
				} else {
					map[i][j] = Double.parseDouble(item);
				}
			}
		}
		return map;
	}


	public static BufferedReader getReader(String file) throws IOException {
		BufferedReader reader;
		if (file.endsWith("gz"))
			reader = new BufferedReader
					(new InputStreamReader
							(new GZIPInputStream
									(new FileInputStream(file))));
		else
			reader = new BufferedReader(new FileReader(file));
		return reader;
	}

	public static double[] readEigenvalues(String file) throws IOException {

		BufferedReader reader = getReader(file);

		final int dim = Integer.parseInt(reader.readLine());
		double[] vec = new double[dim];
		for (int i = 0; i < dim; i++)
			vec[i] = Double.parseDouble(reader.readLine());
		reader.close();

		return vec;

	}

	public static double[][] readEigenvectors(String file) throws IOException {

		BufferedReader reader = getReader(file);

		final int dim = Integer.parseInt(reader.readLine());
		double[][] mat = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++)
				mat[i][j] = Double.parseDouble(reader.readLine());
		}
		reader.close();

		return mat;

	}

//	private double getProbability(int I, int J, double time) {
//		double probability = 0;
//		for (int k = 0; k < order; k++) {
//			probability += eVec[I][k] * Math.exp(time * eVal[k]) * eVec[J][k];
//		}
//		return probability;
//
//	}


	public class GraphEdge {
		private int start;
		private int stop;
		private double weight;

		public GraphEdge(int start, int stop, double weight) {
			this.start = start;
			this.stop = stop;
			this.weight = weight;
		}

		public int getStart() {
			return start;
		}

		public void setStart(int start) {
			this.start = start;
		}

		public int getStop() {
			return stop;
		}

		public void setStop(int stop) {
			this.stop = stop;
		}

		public double getWeight() {
			return weight;
		}

		public void setWeight(double weight) {
			this.weight = weight;
		}
	}

	public double getWeight(int x0, int y0, int x1, int y1) {
		return 1.0 / (100 * Math.abs(map[x0][y0] - map[x1][y1]) + 1.0);
	}


	public static void main(String[] args) {
		double[][] map = null;
		try {
			map = readGRASSAscii("combined.asc");
		} catch (IOException e) {
			e.printStackTrace();
		}
		new TopographicalMap(map).doStuff();
	}

	private double[][] map;
	private int[] indexByXY;
	private int[] xyByIndex;
	SparseMatrixExponential matrixExp;


	private int xDim;
	private int yDim;
	private int order;
	private int nonZeroElements;

	public static double[][] sillyMap = {
			{0.0, 100.0, 0.0, 100.0, 0.0},
			{0.0, 0.0, 0.0, 100.0, 0.0},
			{100.0, 0.0, 100.0, 100.0, 0.0},
			{100.0, 0.0, 100.0, 100.0, 0.0},
			{100.0, 0.0, 0.0, 0.0, 0.0}
	};

}
