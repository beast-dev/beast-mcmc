package dr.app.tools;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.inference.model.MatrixParameter;
import dr.math.matrixAlgebra.Matrix;
import dr.util.TIFFWriter;

import java.io.DataOutputStream;
import java.io.FileOutputStream;

/**
 * @author Marc Suchard
 *         Based on the class CalculateSplitRates in JEBL
 */
class DensityMap {
	private final String SEP = "\t";
	private final String DBL = "%5.4f";

	private int binX;
	private int binY;

	private int[][] data;
	private int[] counts;
	private int count;
	private double startX;
	private double endX;
	private double startY;
	private double endY;
	private double scaleX;
	private double scaleY;

	private double minX = Double.POSITIVE_INFINITY;
	private double maxX = Double.NEGATIVE_INFINITY;
	private double minY = Double.POSITIVE_INFINITY;
	private double maxY = Double.NEGATIVE_INFINITY;

	private double upperX;
	private double lowerX;
	private double upperY;
	private double lowerY;

	private boolean jointDensity = false;
	private boolean isCalibrated = false;

	private int slice;

	public DensityMap(int slice, int binX, int binY,
	                  double upperX, double lowerX,
	                  double upperY, double lowerY) {
		this.slice = slice;
		this.binX = binX;
		this.binY = binY;
		data = new int[binX][binY];
		counts = new int[binX];
		count = 0;
		this.upperX = upperX;
		this.lowerX = lowerX;
		this.upperY = upperY;
		this.lowerY = lowerY;
	}

	public void calibrate(Tree tree, String attributeName) {
		boolean foundAttribute = false;

		if (isCalibrated) {
			throw new RuntimeException("Already calibrated");
		}

		if (jointDensity) {
			throw new RuntimeException("Already calibrated as a joint density map");
		}

		double height = tree.getNodeHeight(tree.getRoot());
		if (height > maxX) {
			maxX = height;
		}
		minX = 0.0;
		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double) tree.getNodeAttribute(node, attributeName);
				if (value != null) {
					if (value < minY)
						minY = value;
					if (value > maxY)
						maxY = value;
					foundAttribute = true;
				}
			}
		}
		if (!foundAttribute) {
			throw new RuntimeException("Can't find any attributes, " + attributeName + ", in tree " + tree.getId());
		}

	}

	public void calibrate(Tree tree, String attributeName1, String attributeName2) {
		boolean foundAttribute1 = false;
		boolean foundAttribute2 = false;

		jointDensity = true;

		if (isCalibrated) {
			throw new RuntimeException("Already calibrated");
		}

//		double height = tree.getNodeHeight(tree.getRoot());
//		if (height > maxX) {
//			maxX = height;
//		}
		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double) tree.getNodeAttribute(node, attributeName1);
				if (value != null) {
					if (value < minX)
						minX = value;
					if (value > maxX)
						maxX = value;
					foundAttribute1 = true;
				}
				value = (Double) tree.getNodeAttribute(node, attributeName2);
				if (value != null) {
					if (value < minY)
						minY = value;
					if (value > maxY)
						maxY = value;
					foundAttribute2 = true;
				}
			}
		}
		if (!foundAttribute1) {
			throw new RuntimeException("Can't find any attributes, " + attributeName1 + ", in tree " + tree.getId());
		}

		if (!foundAttribute2) {
			throw new RuntimeException("Can't find any attributes, " + attributeName2 + ", in tree " + tree.getId());
		}
//		System.err.printf("Calibrated: minY = %3.2f, maxY = %3.2f, minX = %3.2f, maxX = %3.2f\n",minY,maxY,minX,maxX);
//		System.exit(-1);
	}

	public void addTree(Tree tree, String attributeName) {

		checkCalibration();

		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double) tree.getNodeAttribute(node, attributeName);
				if (value != null) {
					addBranch(tree.getNodeHeight(node), tree.getNodeHeight(tree.getParent(node)), value);
				}
			}
		}
	}

	private void checkCalibration() {
		if (!isCalibrated) {
			startX = minX;
			if (lowerX != Double.NEGATIVE_INFINITY) {
				startX = lowerX;
			}

			endX = maxX;
			if (upperX != Double.POSITIVE_INFINITY) {
				endX = upperX;
			}

			startY = minY;
			if (lowerY != Double.NEGATIVE_INFINITY) {
				startY = lowerY;
			}

			endY = maxY;
			if (upperY != Double.POSITIVE_INFINITY) {
				endY = upperY;
			}

			scaleX = (endX - startX) / (double) (binX - 1);  // -1 necessary to ensure that maxValue falls in the last box
			scaleY = (endY - startY) / (double) (binY - 1);

			isCalibrated = true;
		}


	}

	public void addTree(Tree tree, double sampleTime, String attributeName1, String attributeName2) {

		checkCalibration();

		double[][] variance = null;
		Object[] obj = (Object[]) tree.getAttribute(MultivariateDiffusionModel.PRECISION_TREE_ATTRIBUTE);
		if (obj != null) {
			variance = new Matrix(
					MatrixParameter.parseFromSymmetricDoubleArray(obj).getParameterAsMatrix()
			).inverse().toComponents();
		}
		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				NodeRef parent = tree.getParent(node);
				double t1 = tree.getNodeHeight(node);
				double t2 = tree.getNodeHeight(parent);
				if (t1 <= sampleTime && t2 >= sampleTime) {
					Double valueX1 = (Double) tree.getNodeAttribute(node, attributeName1);
					Double valueY1 = (Double) tree.getNodeAttribute(node, attributeName2);
					Double valueX2 = (Double) tree.getNodeAttribute(parent, attributeName1);
					Double valueY2 = (Double) tree.getNodeAttribute(parent, attributeName2);
					if (valueX1 != null && valueY1 != null && valueX2 != null && valueY2 != null) {
						addPoint(sampleTime, t1, t2, valueX1, valueY1, valueX2, valueY2, variance);
					}
				}
			}
		}
	}

	public int[][] getDensityMap() {
		return data;
	}

	private void addBranch(double start, double end, double y) {
		if (start >= endX || end <= startX) {
			// branch is outside bounds...
			return;
		}
		if (y > endY || y < startY) {
			// value is outside bounds...
			return;
		}

		// clip the branch to the bounds
		if (start < startX) {
			start = startX;
		}
		if (end > endX) {
			end = endX;
		}

		// determine bin for y
		int Y = (int) ((y - startY) / scaleY);
		// determine start and end bin for x
		int START = (int) ((start - startX) / scaleX);
		int END = (int) ((end - startX) / scaleX);

//			System.out.println(start+":"+end+" -> "+START+":"+END);
		for (int i = START; i <= END; i++) {
			data[i][Y] += 1;
			counts[i] += 1;
		}
	}

	private void addPoint(double t, double startTime, double endTime, double x0, double y0, double x1, double y1, double[][] variance) {
		double t0 = t - startTime;
		double t1 = endTime - t;
		double x, y;
		if (t0 == 0) {
			x = x0;
			y = y0;
		} else if (t1 == 0) {
			x = x1;
			y = y1;
		} else {
			x = ((x0 / t0) + (x1 / t1)) / ((1.0 / t0) + (1.0 / t1));
			y = ((y0 / t0) + (y1 / t1)) / ((1.0 / t0) + (1.0 / t1));

			if (variance != null) {
				// todo add stochastic noise
			}
		}

		if (x > endX || x < startX || y > endY || y < startY) {
			// point is outside bounds...
			return;
		}

		// determine bin for x
		int X = (int) ((x - startX) / scaleX);
		// determine bin for y
		int Y = (int) ((y - startY) / scaleY);

		data[X][Y] += 1;
		count += 1;
	}

	public String toString() {
		return toString(true);
	}

	public void writeAsTIFF(String fileName) {

		double[][] matrix = normalize(255);
		try {
			DataOutputStream tiffOut = new DataOutputStream(new FileOutputStream(fileName));
			TIFFWriter.writeDoubleArray(tiffOut, matrix);
			tiffOut.close();
		} catch (Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}


	private double[][] normalize(double max) {

		double[][] matrix = new double[binX][binY];

		double maxValue = 0;
		for (int i = 0; i < binY; i++) {
			for (int j = 0; j < binX; j++) {
				if (data[j][i] > maxValue) {
					maxValue = data[j][i];
				}
			}
		}

		for (int i = 0; i < binY; i++) {
			for (int j = 0; j < binX; j++) {
				matrix[j][i] = ((double) data[j][i] / maxValue) * max;
//				double dblCount;
//				if (jointDensity) {
//					dblCount = (double) count;
//				} else {
//					dblCount = (double) counts[j];
//				}
//				if (dblCount > 0) {
//					matrix[j][i] = (double) data[j][i] / dblCount * max;
//				} else {
//					matrix[j][i] = 0.0;
//				}
			}
		}

		return matrix;

	}

	public String toString(boolean printHeaders) {

		StringBuilder sb = new StringBuilder();
		if (printHeaders) {
			sb.append(String.format("%7.5f", (double) slice));  // todo should return 3rd dimension coordinate
			for (int i = 0; i < binX; i++) {
				sb.append(SEP);
				sb.append(String.format("%7.5f", startX + scaleX * i));
			}
			sb.append("\n");
		}

		double[][] matrix = normalize(1.0);

		for (int i = 0; i < binY; i++) {
			if (printHeaders)
				sb.append(String.format("%7.5f", startY + scaleY * i));

			for (int j = 0; j < binX; j++) {

				if (j > 0 || printHeaders)
					sb.append(SEP);

//				double dblCount;
//				if (jointDensity) {
//					dblCount = (double) count;
//				} else {
//					dblCount = (double) counts[j];
//				}
//				if (dblCount > 0) {
//					sb.append(String.format(DBL,
//							(double) data[j][i] / dblCount
//					));
//				} else {
//					sb.append(String.format(DBL, 0.0));
//				}
				sb.append(String.format(DBL, matrix[j][i]));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

}
