package dr.app.tools;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

/**
 * @author Marc Suchard
 * Based on the class CalculateSplitRates in JEBL
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

	public DensityMap(int binX, int binY,
	                  double upperX, double lowerX,
	                  double upperY, double lowerY) {
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
				Double value = (Double)tree.getNodeAttribute(node, attributeName);
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

		double height = tree.getNodeHeight(tree.getRoot());
		if (height > maxX) {
			maxX = height;
		}
		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double)tree.getNodeAttribute(node, attributeName1);
				if (value != null) {
					if (value < minY)
						minY = value;
					if (value > maxY)
						maxY = value;
					foundAttribute1 = true;
				}
				value = (Double)tree.getNodeAttribute(node, attributeName2);
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
	}

	public void addTree(Tree tree, String attributeName) {
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

			scaleX = (endX - startX) / (double) (binX);
			scaleY = (endY - startY) / (double) (binY);

			isCalibrated = true;
		}

		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double)tree.getNodeAttribute(node, attributeName);
				if (value != null) {
					addBranch(tree.getNodeHeight(node), tree.getNodeHeight(tree.getParent(node)), value);
				}
			}
		}
	}

	public void addTree(Tree tree, double sampleTime, String attributeName1, String attributeName2) {
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

			scaleX = (endX - startX) / (double) (binX);
			scaleY = (endY - startY) / (double) (binY);

			isCalibrated = true;
		}

		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				NodeRef parent = tree.getParent(node);
				double t1 = tree.getNodeHeight(node);
				double t2 = tree.getNodeHeight(parent);
				if (t1 <= sampleTime && t2 >= sampleTime) {
					Double valueX1 = (Double)tree.getNodeAttribute(node, attributeName1);
					Double valueY1 = (Double)tree.getNodeAttribute(node, attributeName2);
					Double valueX2 = (Double)tree.getNodeAttribute(parent, attributeName1);
					Double valueY2 = (Double)tree.getNodeAttribute(parent, attributeName2);
					if (valueX1 != null && valueY1 != null && valueX2 != null && valueY2 != null) {
						addPoint(sampleTime, t1, t2, valueX1, valueY1, valueX2, valueY2);
					}
				}
			}
		}
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

	private void addPoint(double t, double startTime, double endTime, double x0, double y0, double x1, double y1) {
		double t0 = t - startTime;
		double t1 = endTime - t;
		double x = ((x0/t0) + (x1/t1)) / ((1.0/t0) + (1.0/t1));
		double y = ((y0/t0) + (y1/t1)) / ((1.0/t0) + (1.0/t1));

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
		StringBuilder sb = new StringBuilder();
		sb.append("0.0");
		for (int i = 0; i < binX; i++) {
			sb.append(SEP);
			sb.append(String.format("%7.5f", startX + scaleX * i));
		}
		sb.append("\n");
		for (int i = 0; i < binY; i++) {
			sb.append(String.format("%7.5f", startY + scaleY * i));

			for (int j = 0; j < binX; j++) {
				sb.append(SEP);
				double dblCount;
				if (jointDensity) {
					dblCount = (double) count;
				} else {
					dblCount = (double) counts[j];
				}
				if (dblCount > 0) {
					sb.append(String.format(DBL,
							(double) data[j][i] / dblCount
					));
				} else {
					sb.append(String.format(DBL, 0.0));
				}
			}
			sb.append("\n");
		}
		return sb.toString();
	}

}
