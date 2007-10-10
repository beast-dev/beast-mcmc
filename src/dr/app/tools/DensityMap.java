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
	private double startX;
	private double startY;
	private double scaleX;
	private double scaleY;

	private double maxTreeHeight = 0;
	private double minValue = Double.POSITIVE_INFINITY;
	private double maxValue = Double.NEGATIVE_INFINITY;

	private double edgeFraction;

	private boolean isCalibrated = false;

	public DensityMap(int binX, int binY, double edgeFraction) {
		this.binX = binX;
		this.binY = binY;
		data = new int[binX][binY];
		counts = new int[binX];
		this.edgeFraction = edgeFraction;
	}

	public void calibrate(Tree tree, String attributeName) {
		if (isCalibrated) {
			throw new RuntimeException("Already calibrated");
		}

		double height = tree.getNodeHeight(tree.getRoot());
		if (height > maxTreeHeight) {
			maxTreeHeight = height;
		}
		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double)tree.getNodeAttribute(node, attributeName);
				if (value == null) {
					throw new RuntimeException("Node is missing attribute, " + attributeName);
				}
				if (value < minValue)
					minValue = value;
				if (value > maxValue)
					maxValue = value;
			}
		}
	}

	public void addTree(Tree tree, String attributeName) {
		if (!isCalibrated) {
			double spread = maxValue - minValue;
			minValue -= spread * edgeFraction;
			maxValue += spread * edgeFraction;

			startX = 0.0;
			startY = minValue;
			double endX = maxTreeHeight * (1.0 + edgeFraction);
			double endY = maxValue;

			scaleX = (endX - startX) / (double) binX;
			scaleY = (endY - startY) / (double) binY;

			isCalibrated = true;
		}

		for (int i = 0; i < tree.getNodeCount(); i++) {
			NodeRef node = tree.getNode(i);
			if (node != tree.getRoot()) {
				Double value = (Double)tree.getNodeAttribute(node, attributeName);
				addBranch(tree.getNodeHeight(node),
						tree.getNodeHeight(tree.getParent(node)), value);
			}
		}
	}

	private void addBranch(double start, double end, double y) {
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

	public String toString() {
//			double dblTotal = (double) total;
		StringBuilder sb = new StringBuilder();
		sb.append("0.0");
		for (int i = 0; i < binX; i++) {
			sb.append(SEP);
			sb.append(String.format("%7.5f", startX + scaleX * i));
		}
		sb.append("\n");
		for (int i = 0; i < binY; i++) {
			sb.append(String.format("%7.5f", startY + scaleY * i));
			//double dblCounts = (double)counts[i];
			for (int j = 0; j < binX; j++) {
				sb.append(SEP);
				double dblCounts = (double) counts[j];
				if (dblCounts > 0)
					sb.append(String.format(DBL,
							(double) data[j][i] / (double) counts[j]
					));
				else
					sb.append(String.format(DBL, 0.0));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

}
