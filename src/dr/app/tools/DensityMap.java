/*
 * DensityMap.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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
	private double[] average;
	private double[] testStatistic;
	private int[] counts;
	private int[] singleTreeCounts;
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
	                  double upperY, double lowerY, boolean logScale) {
		setUp(slice, binX, binY, upperX, lowerX, upperY, lowerY, logScale);
	}

	private void setUp(int slice, int binX, int binY,
	                   double upperX, double lowerX,
	                   double upperY, double lowerY, boolean logScale) {
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
		this.logScale = logScale;
	}

	public DensityMap(int binX, int binY, double[] ptX, double[] ptY) {

//		double maxX = Double.MIN_VALUE;
//		double maxY = Double.MIN_VALUE;
//		double minX = Double.MAX_VALUE;
//		double minY = Double.MAX_VALUE;

		if (ptX.length != ptY.length) {
			throw new RuntimeException("Attempting to construct an unbalanced DensityMap");
		}

		int len = ptX.length;

		for (int i = 0; i < len; i++) {
			if (ptX[i] < minX) minX = ptX[i];
			if (ptY[i] < minY) minY = ptY[i];
			if (ptX[i] > maxX) maxX = ptX[i];
			if (ptY[i] > maxY) maxY = ptY[i];
		}

		double scaleX = (maxX - minX) / (double) (binX - 1);  // -1 necessary to ensure that maxValue falls in the last box
		double scaleY = (maxY - minY) / (double) (binY - 1);

		minX -= 1.5 * scaleX; // Ensures that all boarder cells have zero mass
		maxX += 0.5 * scaleX;
		minY -= 1.5 * scaleY;
		maxY += 0.5 * scaleY;


		setUp(0, binX, binY, maxX, minX, maxY, minY, false);
//		this.minX =
		checkCalibration();

		for (int i = 0; i < len; i++) {
			addPoint(ptX[i], ptY[i]);
//			System.err.println("Added: "+ptX[i]+":"+ptY[i]);
		}

	}

	public double[] getXMidPoints() {
		if (!isCalibrated)
			throw new RuntimeException("Density map is not calibrated");
		double[] pts = new double[binX];
		pts[0] = startX + 0.5 * scaleX;
		for (int i = 1; i < binX; i++)
			pts[i] = pts[i - 1] + scaleX;
		return pts;
	}

	public double[] getYMidPoints() {
		if (!isCalibrated)
			throw new RuntimeException("Density map is not calibrated");
		double[] pts = new double[binY];
		pts[0] = startY + 0.5 * scaleY;
		for (int i = 1; i < binY; i++)
			pts[i] = pts[i - 1] + scaleY;
		return pts;
	}

	private boolean logScale = false;

	public void setLogScale(boolean logscale) {
		this.logScale = logscale;
	}

	private double transform(double d) {
		if (logScale)
			return Math.log(d);
		return d;
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
					value = transform(value);
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
//		System.out.println("Calibrated to:");
//		System.out.println("\tminX = "+minX);
//		System.out.println("\tmaxX = "+maxX);
//		System.out.println("\tminY = "+minY);
//		System.out.println("\tmaxY = "+maxY);

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
					value = transform(value);
					if (value < minX)
						minX = value;
					if (value > maxX)
						maxX = value;
					foundAttribute1 = true;
				}
				value = (Double) tree.getNodeAttribute(node, attributeName2);
				if (value != null) {
					value = transform(value);
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
					value = transform(value);
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
					Double valueX1 = transform((Double) tree.getNodeAttribute(node, attributeName1));
					Double valueY1 = transform((Double) tree.getNodeAttribute(node, attributeName2));
					Double valueX2 = transform((Double) tree.getNodeAttribute(parent, attributeName1));
					Double valueY2 = transform((Double) tree.getNodeAttribute(parent, attributeName2));
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
//			average[i] += y;
		}
	}

	private void addPoint(double x, double y) {
		if (x > endX || x < startX || y > endY || y < startY)
			return;
		int X = (int) ((x - startX) / scaleX);
		int Y = (int) ((y - startY) / scaleY);

		data[X][Y] += 1;
		count += 1;
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


	public double[][] getNormalizedDensity(double max) {
		if (!isCalibrated)
			throw new RuntimeException("Density map is not yet calibrated");
		return normalize(max);
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
