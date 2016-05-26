/*
 * SparseMatrixExponential.java
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

package dr.math;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 28, 2007
 * Time: 11:05:44 AM
 * To change this template use File | Settings | File Templates.
 */
public class SparseMatrixExponential {

	public static int maxKrylovBasisSize = 50;
	public static double tolerance = 1E-7;


	static {
		System.loadLibrary("ExpoKit");
	}

	private native void executeDGEXPV(int order, int nz, int maxBasis, double time,
	                                  double[] operandVector, double[] resultVector,
	                                  double tolerance, double matrixNorm,
	                                  double[] rate, int[] indexX, int[] indexY,
	                                  double[] workspace, int lengthWorkspace,
	                                  int[] intWorkspace, int lengethIntWorkspace,
	                                  int flag
	);

	private int order;
	private int nonZeroEntries;
	private int krylovBasisSize;

	private double[] workSpace;
	private int[] intWorkSpace;
	private int lengthWorkspace;
	private int lengthIntWorkspace;


	private int[] fortranIndexX;
	private int[] fortranIndexY;
	private double[] rate;

	private double[] start;

	private double[] stop;

	private double norm;

	int index = 0;

	public SparseMatrixExponential(int order, int nonZeroEntries) {
		this.order = order;
		this.nonZeroEntries = nonZeroEntries;
		setUpWorkspace();
	}

	private void setUpWorkspace() {

		if (order < maxKrylovBasisSize)
			krylovBasisSize = order - 10; // todo determine correct cut-off, consider Pade approximation for small order
		else
			krylovBasisSize = maxKrylovBasisSize - 1;

		lengthIntWorkspace = krylovBasisSize + 2;
		lengthWorkspace = order * (lengthIntWorkspace)
				+ 5 * (lengthIntWorkspace) * (lengthIntWorkspace) + 7;

		intWorkSpace = new int[lengthIntWorkspace];
		workSpace = new double[lengthWorkspace];

		start = new double[order];
		stop = new double[order];

		fortranIndexX = new int[nonZeroEntries];
		fortranIndexY = new int[nonZeroEntries];
		rate = new double[nonZeroEntries];

	}

	public void addEntry(int i, int j, double value) {
		fortranIndexX[index] = i + 1;
		fortranIndexY[index] = j + 1;
		rate[index] = value;
		index++;
	}

	public void setNorm(double norm) {
		this.norm = norm;
	}

	public void calculateInfinityNorm() {
		// todo
	}

	public double getExponentialEntry(int x, int y, double time) {

		start[x] = 1.0;
		int flag = 0;

		executeDGEXPV(order, nonZeroEntries, krylovBasisSize, time, start, stop, tolerance, norm,
				rate, fortranIndexX, fortranIndexY,
				workSpace, lengthWorkspace, intWorkSpace, lengthIntWorkspace, flag
		);

		start[x] = 0.0; // recycle
		return stop[y];  // stop gets overwritten with each call, no need to reset values
	}

	public String sparseRepresentation() {
		StringBuffer sb = new StringBuffer();
		sb.append(order + " " + nonZeroEntries + "\n");
		for (int i = 0; i < nonZeroEntries; i++) {
			sb.append(fortranIndexX[i]);
			sb.append(" ");
			sb.append(fortranIndexY[i]);
			sb.append(" ");
			sb.append(rate[i]);
			sb.append("\n");
		}
		return sb.toString();
	}

}
