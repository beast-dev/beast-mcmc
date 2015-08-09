/*
 * MatrixEntryColumn.java
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

package dr.inference.loggers;

/**
 * @author Marc A. Suchard
 */
public class MatrixEntryColumn extends NumberColumn {

	private int indexI;
	private int indexJ;
	private double[][] mat;
	private Double normalization;

	public MatrixEntryColumn(String name, int indexI, int indexJ, double[][] mat, Double normalization) {
		super(name + "_" + indexI + "/" + indexJ);

		this.mat = mat;

		if (indexI < 0 || indexJ < 0 || indexI >= mat.length || indexJ >= mat[0].length)
			throw new RuntimeException("Out of bounds");

		this.indexI = indexI;
		this.indexJ = indexJ;
		this.normalization = normalization;
	}

	public MatrixEntryColumn(String name, int indexI, int indexJ, double[][] mat) {
		this(name, indexI, indexJ, mat, 1.0);
	}

	public double getDoubleValue() {
		return mat[indexI][indexJ] / normalization;
	}
}
