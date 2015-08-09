/*
 * DecomposedMatrix.java
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

package dr.inference.model;

import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class DecomposedMatrix extends MatrixParameter {

	public static final String DIM = "dim";
	public static final String MATRIX_PARAMETER = "decomposedMatrix";

	public DecomposedMatrix(String name) {
		super(name);
	}

	public DecomposedMatrix(String name, int dim, Parameter parameter) {
		super(name);
		this.dim = dim;
		this.decomposition = parameter;
		matrix = new double[dim][dim];
		savedMatrix = new double[dim][dim];
		composeMatrix();
		addParameter(parameter);
	}

	public void parameterChangedEvent(Parameter parameter, int index) {
		compositionKnown = false;
//		System.err.println("called");
		fireParameterChangedEvent();
	}


	void composeMatrix() {

		// L_{ij} = index[i]+j
		// M_{ij} = \sum_k L_{ik} L_{jk}

		double[] values = decomposition.getParameterValues();

		for (int i = 0; i < dim; i++) {
			for (int j = i; j < dim; j++) {
				matrix[i][j] = 0.0;
				for (int k = 0; k <= i; k++)
					matrix[i][j] += values[index[i] + k] * values[index[j] + k];
				matrix[j][i] = matrix[i][j];
			}
		}

		compositionKnown = true;
	}

	public int getDimension() {
		return dim * dim;
	}

	protected void storeValues() {
		super.storeValues();
		for (int i = 0; i < dim; i++) {
			System.arraycopy(matrix[i], 0, savedMatrix[i], 0, dim);
		}

	}

	public double getParameterValue(int index) {
		int x = index / dim;
		int y = index - x * dim;
		return matrix[x][y];
	}

	protected void restoreValues() {
		super.restoreValues();
		for (int i = 0; i < dim; i++) {
			System.arraycopy(savedMatrix[i], 0, matrix[i], 0, dim);
		}
	}

	public double getParameterValue(int row, int col) {
//		System.err.println("row-col");
		if (!compositionKnown)
			composeMatrix();
		return matrix[row][col];

	}

	public double[][] getParameterAsMatrix() {
//		System.err.println("as-matrix");
		if (!compositionKnown)
			composeMatrix();

		return matrix;
	}


	public int getColumnDimension() {
		return dim;
	}

	public int getRowDimension() {
		return dim;
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MATRIX_PARAMETER;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Parameter parameter = (Parameter) xo.getChild(Parameter.class);
			int dim = xo.getIntegerAttribute(DIM);

			if (dim * (dim + 1) / 2 != parameter.getDimension())
				throw new XMLParseException("Dim attribute and parameter dimension do not match");

			return new DecomposedMatrix(MATRIX_PARAMETER, dim, parameter);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A diagonal matrix parameter constructed from its diagonals.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(Parameter.class),
				AttributeRule.newIntegerArrayRule(DIM, false)
		};

		public Class getReturnType() {
			return DecomposedMatrix.class;
		}
	};

	private boolean compositionKnown = false;
	private Parameter decomposition;
	private int dim;
	private double[][] matrix;
	private double[][] savedMatrix;

	private static int[] index = {0, 1, 3, 6, 10, 15, 21, 28, 36};
}
