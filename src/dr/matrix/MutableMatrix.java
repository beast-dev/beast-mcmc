/*
 * MutableMatrix.java
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

package dr.matrix;

/**
 * An inteface for a mutable matrix.
 *
 * @version $Id: MutableMatrix.java,v 1.3 2005/05/24 20:26:01 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public interface MutableMatrix extends Matrix {

	//***************************************************
	// Setter methods
	//***************************************************
	
	/**
	 * sets the dimension of a square matrix
	 */
	void setDimension(int rows);

	/**
	 * sets the dimension of a rectangular matrix
	 */
	void setDimension(int rows, int columns);
	
	/**
	 * sets the matrix from another matrix object
	 */
	void setMatrix(Matrix matrix);
	
	/**
	 * sets the matrix from an 2D array of doubles
	 */
	void setElements(double[][] values);
	
	/**
	 * sets the matrix from an 1D array of doubles
	 */
	void setElements(double[] values);
	
	/**
	 * sets all elements in a matrix with a single value
	 */
	void setElements(double value);
	
	/**
	 * sets a single element
	 */
	void setElement(int row, int column, double value);
	
	/**
	 * sets the ith element
	 */
	void setElement(int index, double value);
	
	/**
	 * sets a single row
	 */
	void setRow(int row, double[] values);
	
	/**
	 * sets a single column
	 */
	void setColumn(int column, double[] values);
	
	/**
	 * sets the upper triangle
	 */
	void setUpperTriangle(double[] values) throws Matrix.NotSquareException;
	
	/**
	 * sets the lower triangle
	 */
	void setLowerTriangle(double[] values) throws Matrix.NotSquareException;
	
	/**
	 * sets the diagonal
	 */
	void setDiagonal(double[] values) throws Matrix.NotSquareException;
	
	//***************************************************
	// Manipulation methods
	//***************************************************
	
	/**
	 * transposes this matrix
	 */
	void makeTransposed();
	
	/**
	 * makes the matrix symmetrical from the upper triangle
	 */
	void makeSymmetricFromUpperTriangle() throws Matrix.NotSquareException;
	
	/**
	 * makes the matrix symmetrical from the lower triangle
	 */
	void makeSymmetricFromLowerTriangle() throws Matrix.NotSquareException;
	
	/**
	 * makes the this matrix into an identity matrix
	 */
	void makeIdentity() throws Matrix.NotSquareException;
	
	//***************************************************
	// AbstractMutableMatrix class
	//***************************************************
	
	/**
	 * Abstract class providing a base for mutable matrix implementations. At the
	 * very least, classes derived from this need only implement 
	 * setDimension(int, int), getRowCount(), getColumnCount(), 
	 * setElement(int, int) and getElement(int, int).
	 *
	 * By using these access methods, it will not be very efficient.
	 */
	public abstract class AbstractMutableMatrix extends Matrix.AbstractMatrix 
												implements MutableMatrix {
		
		/**
		 * Empty constructor
		 */
		public AbstractMutableMatrix() {
		}
		
		/**
		 * Construct a square matrix
		 */
		public AbstractMutableMatrix(int rows) {
			setDimension(rows);
		}
		
		/**
		 * Construct a rectangular matrix
		 */
		public AbstractMutableMatrix(int rows, int columns) {
			setDimension(rows, columns);
		}
		
		/**
		 * Constructor which copies from matrix
		 */
		public AbstractMutableMatrix(Matrix matrix) {
			setDimension(matrix.getRowCount(), matrix.getColumnCount());
			setMatrix(matrix);
		}
		
		/**
		 * Constructor which copies from 2D array of doubles
		 */
		public AbstractMutableMatrix(double[][] values) {
			setDimension(values.length, values[0].length);
			setElements(values);
		}
		
		/**
		 * Constructor which copies from 1D array of doubles. 
		 */
		public AbstractMutableMatrix(int rows, int columns, double[] values) {
			setDimension(rows, columns);
			setElements(values);
		}
				
		/**
		 * Constructor which initialises every element to value. 
		 */
		public AbstractMutableMatrix(int rows, int columns, double value) {
			setDimension(rows, columns);
			setElements(value);
		}
				
		//***************************************************
		// Setter methods
		//***************************************************
		
		/**
		 * sets the dimension of a square matrix
		 */
		public void setDimension(int rows) {
			setDimension(rows, rows);
		}

		/**
		 * sets the matrix from another matrix object. Both matrices
		 * must be the same size.
		 */
		public void setMatrix(Matrix matrix) {
			int rc = matrix.getRowCount();
			int cc = matrix.getColumnCount();
			for (int r = 0; r < rc; r++) {
				for (int c = 0; c < cc; c++) {
					setElement(r, c, matrix.getElement(r, c));
				}
			}
		}
		
		/**
		 * sets the matrix from an 2D array of doubles. Both matrices
		 * must be the same size.
		 */
		public void setElements(double[][] values) {
			setDimension(values.length, values[0].length);
			for (int r = 0; r < values.length; r++) {
				for (int c = 0; c < values[0].length; c++) {
					setElement(r, c, values[r][c]);
				}
			}
		}
		
		/**
		 * sets the matrix from an 1D array of doubles. The dimension of
		 * the matrix must have already be set appropriately.
		 */
		public void setElements(double[] values) {
			int k = 0;
			int rc = getRowCount();
			int cc = getColumnCount();
			for (int r = 0; r < rc; r++) {
				for (int c = 0; c < cc; c++) {
					setElement(r, c, values[k]);
					k++;
				}
			}
		}
		
		/**
		 * sets all elements in a matrix with a single value
		 */
		public void setElements(double value) {
			int rc = getRowCount();
			int cc = getColumnCount();
			for (int r = 0; r < rc; r++) {
				for (int c = 0; c < cc; c++) {
					setElement(r, c, value);
				}
			}
		}
		
		/**
		 * sets the ith element
		 */
		public void setElement(int index, double value) {
			int r = index / getColumnCount();
			int c = index % getColumnCount();
			setElement(r, c, value);
		}
		
		/**
		 * sets a single row
		 */
		public void setRow(int row, double[] values) {
			int cc = getColumnCount();
			for (int c = 0; c < cc; c++) {
				setElement(row, c, values[c]);
			}
		}
		
		
		/**
		 * sets a single column
		 */
		public void setColumn(int column, double[] values) {
			int rc = getRowCount();
			for (int r = 0; r < rc; r++) {
				setElement(r, column, values[r]);
			}
		}
		
		/**
		 * sets the upper triangle
		 */
		public void setUpperTriangle(double[] values) throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
			
			int k = 0;
			int dim = getRowCount();
			for (int r = 0; r < dim; r++) {
				for (int c = r + 1; c < dim; c++) {
					setElement(r, c, values[k]);
					k++;
				}
			}
		}
		
		/**
		 * sets the lower triangle
		 */
		public void setLowerTriangle(double[] values) throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
			
			int k = 0;
			int dim = getRowCount();
			for (int r = 0; r < dim; r++) {
				for (int c = 0; c < r; c++) {
					setElement(r, c, values[k]);
					k++;
				}
			}
		}
		
		/**
		 * sets the diagonal
		 */
		public void setDiagonal(double[] values) throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			int dim = getRowCount();
			for (int r = 0; r < dim; r++) {
				setElement(r, r, values[r]);
			}
		}
		
		//***************************************************
		// Manipulation methods
		//***************************************************
		
		/**
		 * transposes this matrix
		 */
		public void makeTransposed() {
			double[][] values = getElements2D();
			
			int cc = getRowCount();
			int rc = getColumnCount();
			setDimension(rc, cc);
			for (int r = 0; r < rc; r++) {
				for (int c = 0; c < cc; c++) {
					setElement(r, c, values[c][r]);
				}
			}
		}
		
		/**
		 * makes the matrix symmetrical from the upper triangle
		 */
		public void makeSymmetricFromUpperTriangle() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			double[] values = getUpperTriangle();
			setLowerTriangle(values);
		}
		
		/**
		 * makes the matrix symmetrical from the lower triangle
		 */
		public void makeSymmetricFromLowerTriangle() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			double[] values = getUpperTriangle();
			setLowerTriangle(values);
		}
		
		/**
		 * makes the this matrix into an identity matrix
		 */
		public void makeIdentity() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
		
			setElements(0.0);
			for (int i = 0; i < getRowCount(); i++)
				setElement(i, i, 1.0);
		}
		
	}

}