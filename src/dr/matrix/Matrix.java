/*
 * Matrix.java
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
 * An inteface for an immutable matrix.
 *
 * @version $Id: Matrix.java,v 1.8 2006/06/18 16:20:58 alexei Exp $
 *
 * @author Andrew Rambaut
 */
public interface Matrix {

	//***************************************************
	// Getter methods
	//***************************************************

	/**
	 * @return the number of rows
	 */
	int getRowCount();
	
	/**
	 * @return the number of columns
	 */
	int getColumnCount();
	
	/**
	 * @return the number of elements
	 */
	int getElementCount();
	
	/**
	 * @return the number of elements in 1 triangle. Square matrices only.
	 */
	int getTriangleCount() throws Matrix.NotSquareException;
	
	/**
	 * @return the number of elements in the diagonal. Square matrices only.
	 */
	int getDiagonalCount() throws Matrix.NotSquareException;
	
	/**
	 * @return the matrix elements as a 1D array
	 */
	double[] getElements();
	
	/**
	 * @return the matrix elements as a 2D array
	 */
	double[][] getElements2D();
	
	/**
	 * @return the upper triangle elements as a 1D array
	 */
	double[] getUpperTriangle() throws Matrix.NotSquareException;
	
	/**
	 * @return the upper triangle elements as a 1D array
	 */
	double[] getLowerTriangle() throws Matrix.NotSquareException;
	
	/**
	 * @return the diagonal elements as a 1D array
	 */
	double[] getDiagonal() throws Matrix.NotSquareException;
	
	/**
	 * @return an element
	 */
	double getElement(int row, int column);
	
	/**
	 * @return the ith element
	 */
	double getElement(int index);
	
	/**
	 * @return a single row
	 */
	double[] getRow(int row);
	
	/**
	 * @return a single column
	 */
	double[] getColumn(int column);
	
	/**
	 * @return the minimum value in the matrix
	 */
	double getMinValue();

	/**
	 * @return the minimum value in the matrix
	 */
	double getMaxValue();

	/**
	 * @return whether matrix is square
	 */
	boolean getIsSquare();
	
	/**
	 * @return whether matrix is symmetric
	 */
	boolean getIsSymmetric() throws Matrix.NotSquareException;
	
	/**
	 * @return an id for row
	 */
	String getRowId(int row);
	
	/**
	 * @return an id for row
	 */
	String getColumnId(int column);
	

	//***************************************************
	// Exception classes
	//***************************************************
	
	/**
	 * The base matrix exception
	 */
	public class MatrixException extends Exception { 
		/**
		 * 
		 */
		private static final long serialVersionUID = -5904166681730282246L;
		MatrixException() { super(); }
		MatrixException(String  message) { super(message); }
	}
	
	/**
	 * Thrown when the matrix should be square but it's not
	 */
	public class NotSquareException extends MatrixException {

		/**
		 * 
		 */
		private static final long serialVersionUID = 5121968928197320497L; }
	
	/**
	 * Thrown when the dimensions of a matrix are wrong
	 */
	public class WrongDimensionException extends MatrixException { 
		/**
		 * 
		 */
		private static final long serialVersionUID = -1799942797975356399L;
		WrongDimensionException() { super(); }
		WrongDimensionException(String message) { super(message); }
	}

	//***************************************************
	// Util class
	//***************************************************
	
	/**
	 * Static class providing maths and utility functions
	 */
	public static class Util {

		/**
		 * @return the dot product of a & b
		 */
		public static double dotProduct(Matrix a, Matrix b) {
			throw new RuntimeException("not implemented yet");
		}
		
		/**
		 * the cross product of a & b storing in result
		 */
		public static void product(Matrix a, Matrix b, MutableMatrix result) throws Matrix.WrongDimensionException {
			int rca = a.getRowCount();
			int cca = a.getColumnCount();
			
			int rcb = b.getRowCount();
			int ccb = b.getColumnCount();
			
			if (cca != rcb)
				throw new Matrix.WrongDimensionException("column count of matrix a = " + cca + ", row count of matrix b = " + rcb);
			
			result.setDimension(rca, ccb);
			
			for (int r = 0; r < rca; r++) {
				for (int c = 0; c < ccb; c++) {
					double sum = 0;
					for (int i = 0; i < cca; i++) {
						sum += a.getElement(r, i) * b.getElement(i, c);
					}
					result.setElement(r, c, sum);
				}
			}
		}
		
		/**
		 * the kronecker product of a & b storing in result
		 */
		public static void kroneckerProduct(Matrix a, Matrix b, MutableMatrix result) {
			int ra, ca, rb, cb;
			int rca = a.getRowCount();
			int cca = a.getColumnCount();
			int rcb = b.getRowCount();
			int ccb = b.getColumnCount();
			
			for (rb = 0; rb < rcb; rb++) {
				for (cb = 0; cb < ccb; cb++) {
					for (ra = 0; ra < rca; ra++) {
						for (ca = 0; ca < cca; ca++) {
							result.setElement(rb*rca + ra, cb * cca + ca, 
								a.getElement(ra, ca) * b.getElement(rb, cb));
						}
					}
							
				}
			}
		}
		
		/**
		 * the adds a to b storing in result 
		 */
		public static void add(Matrix a, Matrix b, MutableMatrix result) throws Matrix.WrongDimensionException {
			int rca = a.getRowCount();
			int cca = a.getColumnCount();
			
			if (rca != b.getRowCount() || cca != b.getColumnCount())
				throw new Matrix.WrongDimensionException();
			
			result.setDimension(rca, cca);
			
			for (int r = 0; r < rca; r++) {
				for (int c = 0; c < cca; c++) {
					result.setElement(r, c, a.getElement(r, c) + b.getElement(r, c));
				}
			}
		}
		
		/**
		 * the subtracts b from a storing in result
		 */
		public static void subtract(Matrix a, Matrix b, MutableMatrix result) throws Matrix.WrongDimensionException {
			int rca = a.getRowCount();
			int cca = a.getColumnCount();
			
			if (rca != b.getRowCount() || cca != b.getColumnCount())
				throw new Matrix.WrongDimensionException();
			
			result.setDimension(rca, cca);
			
			for (int r = 0; r < rca; r++) {
				for (int c = 0; c < cca; c++) {
					result.setElement(r, c, a.getElement(r, c) - b.getElement(r, c));
				}
			}
		}
		
		/**
		 * @return the determinant of matrix
		 */
		public static double det(Matrix matrix) throws Matrix.NotSquareException {
			
			int n = matrix.getRowCount();
			int col = matrix.getColumnCount();
			if (n != col) throw new Matrix.NotSquareException();
			
			double[][] D = new double[n][n];
			
			for (int i =0; i < n; i++) {
				for (int j=0; j < n; j++) {
					D[i][j] = matrix.getElement(i, j);
				}
			}
			
			org.apache.commons.math.linear.RealMatrixImpl RM = new org.apache.commons.math.linear.RealMatrixImpl(D);
			return RM.getDeterminant(); 
		}
		
		/**
		 * @return the log of the determinant of matrix
		 */
		public static double logDet(Matrix matrix) throws Matrix.NotSquareException {
			throw new RuntimeException("not implemented yet");
		}
		
		/**
		 * invert the matrix
		 * code modified from: http://www.nauticom.net/www/jdtaft/JavaMatrix.htm
		 * 
		 * This is the Java source code for a matrix inversion routine which uses 
		 * Gaussian reduction. The matrix to be inverted is assumed square and 
		 * nonsingular. The matrix order is n. NOTE: row 0 and column 0 of the 
		 * array D are NOT used. Matrix D is  sized as D[n+1][2*n+2], with the 
		 * left half (excluding row 0 and column 0) initialized to the matrix to 
		 * be inverted.  Upon completion the left half is reduced to an identity 
		 * matrix and the right half is reduced to the inverse of the original matrix. 
		 * The inversion routine initializes the right half of D to an identity matrix.
		 * @fixme Shouldn't this throw an exception if matrix is singular? 
		 */
		public static void invert(MutableMatrix matrix) throws Matrix.NotSquareException {
			
			int n = matrix.getRowCount();
			int col = matrix.getColumnCount();
			if (n != col) throw new Matrix.NotSquareException();
			
			double[][] D = new double[n+1][n*2+2];
			
			for (int i =0; i < n; i++) {
				for (int j=0; j < n; j++) {
					D[i+1][j+1] = matrix.getElement(i, j);
				}
			}
			
			double alpha;
			double beta;
			int i;
			int j;
			int k;
//			int error;

//			error = 0;
			int n2 = 2*n;

			/* init the reduction matrix  */
			for( i = 1; i <= n; i++ )
			{
				for( j = 1; j <= n; j++ )
				{
					D[i][j+n] = 0.;
				}
				D[i][i+n] = 1.0;
			}

			/* perform the reductions  */
			for( i = 1; i <= n; i++ )
			{
				alpha = D[i][i];
				if( alpha == 0.0 ) /* error - singular matrix */
				{
//					error = 1;
					break;
				}
				else
				{
					for( j = 1; j <= n2; j++ )
					{
						D[i][j] = D[i][j]/alpha;
					}
					for( k = 1; k <= n; k++ )
					{
						if( (k-i) != 0 )
						{
							beta = D[k][i];
							for( j = 1; j <= n2; j++ )
							{
								D[k][j] = D[k][j] - beta*D[i][j];
							}
						}
					}
				}
			}
			
			for (i =0; i < n; i++) {
				for (j=0; j < n; j++) {
					matrix.setElement(i, j, D[i+1][j+n+1]);
				}
			}
		}
		
		/**
		 * invert the product of a and b storing int result
		 */
		public static void invert(Matrix a, Matrix b, MutableMatrix result) 
			throws Matrix.NotSquareException, Matrix.WrongDimensionException
		{
			throw new RuntimeException("not implemented yet");
		}
		
		/**
		 * raise the matrix to scalar d
		 */
		public static void raise(Matrix matrix, double d, MutableMatrix result) throws Matrix.NotSquareException {
			throw new RuntimeException("not implemented yet");
		}
		
		public static Matrix createColumnVector(double[] v) { return new ColumnVector(v); }
		
		public static Matrix createRowVector(double[] v) { return new RowVector(v); }
		
		public static MutableMatrix createMutableMatrix(double[][] values) { return new ConcreteMatrix(values); }
		
	}
	
	//***************************************************
	// AbstractMatrix class
	//***************************************************
	
	/**
	 * Abstract class providing a base for immutable matrix implementations. At the
	 * very least, classes derived from this need only implement 
	 * getRowCount(), getColumnCount(), and getElement(int, int).
	 *
	 * By using these access methods, it will not be very efficient.
	 */
	public abstract class AbstractMatrix implements Matrix {
		
		//***************************************************
		// Getter methods
		//***************************************************
				
		/**
		 * @return the number of elements
		 */
		public int getElementCount() {
			return getRowCount() * getColumnCount();
		}
		
		/**
		 * @return the number of elements in 1 triangle. Square matrices only.
		 */
		public int getTriangleCount() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			int dim = getRowCount();
			return ((dim - 1) * dim) / 2;
		}
		
		/**
		 * @return the number of elements in the diagonal. Square matrices only.
		 */
		public int getDiagonalCount() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			return getRowCount();
		}
		
		/**
		 * @return the matrix elements as a 1D array
		 */
		public double[] getElements() {
			double[] values = new double[getElementCount()];
			int k = 0;
			int rc = getRowCount();
			int cc = getColumnCount();
			for (int r = 0; r < rc; r++) {
				for (int c = 0; c < cc; c++) {
					values[k] = getElement(r, c);
					k++;
				}
			}
			return values;
		}
		
		/**
		 * @return the matrix elements as a 2D array
		 */
		public double[][] getElements2D() {
		
			double[][] values = new double[getRowCount()][getColumnCount()];
			int rc = getRowCount();
			int cc = getColumnCount();
			for (int r = 0; r < rc; r++) {
				for (int c = 0; c < cc; c++) {
					values[r][c] = getElement(r, c);
				}
			}
			
			return values;
		}
		
		/**
		 * @return the ith element
		 */
		public double getElement(int index) {
			int r = index / getColumnCount();
			int c = index % getColumnCount();
			return getElement(r, c);
		}

		/**
		 * @return the upper triangle elements as a 1D array
		 */
		public double[] getUpperTriangle() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
			
			double[] values = new double[getTriangleCount()];
			int k = 0;
			int dim = getRowCount();
			for (int r = 0; r < dim; r++) {
				for (int c = r + 1; c < dim; c++) {
					values[k] = getElement(r, c);
					k++;
				}
			}
			
			return values;
		}
		
		/**
		 * @return the upper triangle elements as a 1D array
		 */
		public double[] getLowerTriangle() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
			
			double[] values = new double[getTriangleCount()];
			int k = 0;
			int dim = getRowCount();
			for (int r = 0; r < dim; r++) {
				for (int c = 0; c < r; c++) {
					values[k] = getElement(r, c);
					k++;
				}
			}

			return values;
		}
		
		
		/**
		 * @return the diagonal elements as a 1D array
		 */
		public double[] getDiagonal() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			int dim = getRowCount();
			double[] values = new double[dim];
			for (int r = 0; r < dim; r++) {
				values[r] = getElement(r, r);
			}
			
			return values;
		}
		
		/**
		 * @return a single row
		 */
		public double[] getRow(int row) {
			int dim = getColumnCount();
			double[] values = new double[dim];
			for (int c = 0; c < dim; c++) {
				values[c] = getElement(row, c);
			}
			
			return values;
		}
		
		/**
		 * @return a single column
		 */
		public double[] getColumn(int column) {
			int dim = getRowCount();
			double[] values = new double[dim];
			for (int r = 0; r < dim; r++) {
				values[r] = getElement(r, column);
			}
			
			return values;
		}

		/**
		 * @return the minimum value in the matrix
		 */
		public double getMinValue() {
			double value, minValue = getElement(0);
			int n = getElementCount();
			
			for (int i = 1; i < n; i++) {
				value = getElement(i);
				if (value < minValue)
					minValue = value;
			}
			
			return minValue;
		}

		/**
		 * @return the minimum value in the matrix
		 */
		public double getMaxValue() {
			double value, maxValue = getElement(0);
			int n = getElementCount();
			
			for (int i = 1; i < n; i++) {
				value = getElement(i);
				if (value > maxValue)
					maxValue = value;
			}
			
			return maxValue;
		}
		
		/**
		 * @return whether matrix is square
		 */
		public boolean getIsSquare() {
			return getRowCount() == getColumnCount();
		}
		
		/**
		 * @return whether matrix is symmetric
		 */
		public boolean getIsSymmetric() throws Matrix.NotSquareException {
			if (!getIsSquare())
				throw new Matrix.NotSquareException();
				
			int dim = getRowCount();
			for (int r = 0; r < dim; r++) {
				for (int c = r + 1; r < dim; r++) {
					if (getElement(r, c) != getElement(c, r))
						return false;
				}
			}
			
			return true;
		}			
		
		/**
		 * @return an id for row
		 */
		public String getRowId(int row) { return null; }
		
		/**
		 * @return an id for row
		 */
		public String getColumnId(int column) { return null; }
		

	}

}

class ColumnVector extends Matrix.AbstractMatrix {
	
	public ColumnVector(double[] v) {
		
		this.values = new double[v.length];
		for (int i = 0; i < values.length; i++) {
			values[i] = v[i];	
		}
	}
	
	public final int getColumnCount() { return 1; } 
	public final int getRowCount() { return values.length; }
	public final double getElement(int i, int j) { return values[i]; }

	double[] values = null;
}

class RowVector extends Matrix.AbstractMatrix {
	
	public RowVector(double[] v) {
		
		this.values = new double[v.length];
		for (int i = 0; i < values.length; i++) {
			values[i] = v[i];	
		}
	}
	
	public final int getRowCount() { return 1; } 
	public final int getColumnCount() { return values.length; }
	public final double getElement(int i, int j) { return values[j]; }
	
	double[] values = null;
}