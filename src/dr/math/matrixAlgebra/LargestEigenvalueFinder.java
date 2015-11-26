/*
 * LargestEigenvalueFinder.java
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

package dr.math.matrixAlgebra;

import dr.math.iterations.IterativeProcess;

/**
 * Object used to find the largest eigen value and the corresponding
 * eigen vector of a matrix by successive approximations.
 *
 * @author Didier H. Besset
 */
public class LargestEigenvalueFinder extends IterativeProcess
{
/**
 * Eigenvalue
 */
	private double eigenvalue;
/**
 * Eigenvector
 */
	private Vector eigenvector;
/**
 * Eigenvector of transposed matrix
 */
	private Vector transposedEigenvector;
/**
 * Matrix.
 */
	private Matrix matrix;

/**
 * Constructor method.
 * @param prec double
 * @param a MatrixAlgebra.Matrix
 */
public LargestEigenvalueFinder ( double prec, Matrix a)
{
	this(a);
	this.setDesiredPrecision ( prec);
}
/**
 * Constructor method.
 * @param a MatrixAlgebra.Matrix
 */
public LargestEigenvalueFinder ( Matrix a) 
{
	matrix = a;
	eigenvalue = Double.NaN;
}
/**
 * Returns the eigen value found by the receiver.
 * @return double
 */
public double eigenvalue ( )
{
	return eigenvalue;
}
/**
 * Returns the normalized eigen vector found by the receiver.
 * @return MatrixAlgebra.Vector
 */
public Vector eigenvector ( )
{
	return eigenvector.product( 1.0 / eigenvector.norm());
}
/**
 * Iterate matrix product in eigenvalue information.
 */
public double evaluateIteration()
{
	double oldEigenvalue = eigenvalue;
	transposedEigenvector = 
						transposedEigenvector.secureProduct( matrix);
	transposedEigenvector = transposedEigenvector.product( 1.0 
							/ transposedEigenvector.components[0]);
	eigenvector = matrix.secureProduct( eigenvector);
	eigenvalue = eigenvector.components[0];
	eigenvector = eigenvector.product( 1.0 / eigenvalue);
	return Double.isNaN( oldEigenvalue)
					? 10 * getDesiredPrecision()
					: Math.abs( eigenvalue - oldEigenvalue);
}
/**
 * Set result to undefined.
 */
public void initializeIterations()
{
	eigenvalue = Double.NaN;
	int n = matrix.columns();
	double [] eigenvectorComponents = new double[ n];
	for ( int i = 0; i < n; i++) { eigenvectorComponents [i] = 1.0;}
	eigenvector = new Vector( eigenvectorComponents);
	n = matrix.rows();
	eigenvectorComponents = new double[ n];
	for ( int i = 0; i < n; i++) { eigenvectorComponents [i] = 1.0;}
	transposedEigenvector = new Vector( eigenvectorComponents);
}
/**
 * Returns a finder to find the next largest eigen value of the receiver's matrix.
 * @return MatrixAlgebra.LargestEigenvalueFinder
 */
public LargestEigenvalueFinder nextLargestEigenvalueFinder ( )
{
	double norm = 1.0 / eigenvector.secureProduct(
											transposedEigenvector);
	Vector v1 = eigenvector.product( norm);
	return new LargestEigenvalueFinder( getDesiredPrecision(),
			matrix.secureProduct(SymmetricMatrix.identityMatrix(
				v1.dimension()).secureSubtract(v1.tensorProduct(
											transposedEigenvector))));
}
/**
 * Returns a string representation of the receiver.
 * @return java.lang.String
 */
public String toString()
{
	StringBuffer sb = new StringBuffer();
	sb.append( eigenvalue);
	sb.append(" (");
	sb.append( eigenvector.toString());
	sb.append(')');
	return sb.toString();
}
}