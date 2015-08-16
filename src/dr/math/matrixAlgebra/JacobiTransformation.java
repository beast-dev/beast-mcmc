/*
 * JacobiTransformation.java
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
 * JacobiTransformation
 *
 * @author Didier H. Besset
 */
public class JacobiTransformation extends IterativeProcess 
{
	double[][] rows;
	double[][] transform;
	int p,q;	//Indices of the largest off-diagonal element
/**
 * Create a new instance for a given symmetric matrix.
 * @param m DhbMatrixAlgebra.SymmetricMatrix
 */
public JacobiTransformation (SymmetricMatrix m)
{
	int n = m.rows();
	rows = new double[n][n];
	for ( int i = 0; i < n; i++)
	{
		for ( int j = 0; j < n; j++)
			rows[i][j] = m.components[i][j];
	}	
}
/**
 * @return double[]
 */
public double[] eigenvalues ( )
{
	int n = rows.length;
	double[] eigenvalues = new double[n];
	for ( int i = 0; i < n; i++ )
		eigenvalues[i] = rows[i][i]; 
	return eigenvalues;
}
/**
 * @return DhbMatrixAlgebra.SymmetricMatrix
 */
public Vector[] eigenvectors ( )
{
	int n = rows.length;
	Vector[] eigenvectors = new Vector[n];
	double[] temp = new double[n];
	for ( int i = 0; i < n; i++ )
	{
		for ( int j = 0; j < n; j++)
			temp[j] = transform[j][i];
		eigenvectors[i] = new Vector( temp);
	}	
	return eigenvectors;
}
public double evaluateIteration()
{
	double offDiagonal = largestOffDiagonal();
	transform();
	return offDiagonal;
}
/**
 * @param m int
 */
private void exchange ( int m)
{
	int m1 = m + 1;
	double temp = rows[m][m];
	rows[m][m] = rows[m1][m1];
	rows[m1][m1] = temp;
	int n = rows.length;
	for ( int i = 0; i < n; i++ )
	{
		temp = transform[i][m];
		transform[i][m] = transform[i][m1];
		transform[i][m1] = temp;
	}	
}
public void finalizeIterations ( ) 
{
	int n = rows.length;
	int bound = n - 1;
	int i, m;
	while ( bound >= 0 )
	{
		m = -1;
		for ( i = 0; i < bound; i++ )
		{
			if ( Math.abs( rows[i][i]) < Math.abs( rows[i+1][i+1]) )
			{
				exchange( i);
				m = i;
			}	
		}	
		bound = m;
	}	
	return;
}
public void initializeIterations()
{
	transform = SymmetricMatrix.identityMatrix( rows.length).components;
}
/**
 * @return double	absolute value of the largest off diagonal element
 */
private double largestOffDiagonal( )
{
	double value = 0;
	double r;
	int n = rows.length;
	for (int i = 0; i < n; i++)
	{
		for ( int j = 0; j < i; j++)
		{
			r = Math.abs( rows[i][j]);
			if ( r > value )
			{
				value = r;
				p = i;
				q = j;
			}	
		}	
	}	
	return value;
}
/**
 * Returns a string representation of the system.
 * @return java.lang.String
 */
public String toString()
{
	StringBuffer sb = new StringBuffer();
	char[] separator = { '[', ' '};
	int n = rows.length;
	for ( int i = 0; i < n; i++)
	{
		separator[0] = '{';
		for ( int j = 0; j <= i; j++)
		{
			sb.append( separator);
			sb.append( rows[i][j]);
			separator[0] = ' ';
		}
	sb.append('}');
	sb.append('\n');
	}
	return sb.toString();
}
/**
 * @return DhbMatrixAlgebra.SymmetricMatrix
 */
private void transform ( )
{
	double apq = rows[p][q];
	if ( apq == 0 )
		return;
	double app = rows[p][p];
	double aqq = rows[q][q];
	double arp = ( aqq - app) * 0.5 / apq;
	double t = arp > 0 ? 1 / ( Math.sqrt( arp * arp + 1) + arp)
					   : 1 / ( arp - Math.sqrt( arp * arp + 1));
	double c = 1 / Math.sqrt( t * t + 1);
	double s = t * c;
	double tau = s / ( 1 + c);
	rows[p][p] = app - t * apq;
	rows[q][q] = aqq + t * apq;
	rows[p][q] = 0;
	rows[q][p] = 0;
	int n = rows.length;
	for ( int i = 0; i < n; i++ )
	{
		if ( i != p && i != q )
		{
			rows[p][i] = rows[i][p] - s *( rows[i][q]
												+ tau * rows[i][p]);
			rows[q][i] = rows[i][q] + s *( rows[i][p]
												- tau * rows[i][q]);
			rows[i][p] = rows[p][i];
			rows[i][q] = rows[q][i];
		}
		arp = transform[i][p];
		aqq = transform[i][q];
		transform[i][p] = arp - s * ( aqq + tau * arp);
		transform[i][q] = aqq + s * ( arp - tau * aqq);
	}	
}
}