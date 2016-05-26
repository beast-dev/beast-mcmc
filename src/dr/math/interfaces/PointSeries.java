/*
 * PointSeries.java
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

package dr.math.interfaces;

/**
 * PointSeries is an interface used by many classes of the package numericalMethods.
 *
 * A PointSeries has the responsibility of handling mathematical
 * points in 2-dimensional space.
 * It is a BRIDGE to a vector containing the points.
 *
 * @author Didier H. Besset
 */
public interface PointSeries
{

/**
 * Returns the number of points in the series.
 */
public int size();
/**
 * Returns the x coordinate of the point at the given index.
 * @param index the index of the point.
 * @return x coordinate
 */
public double xValueAt( int index);
/**
 * Returns the y coordinate of the point at the given index.
 * @param index the index of the point.
 * @return y coordinate
 */
public double yValueAt( int index);
}