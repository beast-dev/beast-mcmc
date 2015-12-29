/*
 * CurveMouseClickListener.java
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
 * A CurveMouseClickListener defines an interface to handle mouse click
 * performed on a scatterplot and passed to a curve definition.
 * @see Scatterplot
 * @see CurveDefinition
 *
 *
 * @author Didier H. Besset
 */
public interface CurveMouseClickListener
{


/**
 * Processes the mouse clicks received from the scatterplot.
 * If a mouse listener has been defined, each point is tested for mouse click.
 * If the mouse click falls within the symbol size of a point, the index of
 * that point is passed to the mouse listener, along with the defined parameter.
 * @see #setMouseListener
 * @param index index of the curve point on which the mouse was clicked.
 * @param param the curve identifier.
 */
public boolean handleMouseClick( int index, Object param);
}