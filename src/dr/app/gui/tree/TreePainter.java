/*
 * TreePainter.java
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

package dr.app.gui.tree;

import dr.evolution.tree.Tree;

import java.awt.*;
import java.awt.geom.Point2D;

public interface TreePainter  {

	/**
	*	Set line style
	*/
	void setLineStyle(Stroke lineStroke, Paint linePaint);

	/**
	*	Set hilight style
	*/
	void setHilightStyle(Stroke hilightStroke, Paint hilightPaint);

	/**
	 *	Set label style.
	 */
	void setLabelStyle(Font labelFont, Paint labelPaint);

	/**
	 *	Set hilight label style.
	 */
	void setHilightLabelStyle(Font hilightLabelFont, Paint hilightLabelPaint);

	/**
	 * Do the actual painting.
	 */
	void paintTree(Graphics2D g, Dimension size, Tree tree);

	/**
	*	Find the node under point. Returns -1 if not found.
	*/
	public int findNodeAtPoint(Point2D point);

}
