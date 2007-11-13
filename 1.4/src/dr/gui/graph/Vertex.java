/*
 * Vertex.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.gui.graph;

import java.awt.*;
import java.awt.geom.Dimension2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

public interface Vertex {
	
	Object getUserObject();

	void setColumn(int column);
	int getColumn();
	
	void setRow(int row);
	int getRow();
	
	int getHeightInRows();
	
	int getOutgoingEdgeCount();
	Edge getOutgoingEdge(int index);
	int getOutgoingEdgeIndex(Edge edge);

	int getIncomingEdgeCount();
	Edge getIncomingEdge(int index);
	int getIncomingEdgeIndex(Edge edge);

	Point2D getLocation();
	void setLocation(Point2D location);
	
	Dimension2D getDimension();
	Rectangle2D getBounds();

	void paint(Graphics2D g2);

	int getOutgoingEdgePointCount();
	Point2D getOutgoingEdgePoint(int index);
	Rectangle2D getOutgoingEdgeBounds(int index);

	int getIncomingEdgePointCount();
	Point2D getIncomingEdgePoint(int index);
	Rectangle2D getIncomingEdgeBounds(int index);
}