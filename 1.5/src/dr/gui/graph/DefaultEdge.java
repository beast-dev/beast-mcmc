/*
 * DefaultEdge.java
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
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;

public class DefaultEdge implements Edge {
	public DefaultEdge(Vertex sourceVertex, int sourceIndex, 
						Vertex destinationVertex, int destinationIndex) {
		this.sourceVertex = sourceVertex;
		this.sourceIndex = sourceIndex;
		this.destinationVertex = destinationVertex;
		this.destinationIndex = destinationIndex;
	}
	
	public Vertex getSourceVertex() {
		return sourceVertex;
	}
	
	public int getSourceIndex() {
		return sourceIndex;
	}
		
	public Vertex getDestinationVertex() {
		return destinationVertex;
	}

	public int getDestinationIndex() {
		return destinationIndex;
	}
		
	public void paint(Graphics2D g2) {

		Shape s = new Line2D.Double(getSourcePoint(), getDestinationPoint());
//		Path p = new Cub();
//		p.
		g2.setPaint(Color.black);
		g2.setStroke(new BasicStroke(2.0F));
		g2.draw(s);
	}

	public Point2D getSourcePoint() { return sourceVertex.getOutgoingEdgePoint(sourceIndex); }
	public Point2D getDestinationPoint() { return destinationVertex.getIncomingEdgePoint(destinationIndex); }
	
	private Vertex sourceVertex = null;
	private int sourceIndex = 0;
	
	private Vertex destinationVertex = null;
	private int destinationIndex = 0;
}