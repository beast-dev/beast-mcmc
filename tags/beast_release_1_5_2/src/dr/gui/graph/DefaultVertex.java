/*
 * DefaultVertex.java
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
import java.util.ArrayList;


public class DefaultVertex implements Vertex {

	public DefaultVertex(Object userObject, int width, int height) {
		this.userObject = userObject;
		this.location = new Point2D.Double(0,0);
		this.dimension = new Dimension(width, height);
	}
	
	public Object getUserObject() { return userObject; }
	
	public void setColumn(int column) { this.column = column; }
	public int getColumn() { return column; }
	
	public void setRow(int row) { this.row = row; }
	public int getRow() { return row; }
	
	public int getHeightInRows() { return 1; }
	
	public int getOutgoingEdgeCount() { return outgoingEdges.size(); }
	public Edge getOutgoingEdge(int index) { return outgoingEdges.get(index); }
	public int getOutgoingEdgeIndex(Edge edge) { return outgoingEdges.indexOf(edge); }

	public void addOutgoingEdge(Edge edge) { outgoingEdges.add(edge); }

	public int getIncomingEdgeCount() { return incomingEdges.size(); }
	public Edge getIncomingEdge(int index) { return incomingEdges.get(index); }
	public int getIncomingEdgeIndex(Edge edge) { return incomingEdges.indexOf(edge); }

	public void addIncomingEdge(Edge edge) { incomingEdges.add(edge); }

	public Point2D getLocation() { return location; }
	public void setLocation(Point2D location) { this.location = location; }

	public Dimension2D getDimension() { return dimension; }
	public Rectangle2D getBounds() { 
		return new Rectangle2D.Double(
			location.getX(), location.getY(), 
			dimension.getWidth(), dimension.getHeight()); 
	}

	public void paint(Graphics2D g2) {

		g2.setPaint(Color.black);
		g2.setStroke(new BasicStroke(2.0F));
		g2.draw(getBounds());
	}

	public int getOutgoingEdgePointCount() {
		return getOutgoingEdgeCount();
	}

	public Point2D getOutgoingEdgePoint(int index) {
		
		double d = dimension.getHeight() / getOutgoingEdgeCount();		
		return new Point2D.Double(location.getX() + dimension.getWidth(), 
									location.getY() + (d * (0.5 + index)));
	}
	
	public Rectangle2D getOutgoingEdgeBounds(int index) {
		Point2D point = getOutgoingEdgePoint(index);		
		return new Rectangle2D.Double(point.getX() -4, point.getY() -4, 8, 8);
	}

	public int getIncomingEdgePointCount() {
		return getIncomingEdgeCount();
	}
	
	public Point2D getIncomingEdgePoint(int index) {
		
		double d = dimension.getHeight() / getOutgoingEdgeCount();		
		return new Point2D.Double(location.getX(), 
									location.getY() + (d * (0.5 + index)));
	}
	
	public Rectangle2D getIncomingEdgeBounds(int index) {
		Point2D point = getIncomingEdgePoint(index);		
		return new Rectangle2D.Double(point.getX() -4, point.getY() -4, 8, 8);
	}

	private int column = 0;
	private int row = 0;
	
	private final ArrayList<Edge> outgoingEdges = new ArrayList<Edge>();
	private final ArrayList<Edge> incomingEdges = new ArrayList<Edge>();
	
	private Point2D location;
	private final Dimension2D dimension;
	
	private Object userObject = null;
}