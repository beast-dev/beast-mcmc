/*
 * HierarchicalGraphLayout.java
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
import java.util.Iterator;

public class HierarchicalGraphLayout implements GraphLayout {

	public HierarchicalGraphLayout(double hInset, double vInset, double hGap, double vGap) {
		this.hInset = hInset;
		this.vInset = vInset;
		this.hGap = hGap;
		this.vGap = vGap;
	}
		
	/**
	 * Layout the graph on the canvas.
	 */
	public void layoutGraph(JGraph graph) {

		graphModel = graph.getGraphModel();
		
		if (graphModel.getVertexCount() == 0) return;
		
		columnCount = 0;
		rowCount = 0;
		
		// reset the column and row of each vertex
		Iterator iter = graphModel.getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			vertex.setColumn(Integer.MIN_VALUE);
			vertex.setRow(Integer.MIN_VALUE);
		}
		
		// set the column of each vertex starting at each 'root' vertex
		iter = graphModel.getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			if (vertex.getIncomingEdgeCount() == 0) {
				setColumn(vertex, 0);
			}
		}
		
		// then check to see what the first and last column is
		int minColumn = 0;
		int maxColumn = 0;
		iter = graphModel.getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			int column = vertex.getColumn();
			if (column < minColumn) minColumn = column;
			if (column > maxColumn) maxColumn = column;
		}

		// adjust the column of each node to start at 0 on the left
		if (minColumn != 0) {
			iter = graphModel.getVertices();
			while (iter.hasNext()) {
				Vertex vertex = (Vertex)iter.next();
				vertex.setColumn(vertex.getColumn() - minColumn);
			}
			
			maxColumn -= minColumn;
		}
		
		columnCount = maxColumn + 1;
		
		// set the rowCounts of each column
		int[] rowCounts = new int[columnCount];	
		
		iter = graphModel.getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			int column = vertex.getColumn();
			vertex.setRow(rowCounts[column]);
//			rowCounts[column] += vertex.getHeightInRows();
			rowCounts[column] += 1;
			if (rowCounts[column] > rowCount) rowCount = rowCounts[column];
		}
		
		reorderColumns();
		
		positionComponents();
		
		graph.setPreferredSize(new Dimension((int)sumWidth, (int)sumHeight));
		graph.invalidate();
	}
	
	/**
	 * Recursive function to set the column of a vertex.
	 * If we wish to layout the graph in a vertical orientation
	 * then rows and columns are swapped after laying out.
	 */
	public int setColumn(Vertex vertex, int column) {

		if (vertex.getColumn() == Integer.MIN_VALUE) {

			vertex.setColumn(column);
						
			for (int i = 0; i < vertex.getOutgoingEdgeCount(); i++) {
				Edge edge = vertex.getOutgoingEdge(i);
				int newColumn = setColumn(edge.getDestinationVertex(), column+1);
				if (newColumn < vertex.getColumn()) vertex.setColumn(newColumn);
			}
		}
		
		return vertex.getColumn() - 1;
	}
				
	/**
	 * Heuristically search for an ordering of vertices in columns that 
	 * minimizes the number of crossing edges.
	 */
	public void reorderColumns() {
	
		boolean improved = false;
		int minCrossings = Integer.MAX_VALUE;
		double minSumLength = Double.MAX_VALUE;
		double minSumHeight = Double.MAX_VALUE;
		
		Integer[][] columnOrders = new Integer[columnCount][];
		for (int i = 0; i < columnCount; i++) {
			columnOrders[i] = new Integer[rowCount];
		}
		
		for (int i = 0; i < graphModel.getVertexCount(); i++) {
			Vertex vertex = graphModel.getVertex(i);
			columnOrders[vertex.getColumn()][vertex.getRow()] = new Integer(i);
		}
		
		Integer[] bestOrder = new Integer[rowCount];
		do {	
			improved = false;
			for (int i = columnCount - 1; i >= 0; i--) {
			
				for (int j = 0; j < bestOrder.length; j++) bestOrder[j] = columnOrders[i][j];
				
				dr.util.Permutator perm = new dr.util.Permutator(columnOrders[i]);
				while (perm.hasNext()) {
					Integer[] order = (Integer[])perm.next();
					
					setColumnOrder(i, order);
					positionComponents();
					int crossings = countCrossings();
					double sumLength = sumEdgeLengths();
					
					if (crossings < minCrossings) {
						minCrossings = crossings;
						minSumLength = sumLength;
						improved = true;
						for (int j = 0; j < bestOrder.length; j++) bestOrder[j] = order[j];
					} else if (crossings == minCrossings) {
						if (sumLength < minSumLength) {
							minSumLength = sumLength;
							improved = true;
							for (int j = 0; j < bestOrder.length; j++) bestOrder[j] = order[j];
						} else if (sumLength == minSumLength) {
							if (sumHeight < minSumHeight) {
								minSumHeight = sumHeight;
								improved = true;
								for (int j = 0; j < bestOrder.length; j++) bestOrder[j] = order[j];
							}
						}
					}
				}
					
			
				for (int j = 0; j < bestOrder.length; j++) columnOrders[i][j] = bestOrder[j];
				setColumnOrder(i, bestOrder);
			}
		} while (improved);
	}

	/**
	 * Set the row order of vertices in a column.
	 */
	public void setColumnOrder(int layer, Integer[] order) {

		for (int i = 0; i < order.length; i++) {
			if (order[i] != null) {
				int index1 = order[i].intValue();
				Vertex vertex = graphModel.getVertex(index1);
				vertex.setRow(i);
			}
		}
	}
	
	/**
	 * Count the total number of crossing edges.
	 */
	public int countCrossings() {
		
		int count = 0;
		
		for (int i = 0; i < graphModel.getEdgeCount() - 1; i++) {
			Edge edge1 = graphModel.getEdge(i);
			Point2D source1 = edge1.getSourcePoint();
			Point2D dest1 = edge1.getDestinationPoint();
			
			for (int j = i + 1; j < graphModel.getEdgeCount(); j++) {
				Edge edge2 = graphModel.getEdge(j);
				Point2D source2 = edge2.getSourcePoint();
				Point2D dest2 = edge2.getDestinationPoint();
				if (linesIntersect(source1, dest1, source2, dest2)) {
					count++;
				}
			}
		}
		return count;
	}
	
	/**
	 * Calculate the sum of all the edge lengths.
	 */
	public double sumEdgeLengths() {
		
		double sumLength = 0;
		
		for (int i = 0; i < graphModel.getEdgeCount() - 1; i++) {
			Edge edge = graphModel.getEdge(i);
			Point2D source = edge.getSourcePoint();
			Point2D dest = edge.getDestinationPoint();
			sumLength += source.distance(dest);
		}
		return sumLength;
	}
	
	/**
	 * Do two line segments defined by end points [x1,y1]->[x2,y2] & [u1,v1]->[u2,v2] 
	 * intersect? This assumes that x2 > x1 & u2 > u1.
	 */
	public boolean linesIntersect(Point2D source1, Point2D dest1, Point2D source2, Point2D dest2) {
		
		double a = (dest1.getY() - source1.getY()) / (dest1.getX() - source1.getX());
		double b = source1.getY() - a * source1.getX();
		
		double c = (dest2.getY() - source2.getY()) / (dest2.getX() - source2.getX());
		double d = source2.getY() - c * source2.getX();
		
		double x = (b - d) / (c - a);

		if ( (x > source1.getX()) && (x < dest1.getX()) && 
				(x > source2.getX()) && (x < dest2.getX()) ) return true;

		return false;						
	}
	
	private void positionComponents() {
	
		double[] columnWidths = new double[columnCount];	
		double[][] rowHeights = new double[columnCount][rowCount];	
		double minRowHeight = 10000.0;	
		
		Iterator iter = graphModel.getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			int column = vertex.getColumn();
			int row = vertex.getRow();
			
			Dimension2D d = vertex.getDimension();
			if (columnWidths[column] < d.getWidth()) columnWidths[column] = d.getWidth();

			rowHeights[column][row] = d.getHeight();
			if (minRowHeight > d.getHeight()) minRowHeight = d.getHeight();
		}

		double[] columnPositions = new double[columnCount];	
		
		sumWidth = hInset;
		for (int i = 0; i < columnCount; i++) {
			columnPositions[i] = sumWidth;
			sumWidth += columnWidths[i];
			if (i < columnCount - 1) sumWidth += hGap;
		}
		sumWidth += hInset;

		sumHeight = 0;
		double[][] rowPositions = new double[columnCount][rowCount];	
		for (int i = 0; i < columnCount; i++) {
			double h = vInset;
			for (int j = 0; j < rowCount; j++) {
				rowPositions[i][j] = h;
				if (rowHeights[i][j] == 0) {
					h += minRowHeight;
				} else {
					h += rowHeights[i][j];
				}
				if (i < rowCount - 1) h += vGap;
			}
			h += vInset;
			
			if (h > sumHeight) sumHeight = h;
		}
				
		iter = graphModel.getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			
			Dimension2D d = vertex.getDimension();
			
			double x = columnPositions[vertex.getColumn()]
							+ (columnWidths[vertex.getColumn()] - d.getWidth());
			double y = rowPositions[vertex.getColumn()][vertex.getRow()]; 
//							+ (rowHeights[vertex.getRow()] - d.getHeight()) / 2.0;
			Point2D location = new Point2D.Double(x, y);

			vertex.setLocation(location);
		}
	}
	
	private double hGap;
	private double vGap;
	
	private double hInset;
	private double vInset;
	
	private int columnCount;
	private int rowCount;
	
	private double sumWidth;
	private double sumHeight;
	
	private GraphModel graphModel;
}