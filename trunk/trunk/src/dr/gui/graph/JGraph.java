/*
 * JGraph.java
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

import javax.swing.*;
import java.awt.*;

public class JGraph extends JPanel {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 7654094599593328169L;

	public JGraph(GraphModel graphModel, GraphLayout graphLayout) {

		this.graphModel = graphModel;
		this.graphLayout = graphLayout;
		
		setBackground(Color.white);
	}

	/**
	 * Set the graph model.
	 */
	public void setGraphModel(GraphModel graphModel) {
		this.graphModel = graphModel;
		invalidate();
	}
	
	/**
	 * Get the graph model.
	 */
	public GraphModel getGraphModel() { return graphModel; }
	
	/**
	 * Set the graph layout.
	 */
	public void setGraphLayout(GraphLayout graphLayout) {
		this.graphLayout = graphLayout;
		invalidate();
	}
	
	/**
	 * Get the graph layout.
	 */
	public GraphLayout getGraphLayout() { return graphLayout; }
	
	/**
	 * Layout the graph on the canvas.
	 */
	public void validate() {
		if (graphModel == null || graphLayout == null) return;

		graphLayout.layoutGraph(this);
	}
	
	
	public void paintComponent(Graphics g) {
	
		if (graphModel == null) return;

		Dimension size = getSize();
		
		Graphics2D g2 = (Graphics2D)g;
		g2.setPaint(Color.white);
		g2.fillRect(0,0,size.width,size.height);
		paintGraph(g2, size);
	}
	
	private void paintGraph(Graphics2D g2, Dimension size) {

		for (int i = 0; i < graphModel.getVertexCount(); i++) {
			Vertex vertex = graphModel.getVertex(i);
			vertex.paint(g2);
		}
		
		for (int i = 0; i < graphModel.getEdgeCount(); i++) {
			Edge edge = graphModel.getEdge(i);
			edge.paint(g2);			
		}
	}
		
	/** the graph model */
	private GraphModel graphModel = null;
	
	/** the graph layout */
	private GraphLayout graphLayout = null;
}
