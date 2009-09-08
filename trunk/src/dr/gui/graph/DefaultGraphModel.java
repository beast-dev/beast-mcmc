/*
 * DefaultGraphModel.java
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

import java.util.ArrayList;
import java.util.Iterator;

public class DefaultGraphModel implements GraphModel {

	public void addVertex(Vertex vertex) { vertices.add(vertex); }
	public void removeVertex(Vertex vertex) { vertices.remove(vertex); }

	public int getVertexCount() { return vertices.size(); }
	public Iterator getVertices() { return vertices.iterator(); }
	public Vertex getVertex(int index) { return vertices.get(index); }
		
	public Vertex getVertexWithUserObject(Object userObject) {
		Iterator iter = getVertices();
		while (iter.hasNext()) {
			Vertex vertex = (Vertex)iter.next();
			if (vertex.getUserObject() == userObject) return vertex;
		}
		return null;
	}
	
	public void addEdge(Edge edge) { edges.add(edge); }
	public void removeEdge(Edge edge) { edges.remove(edge); }

	public int getEdgeCount() { return edges.size(); }
	public Iterator getEdges() { return edges.iterator(); }
	public Edge getEdge(int index) { return edges.get(index); }
	
	public void clear() { vertices.clear(); edges.clear(); }
	
	private ArrayList<Vertex> vertices = new ArrayList<Vertex>();
	private ArrayList<Edge> edges = new ArrayList<Edge>();
}