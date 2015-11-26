/*
 * JTreeDisplay.java
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

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.Tree;
import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.print.PageFormat;
import java.awt.print.Printable;

public class JTreeDisplay extends JComponent implements Printable,
														MutableTreeListener,
														MutableTaxonListListener {

	/**
	 *
	 */
	private static final long serialVersionUID = 4845325190293249506L;

	/** the tree */
	protected Tree tree = null;

	/** the tree painter */
	private TreePainter treePainter = null;

	public JTreeDisplay(TreePainter treePainter) {
		this.treePainter = treePainter;
		init();
	}

	/**
	 * @param tree the tree
	 */
	public JTreeDisplay(TreePainter treePainter, Tree tree) {
		this.treePainter = treePainter;

		init();
		setTree(tree);
	}

	/**
	 * Called by all constructors.
	 */
	void init() {

//		setBorder(BorderFactory.createLineBorder(Color.black, 1));

		enableEvents(AWTEvent.MOUSE_EVENT_MASK);

		// adds a mouse listener
		addMouseListener(new MListener());
//		addMouseMotionListener(new MMListener());

	}

	/**
	 * Set the tree.
	 */
	public void setTree(Tree tree) {
		this.tree = tree;

		if (tree != null) {
			if (tree instanceof MutableTree) {
				((MutableTree)tree).addMutableTreeListener(this);
			}

			if (tree instanceof MutableTaxonList) {
				((MutableTaxonList)tree).addMutableTaxonListListener(this);
			}
		}

		repaint();
	}

	/**
	*	Set line style
	*/
	public void setLineStyle(Stroke lineStroke, Paint linePaint) {
		treePainter.setLineStyle(lineStroke, linePaint);
		repaint();
	}

	/**
	*	Set hilight style
	*/
	public void setHilightStyle(Stroke hilightStroke, Paint hilightPaint) {
		treePainter.setHilightStyle(hilightStroke, hilightPaint);
		repaint();
	}

	/**
	 *	Set label style.
	 */
	public void setLabelStyle(Font labelFont, Paint labelPaint) {
		treePainter.setLabelStyle(labelFont, labelPaint);
		repaint();
	}

	/**
	 *	Set hilight label style.
	 */
	public void setHilightLabelStyle(Font hilightLabelFont, Paint hilightLabelPaint) {
		treePainter.setHilightLabelStyle(hilightLabelFont, hilightLabelPaint);
		repaint();
	}

	public void paintComponent(Graphics g) {

		if (tree == null) return;

		Dimension size = getSize();
		Graphics2D g2d = (Graphics2D)g;
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		treePainter.paintTree(g2d, size, tree);

	}

    //********************************************************************
	// Printable interface
	//********************************************************************

	public int print(Graphics g, PageFormat pageFormat, int pageIndex) {
		if (pageIndex > 0) {
			return(NO_SUCH_PAGE);
		} else {
			Graphics2D g2d = (Graphics2D)g;

			double x0 = pageFormat.getImageableX();
			double y0 = pageFormat.getImageableY();

			double w0 = pageFormat.getImageableWidth();
			double h0 = pageFormat.getImageableHeight();

			double w1 = getWidth();
			double h1 = getHeight();

			double scale;

			if (w0 / w1 < h0 / h1) {
				scale = w0 / w1;
			} else {
				scale = h0 /h1;
			}

			g2d.translate(x0, y0);
			g2d.scale(scale, scale);

			// Turn off double buffering
			paint(g2d);
			// Turn double buffering back on
			return(PAGE_EXISTS);
		}
	}

    //********************************************************************
	// MutableTreeListener interface
	//********************************************************************

	public void treeChanged(Tree tree) { repaint(); }

	//********************************************************************
	// MutableTaxonListListener interface
	//********************************************************************

	public void taxonAdded(TaxonList taxonList, Taxon taxon) { repaint(); }
	public void taxonRemoved(TaxonList taxonList, Taxon taxon) { repaint(); }
	public void taxaChanged(TaxonList taxonList) { repaint(); }

	/**
	*	Add a plot listener
	*/
	public void addListener(Listener listener) {

		listeners.add(listener);
	}

	/**
	 * Tells tree listeners that a node has been clicked.
	 */
	protected void fireNodeClickedEvent(int node) {

		for (int i=0; i < listeners.size(); i++) {
			Listener listener = listeners.elementAt(i);
			listener.nodeClicked(node);
		}
	}

	// Listeners

	private final java.util.Vector<Listener> listeners = new java.util.Vector<Listener>();

	public interface Listener {

		public void nodeClicked(int node);

	}

	public class Adaptor implements Listener {

		public void nodeClicked(int node) { }

	}

	public class MListener extends MouseAdapter {

		public void mouseClicked(MouseEvent me) {

			int node = treePainter.findNodeAtPoint(me.getPoint());

			fireNodeClickedEvent(node);
		}
	}
}
