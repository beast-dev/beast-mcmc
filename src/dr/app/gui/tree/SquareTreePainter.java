/*
 * SquareTreePainter.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;

import java.awt.*;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Map;

public class SquareTreePainter implements TreePainter {

	protected Stroke lineStroke = new BasicStroke((float)2.0);
	protected Paint linePaint = Color.black;

	protected Stroke hilightStroke = new BasicStroke((float)2.0);
	protected Paint hilightPaint = Color.blue;

	protected int maxFontSize = 12;

	protected Font labelFont = new Font("Helvetica", Font.PLAIN, 12);
	protected Paint labelPaint = Color.black;

	protected Font hilightLabelFont = new Font("Helvetica", Font.PLAIN, 12);
	protected Paint hilightLabelPaint = Color.blue;

    final int boxPlotSize = 3;

	public SquareTreePainter() {
		this.rememberYPositions = false;
	}

	public SquareTreePainter(boolean rememberYPositions) {
		this.rememberYPositions = rememberYPositions;
	}

	public SquareTreePainter(TaxonList taxonList) {

		if (taxonList != null) {

			// initialize taxa starting y positions based on the given taxon list
			double y = 0.5;
			for (int i = 0; i < taxonList.getTaxonCount(); i++) {
				yPositionMap.put(taxonList.getTaxonId(i), y);
				y += 1.0;
			}
		}
		this.rememberYPositions = true;
	}

	/**
	*	Set line style
	*/
	public void setLineStyle(Stroke lineStroke, Paint linePaint) {
		this.lineStroke = lineStroke;
		this.linePaint = linePaint;
	}

	/**
	*	Set line style
	*/
	public void setLinePaint(Paint linePaint) {
		this.linePaint = linePaint;
	}


	/**
	*	Set hilight style
	*/
	public void setHilightStyle(Stroke hilightStroke, Paint hilightPaint) {
		this.hilightStroke = hilightStroke;
		this.hilightPaint = hilightPaint;
	}

	public void setFontSize(int size) {
		maxFontSize = size;
	}

	/**
	 *	Set label style.
	 */
	public void setLabelStyle(Font labelFont, Paint labelPaint) {
		this.labelFont = labelFont;
		this.labelPaint = labelPaint;
	}

	/**
	 *	Set hilight label style.
	 */
	public void setHilightLabelStyle(Font hilightLabelFont, Paint hilightLabelPaint) {
		this.hilightLabelFont = hilightLabelFont;
		this.hilightLabelPaint = hilightLabelPaint;
	}

	public void setUserDefinedHeight(double height) {
		this.userDefinedHeight = height;
	}

	public void drawLabels(boolean drawLabels) {
		this.drawLabels = drawLabels;
	}

	public void drawHorizontals(boolean drawHorizontals) {
		this.drawHorizontals = drawHorizontals;
	}

	public void drawVerticals(boolean drawVerticals) {
		this.drawVerticals = drawVerticals;
	}

	/**
	 * Do the actual painting.
	 */
	public void paintTree(Graphics2D g2, Dimension size, Tree tree) {
		if (tree == null) return;

		int n = tree.getNodeCount();
		if (nodeRectVert == null || nodeRectVert.length != n) {
			nodeRectVert = new Rectangle2D[n];
			nodeRectHoriz = new Rectangle2D[n];
		}

		scaleY = ((double)size.height) / (tree.getExternalNodeCount());

		double maxLabelHeight = scaleY;

		int fontSize = maxFontSize + 1;
		do {
			fontSize --;
			labelFont = new Font("Helvetica", Font.PLAIN, fontSize);
			g2.setFont(labelFont);
		} while (fontSize > 1 &&
			g2.getFontMetrics().getAscent()
//			+ g2.getFontMetrics().getDescent()
			> maxLabelHeight);

		hilightLabelFont = new Font("Helvetica", Font.PLAIN, fontSize);

		double maxLabelWidth = getMaxLabelWidth(g2, tree);

		currentY = 0.5;

		treeHeight = tree.getNodeHeight(tree.getRoot());
		double height;
		if (userDefinedHeight < 0.0) {
			height = treeHeight;
		} else {
			height = userDefinedHeight;
		}

		scaleX = ((double)size.width - 4 - maxLabelWidth) / (height * 1.02);

        //AffineTransform transform = new AffineTransform(-scaleX, 0, treeHeight*1.02*scaleX, 0, scaleY, 0);
        //g2.transform(transform);

		//paintNode(g2, tree, tree.getRoot(), 0.0, (height * 1.02)-treeHeight, false);
        currentY = 0.5;
		paintBoxPlot(g2, tree, tree.getRoot(), treeHeight, false);
        currentY = 0.5;
		paintBoxPlot(g2, tree, tree.getRoot(), treeHeight, false);
        currentY = 0.5;
		paintNode(g2, tree, tree.getRoot(), (treeHeight * 1.02), treeHeight, false);
	}

	/**
	 * Paint a node.
     * @param x0 the height of the parent node
     * @param x1 the height of the node
	 */
	private double paintNode(Graphics2D g2, Tree tree, NodeRef node,
								double x0, double x1, boolean hilight) {

		double y;

		double ix0 = convertX(x0);
		double ix1 = convertX(x1);
		double iy;

		if (tree.getNodeAttribute(node, "selected") != null) {
			hilight = true;
		}

		//Color color = null;
		//Object colObj = tree.getNodeAttribute(node, "colour");
		//if (colObj != null) {
		//	if (colObj instanceof Color) {
		//		color = (Color)colObj;
		//	}
		//}

		if (tree.isExternal(node)) {

			if (rememberYPositions) {
				// remember the y positions of taxa that you have seen before... AD
				String taxonId = tree.getNodeTaxon(node).getId();
				Double pos = yPositionMap.get(taxonId);
				if (pos != null) {
					y = pos;
				} else {
					y = currentY;
					currentY += 1.0;
					yPositionMap.put(taxonId, y);
				}
			} else {
				y = currentY;
				currentY += 1.0;
			}

			if (hilight) {
				g2.setPaint(hilightLabelPaint);
				g2.setFont(hilightLabelFont);
			} else {
				g2.setPaint(labelPaint);
				g2.setFont(labelFont);
			}


			String label = tree.getTaxonId(node.getNumber());
			double labelWidth = g2.getFontMetrics().stringWidth(label);
			double labelHeight = g2.getFontMetrics().getAscent();
			double labelOffset = labelHeight / 2;

			iy = convertY(y);

			if (label != null && label.length() > 0 && drawLabels) {
				g2.drawString(label, (float)(ix1 + 4), (float)(iy + labelOffset));
			}


			nodeRectVert[node.getNumber()] =
				new Rectangle.Double(ix1 + 4, iy, labelWidth, labelHeight);

			if (hilight) {
				g2.setPaint(hilightPaint);
				g2.setStroke(hilightStroke);
			} else {
                // use tree color attribute if set
				if (colorAttribute != null) {
                    Paint c = (Color)tree.getNodeAttribute(node,colorAttribute);
                    if (c == null) c = linePaint;
                    g2.setPaint(c);
                } else {
                    g2.setPaint(linePaint);
                }
				if (lineAttribute != null) {
                    Stroke stroke = (Stroke)tree.getNodeAttribute(node,lineAttribute);
                    if (stroke == null) stroke = lineStroke;
                    g2.setStroke(stroke);
                } else g2.setStroke(lineStroke);
			}

		} else {
			double y0, y1;

			NodeRef child = tree.getChild(node, 0);
			double length = tree.getNodeHeight(node) - tree.getNodeHeight(child);

			y0 = paintNode(g2, tree, child, x1, x1-length, hilight);
			y1 = y0;

			for (int i = 1; i < tree.getChildCount(node); i++) {
				child = tree.getChild(node, i);
				length = tree.getNodeHeight(node) - tree.getNodeHeight(child);

				y1 = paintNode(g2, tree, child, x1, x1-length, hilight);
			}

			double iy0 = convertY(y0);
			double iy1 = convertY(y1);

			if (hilight) {
				g2.setPaint(hilightPaint);
				g2.setStroke(hilightStroke);
			} else {
				 // use tree color attribute if set
				if (colorAttribute != null) {
                    Paint c = (Color)tree.getNodeAttribute(node,colorAttribute);
                    if (c == null) c = linePaint;
                    g2.setPaint(c);
                } else {
                    g2.setPaint(linePaint);
                }
				if (lineAttribute != null) {
                    Stroke stroke = (Stroke)tree.getNodeAttribute(node,lineAttribute);
                    if (stroke == null) stroke = lineStroke;
                    g2.setStroke(stroke);
                } else g2.setStroke(lineStroke);
			}

			if (drawHorizontals) {
				Line2D line = new Line2D.Double(ix1, iy0, ix1, iy1);
				g2.draw(line);
			}

			nodeRectVert[node.getNumber()] = new Rectangle.Double(ix1-2, iy0-2, 5, (iy1 - iy0) + 4);

			y = (y1 + y0) / 2;
			iy = convertY(y);

		}

		if (drawVerticals) {
			Line2D line = new Line2D.Double(ix0, iy, ix1, iy);
			g2.draw(line);
		}

		nodeRectHoriz[node.getNumber()] = new Rectangle.Double(ix0-2, iy-2, (ix1 - ix0) + 4, 5);


        if (shapeAttribute != null) {
            Shape shape = (Shape)tree.getNodeAttribute(node,shapeAttribute);
            if (shape != null) {
                Rectangle bounds = shape.getBounds();
                double tx = ix1-bounds.getWidth()/2.0;
                double ty = iy-bounds.getHeight()/2.0;
                g2.translate(tx,ty);
                g2.fill(shape);
                g2.translate(-tx,-ty);
            }
        }

        if (labelAttribute != null) {
            Object label = tree.getNodeAttribute(node,labelAttribute);
            if (label != null) {
                Color c = g2.getColor();
                Font f = g2.getFont();
                Font fsmall = f.deriveFont(f.getSize()-1.0f);
                g2.setFont(fsmall);
                String labelString = label.toString();
                int width = g2.getFontMetrics().stringWidth(labelString);
                g2.setColor(textColor);
                g2.drawString(labelString,(float)(ix1-width-1.0),(float)(iy-2.0));

                // recover color and font
                g2.setColor(c);
                g2.setFont(f);
            }
        }

		return y;
	}

    /**
	 * Paint a box plot of the uncertainty in a node height.
     * @param x1 the height of the node
	 */
	private double paintBoxPlot(Graphics2D g2, Tree tree, NodeRef node, double x1, boolean fill) {

		double y;

		double iy;

        if (tree.isExternal(node)) {

			if (rememberYPositions) {
				// remember the y positions of taxa that you have seen before... AD
				String taxonId = tree.getNodeTaxon(node).getId();
				Double pos = yPositionMap.get(taxonId);
				if (pos != null) {
					y = pos;
				} else {
					y = currentY;
					currentY += 1.0;
					yPositionMap.put(taxonId, y);
				}
			} else {
				y = currentY;
				currentY += 1.0;
			}

            iy = convertY(y);

		} else {
			double y0, y1;

			NodeRef child = tree.getChild(node, 0);
			double length = tree.getNodeHeight(node) - tree.getNodeHeight(child);

			y0 = paintBoxPlot(g2, tree, child, x1-length, fill);
			y1 = y0;

			for (int i = 1; i < tree.getChildCount(node); i++) {
				child = tree.getChild(node, i);
				length = tree.getNodeHeight(node) - tree.getNodeHeight(child);

				y1 = paintBoxPlot(g2, tree, child, x1-length, fill);
			}

			y = (y1 + y0) / 2;
			iy = convertY(y);
		}

        // look for node height statistics
        Double mean = (Double)tree.getNodeAttribute(node, "nodeHeight.mean");
        if (mean != null) {
            if (tree.isRoot(node)) {
                System.out.println(mean.doubleValue());
            }

            Double hpdUpper = (Double)tree.getNodeAttribute(node, "nodeHeight.hpdUpper");
            Double hpdLower = (Double)tree.getNodeAttribute(node, "nodeHeight.hpdLower");
            //Double min = (Double)tree.getNodeAttribute(node, "nodeHeight.min");
            //Double max = (Double)tree.getNodeAttribute(node, "nodeHeight.max");

            // plot height statistics as box plot
            //double meanX = convertX(mean.doubleValue());
            //double minX = convertX(min.doubleValue());
            //double maxX = convertX(max.doubleValue());
            double upperX = convertX(hpdUpper);
            double lowerX = convertX(hpdLower);
            //System.out.println(upperX + " " + lowerX);

            g2.setStroke(lineStroke);
            if (fill) {
                g2.setColor(Color.white);
                g2.fill(new Rectangle2D.Double(upperX, iy-boxPlotSize,lowerX-upperX,2*boxPlotSize));
            }
            g2.setColor(Color.gray);
            g2.draw(new Rectangle2D.Double(upperX, iy-boxPlotSize,lowerX-upperX,2*boxPlotSize));
            //g2.draw(new Line2D.Double(meanX, iy-5,meanX,iy+5));

            //g2.setStroke(whiskerStroke);
            //g2.draw(new Line2D.Double(lowerX, iy,minX,iy));
            //g2.draw(new Line2D.Double(maxX, iy,upperX,iy));
        }

		return y;
	}

	private double convertX(double x) {
        return ((treeHeight*1.02)-x) * scaleX;
        //return x;
    }

	private double convertY(double y) {
        return y * scaleY;
        //return y;
    }

	/**
	 * @return the maximum label width
	 */
	private double getMaxLabelWidth(Graphics2D g, Tree tree) {
		double maxLabelWidth = 0.0;
		for (int i = 0; i < tree.getTaxonCount(); i++) {
			String label = tree.getTaxonId(i);
			double labelWidth = g.getFontMetrics().stringWidth(label);
			if (labelWidth > maxLabelWidth)
				maxLabelWidth = labelWidth;
		}
		return maxLabelWidth;
	}

	/**
	*	Find the node under point. Returns -1 if not found.
	*/
	public final int findNodeAtPoint(Point2D point) {
		if (nodeRectVert == null || nodeRectHoriz == null) return -1;

		int i = 0;
		double x = point.getX();
		double y = point.getY();
		while (i < nodeRectVert.length) {
			if (nodeRectVert[i].contains(x, y) ||
					nodeRectHoriz[i].contains(x, y) )
				return i;
			else
				i++;
		}

		return -1;
	}

    public void setColorAttribute(String s) {
        colorAttribute = s;
    }

    public void setLineAttribute(String s) {
        lineAttribute = s;
    }

    public void setShapeAttribute(String s) {
        shapeAttribute = s;
    }

    public void setLabelAttribute(String s) {
        labelAttribute = s;
    }

	// PRIVATE MEMBERS

	private double scaleX, scaleY, currentY;
    private double treeHeight;

	private Rectangle2D[] nodeRectVert;
	private Rectangle2D[] nodeRectHoriz;

	private final Color textColor = Color.black;
	private final Map<String, Double> yPositionMap = new HashMap<String, Double>();
	private boolean rememberYPositions = false;
	private double userDefinedHeight = -1.0;
	private boolean drawLabels = true;
	private boolean drawHorizontals = true;
	private boolean drawVerticals = true;

    private String colorAttribute = null;
    private String lineAttribute = null;
    private String shapeAttribute = null;
    private String labelAttribute = null;

}
