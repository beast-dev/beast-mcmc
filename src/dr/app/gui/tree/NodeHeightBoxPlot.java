/*
 * NodeHeightBoxPlot.java
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
import dr.evolution.tree.NodeRef;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Line2D;

/**
 * Created by IntelliJ IDEA.
 * User: alexei
 * Date: Dec 3, 2004
 * Time: 10:09:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class NodeHeightBoxPlot implements NodeDecorator {

    public NodeHeightBoxPlot() {

        this(new BasicStroke(1.0f), new BasicStroke(0.5f), true);
    }

    public NodeHeightBoxPlot(Stroke lineStroke, Stroke whiskerStroke, boolean fill) {
        this.lineStroke = lineStroke;
        this.whiskerStroke = whiskerStroke;
        this.fill = fill;
    }

    public boolean isDecoratable(Tree tree, NodeRef node) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void decorateNode(Tree tree, NodeRef node, Graphics2D g2, CoordinateTransform transform) {

        // look for node height statistics
        Double mean = (Double)tree.getNodeAttribute(node, "nodeHeight.mean");
        if (mean != null) {
            Double hpdUpper = (Double)tree.getNodeAttribute(node, "nodeHeight.hpdUpper");
            Double hpdLower = (Double)tree.getNodeAttribute(node, "nodeHeight.hpdLower");
            Double min = (Double)tree.getNodeAttribute(node, "nodeHeight.min");
            Double max = (Double)tree.getNodeAttribute(node, "nodeHeight.max");

            // plot height statistics as box plot
            double meanX = transform.xCoordinate(mean.doubleValue());
            double minX = transform.xCoordinate(min.doubleValue());
            double maxX = transform.xCoordinate(max.doubleValue());
            double upperX = transform.xCoordinate(hpdUpper.doubleValue());
            double lowerX = transform.xCoordinate(hpdLower.doubleValue());
            //System.out.println(upperX + " " + lowerX);

            double y = transform.yCoordinate(tree, node);

            g2.setStroke(lineStroke);
            if (fill) {
                g2.setColor(Color.white);
                g2.fill(new Rectangle2D.Double(upperX, y-5,lowerX-upperX,10));
            }
            g2.setColor(Color.gray);
            g2.draw(new Rectangle2D.Double(upperX, y-5,lowerX-upperX,10));
            g2.draw(new Line2D.Double(meanX, y-5,meanX,y+5));

            g2.setStroke(whiskerStroke);
            g2.draw(new Line2D.Double(lowerX, y,minX,y));
            g2.draw(new Line2D.Double(maxX, y,upperX,y));
        }

	}

    Stroke lineStroke = new BasicStroke(1.0f);
    Stroke whiskerStroke = new BasicStroke(0.5f);
    boolean fill = true;
}
