/*
 * ChartSelector.java
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

package dr.app.gui.chart;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

/**
 * @author Andrew Rambaut
 * @version $Id: TreePaneSelector.java 822 2007-10-26 13:50:26Z rambaut $
 */
public class ChartSelector implements MouseListener, MouseMotionListener {

    public ChartSelector(JChart chart) {
        this.chart = chart;
        chart.addMouseListener(this);
        chart.addMouseMotionListener(this);

    }

    public boolean isCrossHairCursor() {
        return crossHairCursor;
    }

    public void setCrossHairCursor(boolean crossHairCursor) {
        this.crossHairCursor = crossHairCursor;
    }

    private void setupCursor() {
        if (crossHairCursor) {
            chart.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
        } else {
            chart.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
        }

        chart.repaint();
    }

    public void mouseClicked(MouseEvent mouseEvent) {
//        Node selectedNode = treePane.getNodeAt((Graphics2D) treePane.getGraphics(), mouseEvent.getPoint());
//        if (!mouseEvent.isShiftDown()) {
//            treePane.clearSelection();
//        }
//
//        SelectionMode mode = selectionMode;
//        if (mouseEvent.isAltDown()) {
//            if (mode == SelectionMode.NODE) {
//                mode = SelectionMode.CLADE;
//            } else if (mode == SelectionMode.CLADE) {
//                mode = SelectionMode.NODE;
//            }
//        }
//
//        treePane.addSelectedNode(selectedNode);

    }

    public void mousePressed(MouseEvent mouseEvent) {
        // This is used for dragging in combination with mouseDragged
        // in the MouseMotionListener, below.
        dragPoint = new Point2D.Double(mouseEvent.getPoint().getX(), mouseEvent.getPoint().getY());
        if (!mouseEvent.isShiftDown()) {
            chart.clearSelection();
        }

    }

    public void mouseReleased(MouseEvent mouseEvent) {
        if (chart == null) {
            return;
        }

        double x1 = Math.min(dragPoint.getX(), mouseEvent.getPoint().getX());
        double y1 = Math.min(dragPoint.getY(), mouseEvent.getPoint().getY());
        double x2 = Math.max(dragPoint.getX(), mouseEvent.getPoint().getX());
        double y2 = Math.max(dragPoint.getY(), mouseEvent.getPoint().getY());
        chart.selectPoints(new Rectangle2D.Double(x1, y1, x2 - x1, y2 - y1), mouseEvent.isShiftDown());

        chart.setDragRectangle(null);
    }

    public void mouseEntered(MouseEvent mouseEvent) {
    }

    public void mouseExited(MouseEvent mouseEvent) {
    }

    public void mouseMoved(MouseEvent mouseEvent) {
    }

    public void mouseDragged(MouseEvent mouseEvent) {

        if (dragPoint == null) {
            return;
        }

        double x1 = Math.min(dragPoint.getX(), mouseEvent.getPoint().getX());
        double y1 = Math.min(dragPoint.getY(), mouseEvent.getPoint().getY());
        double x2 = Math.max(dragPoint.getX(), mouseEvent.getPoint().getX());
        double y2 = Math.max(dragPoint.getY(), mouseEvent.getPoint().getY());
        chart.setDragRectangle(new Rectangle2D.Double(x1, y1, x2 - x1, y2 - y1));
//        chart.scrollPointToVisible(mouseEvent.getPoint());
    }

    private JChart chart;
    private Point2D dragPoint = null;

    private boolean crossHairCursor = false;
}