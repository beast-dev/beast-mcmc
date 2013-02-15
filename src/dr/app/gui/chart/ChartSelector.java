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
        if (chart == null) {
            return;
        }

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
        chart.selectPoints(new Rectangle2D.Double(x1, y1, x2 - x1, y2 - y1));

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