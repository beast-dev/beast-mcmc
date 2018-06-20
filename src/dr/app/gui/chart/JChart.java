/*
 * JChart.java
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

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.*;
import java.util.ArrayList;
import java.util.List;

public class JChart extends JPanel {

    private static final boolean DEBUG = true;
    /**
     *
     */
    private static final long serialVersionUID = -7064065852204509247L;
    private Axis yAxis, xAxis;
    private List<Plot> plots = new ArrayList<Plot>();

    private Paint plotBackgroundPaint = Color.white;

    private Stroke originStroke = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL);
    private Paint originPaint = Color.lightGray;

    private Stroke frameStroke = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER);
    private Paint framePaint = Color.black;

    private Stroke axisStroke = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER);
    private Paint axisPaint = Color.black;

    private Font labelFont = new Font("Helvetica", Font.PLAIN, 12);
    private Paint labelPaint = Color.black;

    private Rectangle2D plotBounds = null;

    private Rectangle2D dragRectangle = null;

    private double majorTickSize = 4;
    private double minorTickSize = 2;

    private double xTickLabelOffset;
    private double yTickLabelOffset;

    private boolean showLegend = false;
    private int legendAlignment = SwingConstants.NORTH_EAST;

    private double xScale, yScale, xOffset, yOffset;

    private double aspectRatio;
    private boolean useAspectRatio = false;

    public JChart(Axis xAxis, Axis yAxis) {
        this(null, xAxis, yAxis);
    }

    public JChart(Axis xAxis, Axis yAxis, double aspectRatio) {
        this(null, xAxis, yAxis);
        this.aspectRatio = aspectRatio;
        this.useAspectRatio = true;
    }

    public JChart(Plot plot, Axis xAxis, Axis yAxis) {

        setOpaque(false);

        this.xAxis = xAxis;
        this.yAxis = yAxis;

        if (plot != null) {
            addPlot(plot);
        }

        this.aspectRatio = 1.0;
        this.useAspectRatio = false;

        addMouseListener(new MListener());
        addMouseMotionListener(new MMListener());

    }

    public void setXAxis(Axis xAxis) {
        this.xAxis = xAxis;
        for (Object plot : plots) {
            Plot p = (Plot) plot;
            p.setAxes(xAxis, yAxis);
        }
        recalibrate();
        repaint();
    }

    public Axis getXAxis() {
        return xAxis;
    }

    public void setYAxis(Axis yAxis) {
        this.yAxis = yAxis;
        for (Plot p : plots) {
            p.setAxes(xAxis, yAxis);
        }
        recalibrate();
        repaint();
    }

    public Axis getYAxis() {
        return yAxis;
    }

    public void setFontSize(int size) {
        labelFont = new Font("Helvetica", Font.PLAIN, size);
    }

    public void addPlot(Plot plot) {
        plot.setAxes(xAxis, yAxis);
        plot.setChart(this);
        plot.setPlotNumber(plots.size());
        plots.add(plot);
        recalibrate();
        repaint();
    }

    public void removePlot(Plot plot) {
        plots.remove(plot);
        xAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
        yAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
        for (Plot p : plots) {
            p.setAxes(xAxis, yAxis);
        }
        recalibrate();
        repaint();
    }

    public void removeAllPlots() {
        plots.clear();
        xAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
        yAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
        recalibrate();
        repaint();
    }

    public List<Plot> getPlots() {
        return plots;
    }

    public void setPlots(List<Plot> plots) {
        this.plots = plots;
    }

    public int getPlotCount() {
        return plots.size();
    }

    public Plot getPlot(int index) {
        return plots.get(index);
    }

    public Rectangle2D getDragRectangle() {
        return dragRectangle;
    }

    public void setDragRectangle(Rectangle2D dragRectangle) {
        this.dragRectangle = dragRectangle;
        repaint();
    }

    public void selectPoints(Rectangle2D dragRectangle, boolean addToSelection) {
        for (Plot plot : plots) {
            plot.selectPoints(dragRectangle, addToSelection);
        }
        repaint();
    }

    public void clearSelection() {
        for (Plot plot : plots) {
            plot.clearSelection();
        }
        repaint();
    }

    private void resetPlots() {
//		xAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
//		yAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
        for (Plot p : plots) {
            p.resetAxes();
        }
    }

    /**
     *	Set origin style. Use a stroke of null to not draw origin
     */
    public void setOriginStyle(Stroke originStroke, Paint originPaint) {
        this.originStroke = originStroke;
        this.originPaint = originPaint;
    }

    /**
     *	Set axis style.
     */
    public void setAxisStyle(Stroke axisStroke, Paint axisPaint) {
        this.axisStroke = axisStroke;
        this.axisPaint = axisPaint;
    }

    /**
     *	Get axis stroke.
     */
    public Stroke getAxisStroke() {
        return axisStroke;
    }

    /**
     *	Get axis paint.
     */
    public Paint getAxisPaint() {
        return axisPaint;
    }

    /**
     *	Set frame style. Use a stroke of null to not draw frame
     */
    public void setFrameStyle(Stroke frameStroke, Paint framePaint) {
        this.frameStroke = frameStroke;
        this.framePaint = framePaint;
    }

    /**
     *	Set label style.
     */
    public void setLabelStyle(Font labelFont, Paint labelPaint) {
        this.labelFont = labelFont;
        this.labelPaint = labelPaint;
    }

    /**
     *	Get label font.
     */
    public Font getLabelFont() {
        return labelFont;
    }

    /**
     *	Get label paint.
     */
    public Paint getLabelPaint() {
        return labelPaint;
    }

    /**
     *	Set legend on/off.
     */
    public void setShowLegend(boolean showLegend) {
        this.showLegend = showLegend;
        repaint();
    }

    /**
     *	Set legend alignment (use SwingConstants.NORTH_EAST etc.).
     */
    public void setLegendAlignment(int alignment) {
        this.legendAlignment = alignment;
        repaint();
    }

    public Dimension getMinimumSize() {
        return new Dimension(128, 128);
    }

    public Rectangle2D getPlotBounds() {
        return plotBounds;
    }

    public double getMajorTickSize() {
        return majorTickSize;
    }

    public double getMinorTickSize() {
        return minorTickSize;
    }

    public double getXTickLabelOffset() {
        return xTickLabelOffset;
    }

    public double getYTickLabelOffset() {
        return yTickLabelOffset;
    }

    private boolean calibrated = true;
    public void recalibrate() { calibrated = false; }
    protected void calibrate(Graphics2D g2, Dimension size) {
        resetPlots();
    }
    protected boolean hasContents() {
        return plots.size() > 0;
    }

    public void paintComponent(Graphics g) {

        if (!hasContents()) return;

        Graphics2D g2 = (Graphics2D)g;
        Dimension size = getSize();

//        g2.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        try {

            if (!calibrated) {
                calibrate(g2, size);
                calibrated = true;
            }

            g2.setFont(labelFont);

            double tickLabelHeight = g2.getFontMetrics().getHeight();

            xTickLabelOffset = g2.getFontMetrics().getAscent();
            yTickLabelOffset = g2.getFontMetrics().getAscent() / 2;

            double maxXTickLabelWidth = getMaxTickLabelWidth(g2, xAxis);
            double maxYTickLabelWidth = getMaxTickLabelWidth(g2, yAxis);

            double w = size.width - (majorTickSize * 1.25) - maxYTickLabelWidth - (maxXTickLabelWidth / 2);
            double h = size.height - yTickLabelOffset - (majorTickSize * 1.25) - tickLabelHeight;

            //store current transformation in case of a translate transform
            AffineTransform orig = g2.getTransform();
            if (useAspectRatio) {
                double newH = w / aspectRatio;
                if (newH > h) {
                    w = h * aspectRatio;
                } else {
                    h = newH;
                }
                if (size.width > size.height) {
                    g2.translate((size.width - size.height)/2.0, 0.0);
                }
            }

            plotBounds = new Rectangle2D.Double((majorTickSize * 1.25) + maxYTickLabelWidth, yTickLabelOffset, w, h);

            xOffset = plotBounds.getX();
            yOffset = plotBounds.getMaxY();
            xScale = w / (xAxis.transform(xAxis.getMaxAxis()) - xAxis.transform(xAxis.getMinAxis()));
            yScale = -h / (yAxis.transform(yAxis.getMaxAxis()) - yAxis.transform(yAxis.getMinAxis()));

            g2.setPaint(plotBackgroundPaint);
            g2.fill(plotBounds);
            g2.setClip(plotBounds);

            Stroke oldStroke = g2.getStroke();

            if (originStroke != null && originPaint != null) {
                double minX = xAxis.getMinAxis();
                double maxX = xAxis.getMaxAxis();
                double minY = yAxis.getMinAxis();
                double maxY = yAxis.getMaxAxis();

                g2.setPaint(originPaint);
                g2.setStroke(originStroke);

                if (minX<0.0 && maxX>0.0) {
                    Line2D line = new Line2D.Double(transformX(0.0), transformY(minY),
                            transformX(0.0), transformY(maxY));
                    g2.draw(line);
                }

                if (minY<0.0 && maxY>0.0) {
                    Line2D line = new Line2D.Double(transformX(minX), transformY(0.0),
                            transformX(maxX), transformY(0.0));
                    g2.draw(line);
                }
            }
            g2.setStroke(oldStroke);

            paintContents(g2);

            if (showLegend) {
                paintLegend(g2);
            }

            if (dragRectangle != null) {
                g2.setPaint(new Color(128, 128, 128, 128));
                g2.fill(dragRectangle);
            }

            g2.setClip(null);

            paintFrame(g2);

            paintAxis(g2, xAxis, true);
            paintAxis(g2, yAxis, false);

            //restore transformation in any case
            g2.setTransform(orig);

        } catch (ChartRuntimeException cre) {
            if (DEBUG) {
                System.err.println(cre.getStackTrace());
            }
            // ignore (just won't paint chart)
        }

    }

    protected void paintContents(Graphics2D g2) {
        int plotNumber = 0;
        for (Plot plot : plots) {
            plot.paintPlot(g2, xScale, yScale, xOffset, yOffset, plotNumber, plots.size());
            plotNumber += 1;
        }
    }

    protected void paintFrame(Graphics2D g2) {
        g2.setPaint(framePaint);     
        g2.setStroke(frameStroke);
        g2.draw(plotBounds);
    }


    /**
     *	Transform a chart co-ordinates into a drawing co-ordinates
     */
    protected double transformX(double value) {
        return ((xAxis.transform(value) - xAxis.transform(xAxis.getMinAxis())) * xScale) + xOffset;
    }

    /**
     *	Transform a chart co-ordinates into a drawing co-ordinates
     */
    protected double transformY(double value) {
        return ((yAxis.transform(value) - yAxis.transform(yAxis.getMinAxis())) * yScale) + yOffset;
    }

    /*
     *	Transform a drawing co-ordinate into a chart co-ordinate
      * /
     private double untransformX(double value) {

         return xAxis.untransform(
                         xAxis.transform(xAxis.getMinAxis()) + ((value - xOffset) / xScale));
     }*/

    /*
     *	Transform a drawing co-ordinate into a chart co-ordinate
      * /
     private double untransformY(double value) {

         return yAxis.untransform(
                         yAxis.transform(yAxis.getMinAxis()) + ((value - yOffset) / yScale));
     }*/

    protected void paintLegend(Graphics2D g2)
    {
        float width = 0;
        int itemCount = 0;

        for (int i = 0; i < getPlotCount(); i++) {
            String name = getPlot(i).getName();
            if (name != null) {
                float w = (float)g2.getFontMetrics().stringWidth(name);
                if (width < w) {
                    width = w;
                }
                itemCount++;
            }
        }

        if (width == 0) return; // no plots have names so just return

        width += 32;
        float height = (float)(g2.getFontMetrics().getAscent() + 8) * itemCount;

        float x, y;

        if (legendAlignment == SwingConstants.NORTH_WEST ||
                legendAlignment == SwingConstants.WEST ||
                legendAlignment == SwingConstants.SOUTH_WEST) {
            x = (float)(plotBounds.getX() + 8);
        } else if (legendAlignment == SwingConstants.NORTH_EAST ||
                legendAlignment == SwingConstants.EAST ||
                legendAlignment == SwingConstants.SOUTH_EAST) {
            x = (float)(plotBounds.getX() + plotBounds.getWidth() - width - 8);
        } else { // centered
            x = (float)(plotBounds.getX() + ((plotBounds.getWidth() - width) / 2));
        }

        if (legendAlignment == SwingConstants.NORTH_WEST ||
                legendAlignment == SwingConstants.NORTH ||
                legendAlignment == SwingConstants.NORTH_EAST) {
            y = (float)(plotBounds.getY() + 8);
        } else if (legendAlignment == SwingConstants.SOUTH_EAST ||
                legendAlignment == SwingConstants.SOUTH ||
                legendAlignment == SwingConstants.SOUTH_WEST) {
            y = (float)(plotBounds.getY() + plotBounds.getHeight() - height - 8);
        } else { // centered
            y = (float)(plotBounds.getY() + ((plotBounds.getHeight() - height) / 2));
        }

        Rectangle2D legendBounds = new Rectangle2D.Float(x, y, width, height);
        g2.setPaint(framePaint);
        g2.setStroke(frameStroke);
        g2.draw(legendBounds);

        x += 8;

        float iy = 8 + g2.getFontMetrics().getAscent();
        y += g2.getFontMetrics().getAscent() + 4;
        for (int i = 0; i < getPlotCount(); i++) {
            Plot plot = getPlot(i);
            String name = plot.getName();
            if (name != null) {

                g2.setPaint(plot.getLineColor());
                g2.fill(new Rectangle2D.Float(x, y - 8, 8, 8));

                g2.setPaint(framePaint);
                g2.drawString(name, x + 16, y);

                y += iy;
            }
        }

    }

    /**
     *	Get the maximum width of the labels of an axis
     */
    protected double getMaxTickLabelWidth(Graphics2D g2, Axis axis)
    {
        String label;
        double width;
        double maxWidth = 0;

        if (axis.getLabelFirst()) { // Draw first minor tick as a major one (with a label)
            label = axis.format(axis.getMinorTickValue(0, -1));
            width = g2.getFontMetrics().stringWidth(label);
            if (maxWidth < width)
                maxWidth = width;
        }
        int n = axis.getMajorTickCount();
        for (int i = 0; i < n; i++) {
            double value = axis.getMajorTickValue(i);
            label = axis.format(value);
            width = g2.getFontMetrics().stringWidth(label);
            if (maxWidth < width)
                maxWidth = width;
        }
        if (axis.getLabelLast()) { // Draw first minor tick as a major one (with a label)
            label = axis.format(axis.getMinorTickValue(0, n - 1));
            width = g2.getFontMetrics().stringWidth(label);
            if (maxWidth < width)
                maxWidth = width;
        }

        return maxWidth;
    }

    protected void paintAxis(Graphics2D g2, Axis axis, boolean horizontalAxis) {
        if (axis.getIsDiscrete()) {
            int n1 = axis.getMajorTickCount();
            int n2, i, j;

            n2 = axis.getMinorTickCount(-1);

            for (i = 0; i < n1; i++) {
                paintMajorTick(g2, axis.getMajorTickValue(i), axis.getMajorTickLabel(i), horizontalAxis);
                n2 = axis.getMinorTickCount(i);

                if (i == (n1 - 1) && axis.getLabelLast()) { // Draw last minor tick as a major one

                    paintMajorTick(g2, axis.getMinorTickValue(0, i), axis.format((int) axis.getMinorTickValue(0, i)), horizontalAxis);

                    for (j = 1; j < n2; j++) {
                        paintMinorTick(g2, axis.getMinorTickValue(j, i), horizontalAxis);
                    }
                } else {

                    for (j = 0; j < n2; j++) {
                        paintMinorTick(g2, axis.getMinorTickValue(j, i), horizontalAxis);
                    }
                }
            }
        } else {
            int n1 = axis.getMajorTickCount();
            int n2, i, j;

            n2 = axis.getMinorTickCount(-1);
            if (axis.getLabelFirst()) { // Draw first minor tick as a major one (with a label)

                paintMajorTick(g2, axis.getMinorTickValue(0, -1), horizontalAxis);

                for (j = 1; j < n2; j++) {
                    paintMinorTick(g2, axis.getMinorTickValue(j, -1), horizontalAxis);
                }
            } else {

                for (j = 0; j < n2; j++) {
                    paintMinorTick(g2, axis.getMinorTickValue(j, -1), horizontalAxis);
                }
            }

            for (i = 0; i < n1; i++) {

                paintMajorTick(g2, axis.getMajorTickValue(i), horizontalAxis);
                n2 = axis.getMinorTickCount(i);

                if (i == (n1 - 1) && axis.getLabelLast()) { // Draw last minor tick as a major one

                    paintMajorTick(g2, axis.getMinorTickValue(0, i), horizontalAxis);

                    for (j = 1; j < n2; j++) {
                        paintMinorTick(g2, axis.getMinorTickValue(j, i), horizontalAxis);
                    }
                } else {

                    for (j = 0; j < n2; j++) {
                        paintMinorTick(g2, axis.getMinorTickValue(j, i), horizontalAxis);
                    }
                }
            }
        }
    }

    protected String getXAxisLabel(double value) {
        return xAxis.format(value);
    }

    protected String getYAxisLabel(double value) {
        return yAxis.format(value);
    }

//    protected void paintMajorTick(Graphics2D g2, double value, String label, boolean horizontalAxis) {
//        g2.setPaint(getAxisPaint());
//        g2.setStroke(getAxisStroke());
//
//        if (label == null) {
//            label = " ";
//        }
//
//        if (horizontalAxis) {
//            double pos = transformX(value);
//
//            Line2D line = new Line2D.Double(pos, getPlotBounds().getMaxY(), pos, getPlotBounds().getMaxY() + getMajorTickSize());
//            g2.draw(line);
//
//            g2.setPaint(getLabelPaint());
//            double width = g2.getFontMetrics().stringWidth(label);
//
//            if (label == null) {
//                label = getXAxis().format(value);
//            }
//            g2.drawString(label, (float) (pos - (width / 2)), (float) (getPlotBounds().getMaxY() + (getMajorTickSize() * 1.25) + getXTickLabelOffset()));
//
//        } else {
//            double pos = transformY(value);
//
//            Line2D line = new Line2D.Double(getPlotBounds().getMinX(), pos, getPlotBounds().getMinX() - getMajorTickSize(), pos);
//            g2.draw(line);
//
//            g2.setPaint(getLabelPaint());
//            double width = g2.getFontMetrics().stringWidth(label);
//
//            if (label == null) {
//                label = getYAxis().format(value);
//            }
//            g2.drawString(label, (float)(getPlotBounds().getMinX() - width - (getMajorTickSize() * 1.25)), (float)(pos + getYTickLabelOffset()));
//        }
//    }

    protected void paintMajorTick(Graphics2D g2, double value, boolean horizontalAxis) {
        paintMajorTick(g2, value, null, horizontalAxis);
    }

    protected void paintMajorTick(Graphics2D g2, double value, String label, boolean horizontalAxis)
    {
        g2.setPaint(axisPaint);
        g2.setStroke(axisStroke);


        if (horizontalAxis) {
                if (label == null) {
                label = getXAxisLabel(value);
            }

            if (label != null) {
                double pos = transformX(value);

                Line2D line = new Line2D.Double(pos, plotBounds.getMaxY(), pos, plotBounds.getMaxY() + majorTickSize);
                g2.draw(line);

                g2.setPaint(labelPaint);
                double width = g2.getFontMetrics().stringWidth(label);
            /*if (rotateHorizontal) {
                AffineTransform orig = g2.getTransform();
                g2.rotate(Math.toRadians(90), (float) (pos - width / 2), (float) (plotBounds.getMaxY() + (majorTickSize * 1.25) + xTickLabelOffset));
                g2.drawString(label, (float) (pos - (width / 2)), (float) (plotBounds.getMaxY() + (majorTickSize * 1.25) + xTickLabelOffset - width / 2));
                g2.setTransform(orig);
            } else {
                g2.drawString(label, (float) (pos - (width / 2)), (float) (plotBounds.getMaxY() + (majorTickSize * 1.25) + xTickLabelOffset));
            }*/
                g2.drawString(label, (float) (pos - (width / 2)), (float) (plotBounds.getMaxY() + (majorTickSize * 1.25) + xTickLabelOffset));
            }

        } else {
            if (label == null) {
                label = getYAxisLabel(value);
            }

            if (label != null) {
                double pos = transformY(value);

                Line2D line = new Line2D.Double(plotBounds.getMinX(), pos, plotBounds.getMinX() - majorTickSize, pos);
                g2.draw(line);

                g2.setPaint(labelPaint);
                double width = g2.getFontMetrics().stringWidth(label);
                g2.drawString(label, (float) (plotBounds.getMinX() - width - (majorTickSize * 1.25)), (float) (pos + yTickLabelOffset));
            }
        }
    }

    protected void paintMinorTick(Graphics2D g2, double value, boolean horizontalAxis)
    {

        g2.setPaint(axisPaint);
        g2.setStroke(axisStroke);
        if (horizontalAxis) {
            double pos = transformX(value);

            Line2D line = new Line2D.Double(pos, plotBounds.getMaxY(), pos, plotBounds.getMaxY() + minorTickSize);
            g2.draw(line);
        } else {
            double pos = transformY(value);

            Line2D line = new Line2D.Double(plotBounds.getMinX(), pos, plotBounds.getMinX() - minorTickSize, pos);
            g2.draw(line);
        }
    }

//	private Point2D currentPoint = null;

    public class MMListener extends MouseMotionAdapter {

        public void mouseMoved(MouseEvent me) {
//			currentPoint = me.getPoint();
        }

        public void mouseDragged(MouseEvent me) {

//			selectionEnd = boundedPoint(me.getPoint(), getSize());
//			currentPoint = realPoint(selectionEnd);

//			repaint(200);
        }
    }

    public class MListener extends MouseAdapter {

        Point2D start, finish;

        public void mouseExited(MouseEvent me) {
//			currentPoint = null;
        }

        public void mousePressed(MouseEvent me) {
        }


        public void mouseReleased(MouseEvent me) {
        }

        public void mouseClicked(MouseEvent me) {

            if (plotBounds != null && plotBounds.contains(me.getPoint())) {
                for (Plot plot : plots) {
                    plot.pointClicked(me.getPoint(), me.isShiftDown());
                }


            }
        }
    }
}
