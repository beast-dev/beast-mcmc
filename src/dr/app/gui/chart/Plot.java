/*
 * Plot.java
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

import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.*;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Plot.java
 * <p/>
 * Description:	Provides an interface for plots. Plots are elements of a
 * chart that renders some data in a particular style. Multiple plots can
 * be added to the same chart to get complex results.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Plot.java,v 1.23 2006/01/03 10:26:20 rambaut Exp $
 */

public interface Plot {

    // These constants are used for automatic scaling to select exactly
    // where the axis starts and stops.
    int NO_MARK = 0;
    int POINT_MARK = 1;
    int CROSS_MARK = 2;
    int PLUS_MARK = 3;
    int CIRCLE_MARK = 4;
    int SQUARE_MARK = 5;
    int DIAMOND_MARK = 6;

    /**
     * Set axes
     */
    void setAxes(Axis xAxis, Axis yAxis);

    /**
     * Resets axis ranges (if new data has been added)
     */
    void resetAxes();

    /**
     * Set data
     */
    void setData(List<Double> xData, List<Double> yData);

    /**
     * Set data
     */
    void setData(Variate.N xData, Variate.N yData);

    /**
     * Set line style
     */
    void setLineStyle(Stroke lineStroke, Paint linePaint);

    /**
     * Set line stroke
     */
    void setLineStroke(Stroke lineStroke);

    /**
     * Set line color
     */
    void setLineColor(Paint linePaint);

    /**
     * Get line color
     */
    Paint getLineColor();

    /**
     * Get line stroke
     */
    Stroke getLineStroke();

    /**
     * Set mark style
     */
    void setMarkStyle(int markType, double markSize, Stroke markStroke,
                      Paint markPaint, Paint markFillPaint);

    /**
     * Set mark style
     */
    void setMarkStyle(Shape mark, Stroke markStroke,
                      Paint markPaint, Paint markFillPaint);

    /**
     * Paint actual plot
     */
    void paintPlot(Graphics2D g2, double xScale, double yScale,
                   double xOffset, double yOffset, int plotNumber, int plotCount);

    /**
     * Set name
     */
    void setName(String name);

    /**
     * Get name
     */
    String getName();

    int getPlotNumber();

    void setPlotNumber(int plotNumber);

    double getXLocation();

    double getYLocation();

    void setLocation(double x, double y);

    /**
     * A point on the plot has been clicked
     */
    void pointClicked(Point2D point, boolean isShiftDown);

    void selectPoint(int index, boolean addToSelection);

    void selectPoints(Rectangle2D dragRectangle, boolean addToSelection);

    void clearSelection();

    void setSelectedPoints(final Collection<Integer> selectedPoints);

    Set<Integer> getSelectedPoints();

    Variate getXData();

    Variate getYData();

    void setChart(JChart chart);

    JChart getChart();

    public interface Listener {

        void pointClicked(double x, double y, boolean isShiftDown);

        void markClicked(int index, double x, double y, boolean isShiftDown);

        void selectionChanged(Set<Integer> selectedPoints);

        void rangeXSelected(double lower, double upper);

        void rangeYSelected(double lower, double upper);

        void rangeXYSelected(double lowerX, double lowerY, double upperX, double upperY);

    }

    public class Adaptor implements Listener {
        public void pointClicked(double x, double y, boolean isShiftDown) {
        }

        public void markClicked(int index, double x, double y, boolean isShiftDown) {
        }

        public void selectionChanged(final Set<Integer> selectedPoints) {
        }

        public void rangeXSelected(double lower, double upper) {
        }

        public void rangeYSelected(double lower, double upper) {
        }

        public void rangeXYSelected(double lowerX, double lowerY, double upperX, double upperY) {
        }
    }

    /**
     * AbstractPlot.java
     * <p/>
     * Description:	An abstract base class for plots
     */
    abstract class AbstractPlot implements Plot {

        private JChart chart = null;

        protected Axis xAxis, yAxis;
        protected Variate.N xData = null;
        protected Variate.N yData = null;

        protected List<Color> colours = null;

        protected Shape mark;

        protected Stroke lineStroke = new BasicStroke(1.5f);
        protected Paint linePaint = Color.black;

        protected Stroke markStroke = new BasicStroke(0.5f);
        protected Paint markPaint = Color.black;
        protected Paint markFillPaint = Color.black;

        protected double xScale, yScale, xOffset, yOffset;

        private String name;

        private int plotNumber = -1, plotCount = -1;
        private double xLocation = 0, yLocation = 0;

        private Set<Integer> selectedPoints = new HashSet<Integer>();

        /**
         * Constructor
         */
        public AbstractPlot() {
        }

        /**
         * Constructor
         */
        public AbstractPlot(Variate.N xData, Variate.N yData) {
            setData(xData, yData);
        }


        /**
         * Constructor
         */
        public AbstractPlot(List<Double> xData, List<Double> yData) {
            setData(xData, yData);
        }

        public JChart getChart() {
            return chart;
        }

        public void setChart(JChart chart) {
            this.chart = chart;
        }

        /**
         * Set data
         */
        public void setData(List<Double> xData, List<Double> yData) {
            Variate.D xd = new Variate.D(xData);
            Variate.D yd = new Variate.D(yData);

            this.xData = xd;
            this.yData = yd;
        }

        /**
         * Set data
         */
        public void setData(Variate.N xData, Variate.N yData) {
            this.xData = xData;
            this.yData = yData;
        }

        /**
         * Set data
         */
        public void setColours(List<Color> colours) {
            this.colours = colours;
        }

        /**
         * Set axes
         */
        public void setAxes(Axis xAxis, Axis yAxis) {
            this.xAxis = xAxis;
            this.yAxis = yAxis;
            setupAxis(xAxis, yAxis, xData, yData);
        }

        /**
         * Resets axis ranges (if new data has been added)
         */
        public void resetAxes() {
            setupAxis(xAxis, yAxis, xData, yData);
        }

        public int getPlotNumber() {
            return plotNumber;
        }

        public void setPlotNumber(int plotNumber) {
            this.plotNumber = plotNumber;
        }

        public double getXLocation() {
            return xLocation;
        }

        public double getYLocation() {
            return yLocation;
        }

        @Override
        public void setLocation(double x, double y) {
            xLocation = x;
            yLocation = y;
        }

        /**
         * Set up the axis with some data
         */
        public void setupAxis(Axis xAxis, Axis yAxis, Variate xData, Variate yData) {
            if (xData != null) {
                if (xAxis instanceof LogAxis) {
                    double minValue = java.lang.Double.POSITIVE_INFINITY;

                    for (int i = 0; i < xData.getCount(); i++) {
                        double value = (Double) xData.get(i);
                        if (value > 0.0 && value < minValue)
                            minValue = value;
                    }

                    xAxis.addRange(minValue, (Double) xData.getMax());
                } else {
                    xAxis.addRange((Double) xData.getMin(), (Double) xData.getMax());
                }
            }
            if (yData != null) {
                if (yAxis instanceof LogAxis) {
                    double minValue = java.lang.Double.POSITIVE_INFINITY;

                    for (int i = 0; i < yData.getCount(); i++) {
                        double value = (Double) yData.get(i);
                        if (value > 0.0 && value < minValue)
                            minValue = value;
                    }

                    yAxis.addRange(minValue, (Double) yData.getMax());
                } else {
                    yAxis.addRange((Double) yData.getMin(), (Double) yData.getMax());
                }
            }
        }

        /**
         * Set line style
         */
        public void setLineStyle(Stroke lineStroke, Paint linePaint) {
            this.lineStroke = lineStroke;
            this.linePaint = linePaint;
        }

        /**
         * Set line stroke
         */
        public void setLineStroke(Stroke lineStroke) {
            this.lineStroke = lineStroke;
        }

        /**
         * Set line color
         */
        public void setLineColor(Paint linePaint) {
            this.linePaint = linePaint;
        }

        public final Paint getLineColor() {
            return linePaint;
        }

        public final Stroke getLineStroke() {
            return lineStroke;
        }

        public final void setName(String name) {
            this.name = name;
        }

        public final String getName() {
            return name;
        }

        /**
         * Set mark style
         */
        public void setMarkStyle(int markType, double markSize, Stroke markStroke,
                                 Paint markPaint, Paint markFillPaint) {

            float w = (float) (markSize / 2.0);

            GeneralPath path;

            switch (markType) {
                case POINT_MARK:
                    path = new GeneralPath();
                    path.moveTo(0, 0);
                    path.lineTo(0, 0);
                    setMarkStyle(path, markStroke, markPaint, markFillPaint);
                    break;
                case CROSS_MARK:
                    path = new GeneralPath();
                    path.moveTo(-w, -w);
                    path.lineTo(w, w);
                    path.moveTo(w, -w);
                    path.lineTo(-w, w);
                    setMarkStyle(path, markStroke, markPaint, markFillPaint);
                    break;
                case PLUS_MARK:
                    path = new GeneralPath();
                    path.moveTo(-w, 0);
                    path.lineTo(w, 0);
                    path.moveTo(0, -w);
                    path.lineTo(0, w);
                    setMarkStyle(path, markStroke, markPaint, markFillPaint);
                    break;
                case CIRCLE_MARK:
                    setMarkStyle(new Ellipse2D.Double(0.0, 0.0, markSize, markSize), markStroke, markPaint, markFillPaint);
                    break;
                case SQUARE_MARK:
                    setMarkStyle(new Rectangle2D.Double(-w, -w, markSize, markSize), markStroke, markPaint, markFillPaint);
                    break;
                case DIAMOND_MARK:
                    path = new GeneralPath();
                    path.moveTo(0, -w);
                    path.lineTo(w, 0);
                    path.lineTo(0, w);
                    path.lineTo(-w, 0);
                    path.closePath();
                    setMarkStyle(path, markStroke, markPaint, markFillPaint);
                    break;
            }


        }

        /**
         * Set mark style
         */
        public void setMarkStyle(Shape mark, Stroke markStroke,
                                 Paint markPaint, Paint markFillPaint) {

            this.mark = mark;
            this.markStroke = markStroke;
            this.markPaint = markPaint;
            this.markFillPaint = markFillPaint;
        }

        public Variate getXData() {
            return xData;
        }

        public Variate getYData() {
            return yData;
        }

        /**
         * Transform a chart co-ordinates into a drawing co-ordinates
         */
        protected double transformX(double value) {
            double tx = xAxis.transform(value + xLocation);
            if (Double.isNaN(tx) || tx == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
            return ((tx - xAxis.transform(xAxis.getMinAxis())) * xScale) + xOffset;
        }

        /**
         * Transform a chart co-ordinates into a drawing co-ordinates
         */
        protected double transformY(double value) {
            double ty = yAxis.transform(value + yLocation);
            if (Double.isNaN(ty) || ty == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
            return ((ty - yAxis.transform(yAxis.getMinAxis())) * yScale) + yOffset;
        }

        /**
         * Transform a drawing co-ordinate into a chart co-ordinate
         */
        protected double untransformX(double value) {

            return xAxis.untransform(
                    xAxis.transform(xAxis.getMinAxis()) + ((value - xOffset) / xScale)) - xLocation;
        }

        /**
         * Transform a drawing co-ordinate into a chart co-ordinate
         */
        protected double untransformY(double value) {

            return yAxis.untransform(
                    yAxis.transform(yAxis.getMinAxis()) + ((value - yOffset) / yScale)) - yLocation;
        }

        /**
         * Draw a line transforming co-ordinates to each axis
         */
        protected void drawLine(Graphics2D g2, double x1, double y1, double x2, double y2) {
            Line2D line = new Line2D.Double(transformX(x1), transformY(y1),
                    transformX(x2), transformY(y2));
            g2.draw(line);
        }

        /**
         * Draw a rectangle directly (no transformation)
         */
        protected void drawRect(Graphics2D g2, float tx1, float ty1, float tx2, float ty2) {
            GeneralPath path = new GeneralPath();
            path.moveTo(tx1, ty1);
            path.lineTo(tx1, ty2);
            path.lineTo(tx2, ty2);
            path.lineTo(tx2, ty1);
            path.closePath();
//			Rectangle2D rect = new Rectangle2D.Double(x, y,	w, h);
            g2.draw(path);
        }

        /**
         * Fill a rectangle directly (no transformation)
         */
        protected void fillRect(Graphics2D g2, float tx1, float ty1, float tx2, float ty2) {

            GeneralPath path = new GeneralPath();
            path.moveTo(tx1, ty1);
            path.lineTo(tx1, ty2);
            path.lineTo(tx2, ty2);
            path.lineTo(tx2, ty1);
            path.closePath();
//			Rectangle2D rect = new Rectangle2D.Double(x, y,	w, h);
            g2.fill(path);
        }

        /**
         * Draw a rectangle transforming co-ordinates to each axis
         */
        protected void drawRect(Graphics2D g2, double x1, double y1, double x2, double y2) {

            float tx1 = (float) transformX(x1);
            float ty1 = (float) transformY(y1);
            float tx2 = (float) transformX(x2);
            float ty2 = (float) transformY(y2);
            drawRect(g2, tx1, ty1, tx2, ty2);
        }

        /**
         * Fill a rectangle transforming co-ordinates to each axis
         */
        protected void fillRect(Graphics2D g2, double x1, double y1, double x2, double y2) {

            float tx1 = (float) transformX(x1);
            float ty1 = (float) transformY(y1);
            float tx2 = (float) transformX(x2);
            float ty2 = (float) transformY(y2);
            fillRect(g2, tx1, ty1, tx2, ty2);
        }

        /**
         * Paint actual plot
         */
        @Override
        public void paintPlot(Graphics2D g2, double xScale, double yScale,
                              double xOffset, double yOffset, int plotNumber, int plotCount) {
            if (xAxis == null || yAxis == null)
                return;

            this.xScale = xScale;
            this.yScale = yScale;
            this.xOffset = xOffset;
            this.yOffset = yOffset;

            this.plotNumber = plotNumber;
            this.plotCount = plotCount;

            if (hasData()) {
                paintData(g2, xData, yData);
            }
        }

        /**
         * Paint data series
         */
        abstract protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData);

        protected boolean hasData() {
            return (xData != null && yData != null && xData.getCount() > 0);
        }

        /**
         * A point on the plot has been clicked
         */
        public void pointClicked(Point2D point, boolean isShiftDown) {

            double x = untransformX(point.getX());
            double y = untransformY(point.getY());

            firePointClickedEvent(x, y, isShiftDown);
        }

        public void selectPoints(final Rectangle2D dragRectangle, final boolean addToSelection) {
            if (dragRectangle == null) {
                return;
            }

            if (!addToSelection) {
                selectedPoints.clear();
            }

            double x0 = untransformX(dragRectangle.getX());
            double y0 = untransformY(dragRectangle.getY() + dragRectangle.getHeight());

            double x1 = untransformX(dragRectangle.getX() + dragRectangle.getWidth());
            double y1 = untransformY(dragRectangle.getY());

            for (int i = 0; i < xData.getCount(); i ++) {
                double x = (Double) xData.get(i);
                double y = (Double) yData.get(i);

                if (x >= x0 && x <= x1 && y >= y0 && y <= y1) {
                    selectedPoints.add(i);
                }
            }

            fireSelectionChanged();
        }

        public void selectPoint(final int index, final boolean addToSelection) {

            if (!addToSelection) {
                selectedPoints.clear();
            }

            selectedPoints.add(index);

            fireSelectionChanged();
        }

        public void clearSelection() {
            selectedPoints.clear();
            fireSelectionChanged();
        }

        public void setSelectedPoints(final Collection<Integer> selectedPoints) {
            this.selectedPoints.clear();
            this.selectedPoints.addAll(selectedPoints);
        }


        public Set<Integer> getSelectedPoints() {
            return selectedPoints;
        }

        // Listeners

        private final java.util.Vector<Listener> listeners = new java.util.Vector<Listener>();

        /**
         * Add a plot listener
         */
        public void addListener(Listener listener) {
            listeners.add(listener);
        }

        /**
         * Tells plot listeners that a point has been clicked.
         */
        protected void firePointClickedEvent(double x, double y, boolean isShiftDown) {
            for (int i = 0; i < listeners.size(); i++) {
                final Listener listener = listeners.elementAt(i);
                listener.pointClicked(x, y, isShiftDown);
            }
        }

        /**
         * Tells plot listeners that a point has been clicked.
         */
        protected void fireMarkClickedEvent(int index, double x, double y, boolean isShiftDown) {
            for (int i=0; i < listeners.size(); i++) {
                final Listener listener = listeners.elementAt(i);
                listener.markClicked(index, x, y, isShiftDown);
            }
        }

        protected void fireSelectionChanged() {
            for (int i = 0; i < listeners.size(); i++) {
                final Listener listener = listeners.elementAt(i);
                listener.selectionChanged(selectedPoints);
            }
        }


    }
}

