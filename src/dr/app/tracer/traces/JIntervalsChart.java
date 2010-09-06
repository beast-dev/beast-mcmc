package dr.app.tracer.traces;

import dr.app.gui.chart.Axis;
import dr.app.gui.chart.DiscreteAxis;
import dr.app.gui.chart.JChart;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.ArrayList;

public class JIntervalsChart extends JChart {

    //private static final int MAX_POINTS = 1000;

    // private boolean isLinePlot = false;

    private class Interval {
        String name;
        double value, upper, lower;
        boolean bold;

        Interval(String name, double value, double upper, double lower, boolean bold) {

            this.name = name;
            this.value = value;
            this.upper = upper;
            this.lower = lower;
            this.bold = bold;
        }
    }

    private final ArrayList<Interval> intervals = new ArrayList<Interval>();

    public JIntervalsChart(Axis yAxis) {
        super(new DiscreteAxis(true, true), yAxis);
    }

    public JIntervalsChart(Axis xAxis, Axis yAxis) {
        super(xAxis, yAxis);
        intervals.clear();
    }

    public void addIntervals(String name, double mean, double upper, double lower, boolean bold) {

        intervals.add(new Interval(name, mean, upper, lower, bold));

        xAxis.addRange(1, intervals.size());
        yAxis.addRange(lower, upper);

        recalibrate();
        repaint();
    }

    public void removeAllIntervals() {
        intervals.clear();
        xAxis.setRange(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        yAxis.setRange(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        removeAllPlots();
        recalibrate();
        repaint();
    }

    protected void calibrate(Graphics2D g2, Dimension size) {
    }

    protected boolean hasContents() {
        if (intervals.size() > 0) {
            return intervals.size() > 0;
        } else {
            return super.hasContents();
        }
    }

    protected void paintMajorTick(Graphics2D g2, double value, boolean horizontalAxis) {
        if (intervals.size() > 0 && horizontalAxis) {
            g2.setPaint(getAxisPaint());
            g2.setStroke(getAxisStroke());

            int index = ((int) value) - 1;
            Interval interval = intervals.get(index);
            String label = interval.name;

            double pos = transformX(value);

            Line2D line = new Line2D.Double(pos, getPlotBounds().getMaxY(), pos, getPlotBounds().getMaxY() + getMajorTickSize());
            g2.draw(line);

            g2.setPaint(getLabelPaint());
            double width = g2.getFontMetrics().stringWidth(label);
            g2.drawString(label, (float) (pos - (width / 2)), (float) (getPlotBounds().getMaxY() + (getMajorTickSize() * 1.25) + getXTickLabelOffset()));
        } else {
            super.paintMajorTick(g2, value, horizontalAxis);
        }
    }

    protected void paintContents(Graphics2D g2) {
        if (intervals.size() > 0) {
            for (int i = 0; i < intervals.size(); i++) {

                Interval interval = intervals.get(i);

                float x = (float) transformX(i + 1);
                float xLeft = (float) transformX(((double) i + 1) - 0.1);
                float xRight = (float) transformX(((double) i + 1) + 0.1);
                //float y = (float)transformY(interval.value);
                float yUpper = (float) transformY(interval.upper);
                float yLower = (float) transformY(interval.lower);

                GeneralPath path = new GeneralPath();
                path.moveTo(xLeft, yUpper);
                path.lineTo(xRight, yUpper);
                path.moveTo(x, yUpper);
                path.lineTo(x, yLower);
                path.moveTo(xLeft, yLower);
                path.lineTo(xRight, yLower);

                if (interval.bold) {
                    g2.setStroke(new BasicStroke(2.0f));
                } else {
                    g2.setStroke(new BasicStroke(1.0f));
                }
                g2.setPaint(Color.black);
                g2.draw(path);
            }
        } else {
            super.paintContents(g2);
        }

    }

}