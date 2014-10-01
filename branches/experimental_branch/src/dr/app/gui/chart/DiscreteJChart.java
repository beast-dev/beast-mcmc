package dr.app.gui.chart;

import java.awt.*;
import java.awt.geom.Line2D;
import java.util.Map;

/**
 *
 */
public class DiscreteJChart extends JChart {

    private Map<Integer, String> categoryDataMap;

    public DiscreteJChart(Axis xAxis, Axis yAxis) {
        super(xAxis, yAxis);
    }

    public void setXAxis(boolean isInteger, Map<Integer, String> categoryDataMap) {
        if (categoryDataMap.isEmpty()) {
            if (isInteger) {
                super.setXAxis(new DiscreteAxis(true, true));
            } else {
                super.setXAxis(new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK_PLUS));
            }

        } else {
            super.setXAxis(new DiscreteAxis(true, true));
        }

        this.categoryDataMap = categoryDataMap;
    }

    public void setYAxis(boolean isInteger, Map<Integer, String> categoryDataMap) {
        if (categoryDataMap.isEmpty()) {
            if (isInteger) {
                super.setYAxis(new DiscreteAxis(true, true));
            } else {
                super.setYAxis(new LinearAxis());
            }

        } else {
            super.setYAxis(new DiscreteAxis(true, true));
        }

        this.categoryDataMap = categoryDataMap;
    }

    protected void paintMajorTick(Graphics2D g2, double value, String label, boolean horizontalAxis) {
        g2.setPaint(getAxisPaint());
        g2.setStroke(getAxisStroke());

        if (label == null) label = " ";

        if (horizontalAxis) {
            double pos = transformX(value);

            Line2D line = new Line2D.Double(pos, getPlotBounds().getMaxY(), pos, getPlotBounds().getMaxY() + getMajorTickSize());
            g2.draw(line);

            g2.setPaint(getLabelPaint());
            double width = g2.getFontMetrics().stringWidth(label);

            if (label == null) label = xAxis.format(value);
            g2.drawString(label, (float) (pos - (width / 2)), (float) (getPlotBounds().getMaxY() + (getMajorTickSize() * 1.25) + getXTickLabelOffset()));

        } else {
            double pos = transformY(value);

            Line2D line = new Line2D.Double(getPlotBounds().getMinX(), pos, getPlotBounds().getMinX() - getMajorTickSize(), pos);
            g2.draw(line);

            g2.setPaint(getLabelPaint());
            double width = g2.getFontMetrics().stringWidth(label);

            if (label == null) label = yAxis.format(value);
			g2.drawString(label, (float)(getPlotBounds().getMinX() - width - (getMajorTickSize() * 1.25)), (float)(pos + getYTickLabelOffset()));
        }
    }

    protected void paintAxis(Graphics2D g2, Axis axis, boolean horizontalAxis) {
        if ((!categoryDataMap.isEmpty()) && axis.getIsDiscrete()) {
            int n1 = axis.getMajorTickCount();
            int n2, i, j;

            n2 = axis.getMinorTickCount(-1);

            for (i = 0; i < n1; i++) {
                paintMajorTick(g2, axis.getMajorTickValue(i), categoryDataMap.get((int) axis.getMajorTickValue(i)), horizontalAxis);
                n2 = axis.getMinorTickCount(i);

                if (i == (n1 - 1) && axis.getLabelLast()) { // Draw last minor tick as a major one

                    paintMajorTick(g2, axis.getMinorTickValue(0, i), categoryDataMap.get((int) axis.getMinorTickValue(0, i)), horizontalAxis);

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
            super.paintAxis(g2, axis, horizontalAxis);
        }

    }

}
