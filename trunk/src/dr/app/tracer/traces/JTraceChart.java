/*
 * JTraceChart.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tracer.traces;

import dr.app.gui.chart.Axis;
import dr.app.gui.chart.DiscreteJChart;
import dr.app.gui.chart.Plot;
import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;
import java.util.List;

public class JTraceChart extends DiscreteJChart {

    private static final int SAMPLE_POINTS = 1000;

    private boolean useSample = false;
    private boolean isLinePlot = true;
    private static final int BURNIN_TRANSLUCENCY = 72;

    private class Trace {
        int stateStart;
        int stateStep;

        int sampleCount;
        int sampleStep;

        List<Double> states = new ArrayList<Double>();
        List<Double> values;

        Trace(int stateStart, int stateStep, List<Double> values) {

            this.stateStart = stateStart;
            this.stateStep = stateStep;

            this.values = values;

            sampleCount = values.size();
            sampleStep = 1;
            while (sampleCount > SAMPLE_POINTS) {
                sampleStep *= 2;
                sampleCount = values.size() / sampleStep;
            }

            double ix = stateStart;
            for (int j = 0; j < values.size(); j++) {
                this.states.add(j, ix);
                ix += stateStep;

            }
        }
    }

    private final ArrayList<Trace> traces = new ArrayList<Trace>();
    private final ArrayList<Trace> burninTraces = new ArrayList<Trace>();

    public JTraceChart(Axis xAxis, Axis yAxis) {
        super(xAxis, yAxis);
    }

    public void setUseSample(boolean useSample) {
        this.useSample = useSample;
    }

    public void setIsLinePlot(boolean isLinePlot) {
        this.isLinePlot = isLinePlot;
    }

    public void addTrace(String name, int stateStart, int stateStep, List<Double> values, List<Double> burninValues, Paint paint) {

        Variate.D yd = new Variate.D(values);

        xAxis.addRange(0, stateStart + (values.size() * stateStep) - stateStep);
        yAxis.addRange(yd.getMin(), yd.getMax());

        traces.add(new Trace(stateStart, stateStep, values));
        if (burninValues != null) {
            burninTraces.add(new Trace(0, stateStep, burninValues));
        }

        Plot plot = new Plot.AbstractPlot() { // create a dummy plot to store paint styles
            protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {
            }
        };

        plot.setLineColor(paint);
        plot.setName(name);
        addPlot(plot);

        recalibrate();
        repaint();
    }

    public List<Double> getTraceStates(int index) {
        Trace trace = traces.get(index);
        return trace.states;
    }

    public List<Double> getTraceValues(int index) {
        Trace trace = traces.get(index);
        return trace.values;
    }

    public void removeAllTraces() {
        traces.clear();
        burninTraces.clear();
        xAxis.setRange(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        yAxis.setRange(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        removeAllPlots();
        recalibrate();
        repaint();
    }

    protected void calibrate(Graphics2D g2, Dimension size) {
    }

    protected boolean hasContents() {
        return traces.size() > 0;
    }

    protected void paintContents(Graphics2D g2) {

        if (!isLinePlot) {
            g2.setStroke(new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER));
        } else {
            g2.setStroke(new BasicStroke(0.5f));
        }


        for (int i = 0; i < traces.size(); i++) {
            Trace trace = traces.get(i);

            int sampleFrequency = 1;

            if (useSample) {
                sampleFrequency = trace.states.size() / trace.sampleCount;
            }

            paintTrace(g2, trace, getPlot(i).getLineColor(), sampleFrequency);

            if (i < burninTraces.size()) {
                Trace burninTrace = burninTraces.get(i);
                Color colour = (Color)getPlot(i).getLineColor();
                Color burninColor = new Color(colour.getRed(), colour.getGreen(), colour.getBlue(), BURNIN_TRANSLUCENCY);
                paintTrace(g2, burninTrace, burninColor, sampleFrequency);
            }
        }

    }

    private void paintTrace(Graphics2D g2, Trace trace, Paint paint, int sampleFrequency) {
        float x = (float) transformX(trace.states.get(0));
        float y = (float) transformY(trace.values.get(0));

        GeneralPath path = new GeneralPath();
        path.moveTo(x, y);
        if (!isLinePlot) {
            path.lineTo(x, y);
        }

        int k = sampleFrequency;

        for (int j = 1; j < trace.states.size(); j++) {

            x = (float) transformX(trace.states.get(k));
            y = (float) transformY(trace.values.get(k));

            if (!isLinePlot) {
                path.moveTo(x, y);
            }
            path.lineTo(x, y);

            k += sampleFrequency;
            if (k >= trace.states.size()) {
                break;
            }
        }

        g2.setPaint(paint);
        g2.draw(path);
    }

}