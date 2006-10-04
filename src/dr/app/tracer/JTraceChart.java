/*
 * JTraceChart.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.tracer;

import dr.gui.chart.Axis;
import dr.gui.chart.JChart;
import dr.gui.chart.Plot;
import dr.util.Variate;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;

public class JTraceChart extends JChart {

	private static final int SAMPLE_POINTS = 1000;
	
	private boolean useSample = true;
	private boolean isLinePlot = true;
	
	private class Trace {
		int stateStart;
		int stateStep;

		int sampleCount;
		int sampleStep;

		double[] states;
		double[] values;
		
		Trace(int stateStart, int stateStep, double[] values) {
		
			this.stateStart = stateStart;
			this.stateStep = stateStep;
			
			this.values = values;
			
			sampleCount = values.length;
			sampleStep = 1;
			while (sampleCount > SAMPLE_POINTS) {
				sampleStep *= 2;
				sampleCount = values.length / sampleStep;
			}

			this.states = new double[values.length];
			
			int k = 0;
			double ix = stateStart;
			for (int j = 0; j < values.length; j++) {
				this.states[j] = ix;
				ix += stateStep;
				
			}
		}
	}

	private ArrayList traces = new ArrayList();
	
	public JTraceChart(Axis xAxis, Axis yAxis) {
		super(xAxis, yAxis);
	}
	
	public void setUseSample(boolean useSample) {
		this.useSample = useSample;
	}
	
	public void setIsLinePlot(boolean isLinePlot) {
		this.isLinePlot = isLinePlot;
	}
	
	public void addTrace(String name, int stateStart, int stateStep, double[] values) {
	
		dr.util.Variate.Double yd = new dr.util.Variate.Double(values);
		
		xAxis.addRange(0, stateStart + (values.length * stateStep) - stateStep);
		yAxis.addRange(yd.getMin(), yd.getMax());
		
		traces.add(new Trace(stateStart, stateStep, values));
		
		Plot plot = new Plot.AbstractPlot() { // create a dummy plot to store paint styles
			protected void paintData(Graphics2D g2, dr.util.Variate xData, dr.util.Variate yData) { }
		};
		Paint[] paints = new Paint[] {
			Color.black,
			Color.blue,
			Color.red,
			Color.green,
			Color.magenta,
			Color.gray,
			Color.pink,
			Color.orange,
			Color.lightGray,
			Color.darkGray
		};
		int index = getPlotCount() % paints.length;
		plot.setLineColor(paints[index]);
		plot.setName(name);
		addPlot(plot);
		
		recalibrate();
		repaint();	
	}

    public double[] getTraceStates(int index) {
        Trace trace = (Trace)traces.get(index);
        return trace.states;
    }

    public double[] getTraceValues(int index) {
        Trace trace = (Trace)traces.get(index);
        return trace.values;
	}
	
	public void removeAllTraces() {
		traces.clear();
		xAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
		yAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
		removeAllPlots();
		recalibrate();
		repaint();
	}
 
	protected void calibrate(Graphics2D g2, Dimension size) { }
	
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
					
			Trace trace = (Trace)traces.get(i);

			float x = (float)transformX(trace.states[0]);
			float y = (float)transformY(trace.values[0]);
			
			GeneralPath path = new GeneralPath();
			path.moveTo(x, y);
			if (!isLinePlot) {
				path.lineTo(x, y);
			}
			
			int n = trace.states.length;
			int ik = 1;
			
			if (useSample) {
				n = trace.sampleCount;
				ik = trace.states.length / n;
			}
		
			int k = ik;

			for (int j = 1; j < n; j++) {
			
				x = (float)transformX(trace.states[k]);
				y = (float)transformY(trace.values[k]);
				
				if (!isLinePlot) {
					path.moveTo(x, y);
				}
				path.lineTo(x, y);
				
				k += ik;
			}

			g2.setPaint(getPlot(i).getLineColor());			
			g2.draw(path);
		}		
		
	}

}