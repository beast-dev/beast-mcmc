/*
 * JIntervalsChart.java
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
import dr.gui.chart.DiscreteAxis;
import dr.gui.chart.JChart;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.ArrayList;

public class JIntervalsChart extends JChart {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 8940887369520091972L;

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

	private ArrayList intervals = new ArrayList();
	
	public JIntervalsChart(Axis yAxis) {
		super(new DiscreteAxis(true, true), yAxis);
	}
	
	public void addIntervals(String name, double value, double upper, double lower, boolean bold) {
	
		intervals.add(new Interval(name, value, upper, lower, bold));
		
		xAxis.addRange(1, intervals.size());
		yAxis.addRange(lower, upper);

		recalibrate();
		repaint();	
	}
	
	public void removeAllIntervals() {
		intervals.clear();
		xAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
		yAxis.setRange(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY);
		removeAllPlots();
		recalibrate();
		repaint();
	}

	protected void calibrate(Graphics2D g2, Dimension size) { }
	
	protected boolean hasContents() { 
		return intervals.size() > 0;
	}

	protected void paintMajorTick(Graphics2D g2, double value, boolean horizontalAxis)
	{
		if (horizontalAxis) {
			g2.setPaint(getAxisPaint());
			g2.setStroke(getAxisStroke());
			
			int index = ((int)value) - 1;
			Interval interval = (Interval)intervals.get(index);
			String label = interval.name;
			
			double pos = transformX(value);
			
			Line2D line = new Line2D.Double(pos, getPlotBounds().getMaxY(), pos, getPlotBounds().getMaxY() + getMajorTickSize());
			g2.draw(line);
			
			g2.setPaint(getLabelPaint());
			double width = g2.getFontMetrics().stringWidth(label);
			g2.drawString(label, (float)(pos - (width / 2)), (float)(getPlotBounds().getMaxY() + (getMajorTickSize() * 1.25) + getXTickLabelOffset()));
		} else {
			super.paintMajorTick(g2, value, horizontalAxis);
		}
	}
		
	protected void paintContents(Graphics2D g2) {
		
		for (int i = 0; i < intervals.size(); i++) {
					
			Interval interval = (Interval)intervals.get(i);

			float x = (float)transformX(i + 1);
			float xLeft = (float)transformX(((double)i + 1) - 0.1);
			float xRight = (float)transformX(((double)i + 1) + 0.1);
			float yUpper = (float)transformY(interval.upper);
			float yLower = (float)transformY(interval.lower);
			
			GeneralPath path = new GeneralPath();
			path.moveTo(xLeft , yUpper);
			path.lineTo(xRight, yUpper);
			path.moveTo(x, yUpper);
			path.lineTo(x, yLower);
			path.moveTo(xLeft , yLower);
			path.lineTo(xRight, yLower);
				
			if (interval.bold) {
				g2.setStroke(new BasicStroke(2.0f));
			} else {
				g2.setStroke(new BasicStroke(1.0f));
			}
			g2.setPaint(Color.black);			
			g2.draw(path);
		}		
		
	}

}