/*
 * DensityPlot.java
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

package dr.gui.chart;

import dr.util.Variate;

import java.awt.*;
import java.awt.geom.GeneralPath;

public class DensityPlot extends FrequencyPlot {

	boolean relativeDensity = true;
	int minimumBinCount;
	
	public DensityPlot(Variate data, int minimumBinCount) {
		super(data, minimumBinCount);
		this.minimumBinCount = minimumBinCount;
	}
		
	public DensityPlot(double[] data, int minimumBinCount) {
		super(data, minimumBinCount);
		this.minimumBinCount = minimumBinCount;
	}

	public void setRelativeDensity(boolean relative) {
		relativeDensity = relative;
		setData(raw, minimumBinCount);
	}

   /**	
	*	Set data
	*/
	public void setData(Variate data, int minimumBinCount) {
	
		this.raw = data;
		frequency = getFrequencyDistribution(data, minimumBinCount);
		
		Variate.Double xData = new Variate.Double();
		Variate.Double yData = new Variate.Double();

		double x = frequency.getLowerBound();
		double maxDensity = 0.0;
		for (int i = 0; i < frequency.getBinCount(); i++) {
			double density = frequency.getFrequency(i)/frequency.getBinSize()/data.getCount();
			if (density > maxDensity) maxDensity = density;
		}
		
		for (int i = 0; i < frequency.getBinCount(); i++) {
			xData.add(x + (frequency.getBinSize()/2.0));
			double density = frequency.getFrequency(i)/frequency.getBinSize()/data.getCount();
			if (relativeDensity) {
				yData.add(density/maxDensity);
			} else {
				yData.add(density);
			}
			x += frequency.getBinSize();
		}
		setData(xData, yData);
	}	
	
	/**	
	 *	Set bar fill style. Use a barPaint of null to not fill bar. 
	 *  Bar outline style is set using setLineStyle
	 */
	public void setBarFillStyle(Paint barPaint) {
		throw new IllegalArgumentException();
	}

	/**	
	 *	Paint data series
	 */
	protected void paintData(Graphics2D g2, Variate xData, Variate yData) {

		int n = xData.getCount();
					
		g2.setStroke(lineStroke);
		g2.setPaint(linePaint);
		
		
		float x = (float)transformX(xData.get(0));
		float y = (float)transformY(yData.get(0));
		
		GeneralPath path = new GeneralPath();
		path. moveTo(x, y);

		for (int i = 1; i < n; i++) {
			x = (float)transformX(xData.get(i));
			y = (float)transformY(yData.get(i));
			
			path.lineTo(x, y);
		}
		
		g2.draw(path);
	}	
}
