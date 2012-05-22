/*
 * DiscreteAxis.java
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

package dr.app.gui.chart;



public class DiscreteAxis extends Axis.AbstractAxis {

	private boolean originBetweenCategories;
	private boolean showEveryCategory;

	public DiscreteAxis(boolean originBetweenCategories, boolean showEveryCategory) {
		super(AT_MAJOR_TICK, AT_MAJOR_TICK, true);

		this.originBetweenCategories = originBetweenCategories;
		this.showEveryCategory = showEveryCategory;
	}

	/**
	*	Transform a value
	*/
	public double transform(double value) {
		return value;	// a linear transform !
	}

	/**
	*	Untransform a value
	*/
	public double untransform(double value) {
		return value;	// a linear transform !
	}

	public void calibrate() {
		majorTick=1;
		minorTick=1;

		minTick=minData;
		maxTick=maxData;
		majorTickCount = (int)((maxTick-minTick)/majorTick)+1;
		minorTickCount = 0;

		if (!showEveryCategory) {
			while (majorTickCount > prefMajorTickCount) {
				majorTickCount=(int)((maxTick-minTick)/(majorTick*2))+1;
				if (majorTickCount > prefMajorTickCount) {
					majorTick*=2;
					break;
				}
				majorTickCount=(int)((maxTick-minTick)/(majorTick*4))+1;
				if (majorTickCount > prefMajorTickCount) {
					majorTick*=4;
					break;
				}
				majorTickCount=(int)((maxTick-minTick)/(majorTick*5))+1;
				if (majorTickCount > prefMajorTickCount) {
					majorTick*=5;
					break;
				}

				majorTick*=10;
				majorTickCount=(int)((maxTick-minTick)/majorTick)+1;
			}
		}

		minorTickCount=(int)(majorTick)-1;

		minAxis=minTick;
		maxAxis=maxTick;

		handleAxisFlags();

		if (originBetweenCategories) {
			minAxis-=0.5;
			maxAxis+=0.5;
		}
	}
}
