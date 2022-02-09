/*
 * LinearAxis.java
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


/**
 * Same behavior as LinearAxis,
 * except to remove the improper handling very small range.
 * Tracer issue #167
 * @author Walter Xie
 */
public class LinearAxis2 extends LinearAxis {

	/**
	 * Empty constructor
	 */
	public LinearAxis2() { }

	/**
	 * Axis flag constructor
	 */
	public LinearAxis2(int minAxisFlag, int maxAxisFlag) {
		super(minAxisFlag, maxAxisFlag);
	}


	/**
	 * Plotting small-range values (range < 1.0E-30) is messed up by
	 * the code in its parent method,
	 * but not sure why {@link Axis.AbstractAxis#calibrate()} needs it.
	 * So create a 2nd LinearAxis to overwrite this method.
	 */
	@Override
	public void calibrate() {
		double minValue = minData;
		double maxValue = maxData;

		if( Double.isInfinite(minValue) || Double.isNaN(minValue) ||
				Double.isInfinite(maxValue) || Double.isNaN(maxValue)) {
			// I am not sure which exception is appropriate here.
			throw new ChartRuntimeException("Illegal range values, can't calibrate");
		}

		if (minAxisFlag==AT_ZERO ) {
			minValue = 0;
		} else if (minAxisFlag == AT_VALUE) {
			minValue = this.minValue;
		}

		if (maxAxisFlag==AT_ZERO) {
			maxValue = 0;
		} else if (maxAxisFlag == AT_VALUE) {
			maxValue = this.maxValue;
		}

		double range = maxValue - minValue;
		if (range < 0.0) {
			range = 0.0;
		}

		setEpsilon(range * 1.0E-10);

		if (isAutomatic) {
			// We must find the optimum minMajorTick and maxMajorTick so
			// that they contain the data range (minData to maxData) and
			// are in the right order of magnitude

			if (range == 0) {

				minTick = -1.0;
				maxTick = 1.0;
				majorTick = 1.0;
				minorTick = majorTick;
				majorTickCount = 1;
				minorTickCount = 0;

			} else {
				// First find order of magnitude below the data range...
				majorTick = Math.pow(10.0, Math.floor(log10(range)));

				calcMinTick();
				calcMaxTick();

				calcMajorTick();
				calcMinorTick();
			}
		}

		minAxis = minTick;
		maxAxis = maxTick;

		handleAxisFlags();

		isCalibrated=true;
	}
}

