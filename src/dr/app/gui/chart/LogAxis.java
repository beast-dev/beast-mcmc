/*
 * LogAxis.java
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

public class LogAxis extends Axis.AbstractAxis {

//	private int extraMinorTickCount;


	/**
	 * Empty constructor
	 */
	public LogAxis() {
		super();
        setSignficantFigures(0);
    }

	/**
	 * Axis flag constructor
	 */
	public LogAxis(int minAxisFlag, int maxAxisFlag) {

		setAxisFlags(minAxisFlag, maxAxisFlag);
	}

	/**
	*	Set axis flags
	*/
	public void setAxisFlags(int minAxisFlag, int maxAxisFlag) {

		if (minAxisFlag==AT_DATA || minAxisFlag==AT_ZERO) {
			minAxisFlag = this.minAxisFlag;
		}

		if (maxAxisFlag==AT_DATA || maxAxisFlag==AT_ZERO) {
			maxAxisFlag = this.maxAxisFlag;
		}

		super.setAxisFlags(minAxisFlag, maxAxisFlag);
	}

	/**
	*	Set the range of the data.
	*/
	public void setRange(double minValue, double maxValue) {

		if (maxValue <= 0.0)
			maxValue = 1.0E-100;

		if (minValue <= 0.0)
			minValue = maxValue;
		super.setRange(minValue, maxValue);

	}

	/**
	*	Adds the range of the data.
	*/
	public void addRange(double minValue, double maxValue) {

		if (maxValue <= 0.0)
			maxValue = 1.0E-100;

		if (minValue <= 0.0)
			minValue = maxValue;
 		super.addRange(minValue, maxValue);

	}

	/**
	* Calculate the optimum minimum tick. For a log scale this is always a power of 10
	*/
	public void calcMinTick() {
		minTick=1;
		// Find the nearest multiple of majorTick below minData
		if (minData>1) { // work upwards
			while ((minTick*10)<minData)
				minTick*=10;
		} else if (minData<1) { // work downwards
			while (minTick>minData)
				minTick/=10;
		}
	}

	/**
	* Calculate the optimum maximum tick. For a log scale this is always a power of 10
	*/
	public void calcMaxTick() {
		maxTick=1;
		// Find the nearest multiple of majorTick above maxData
		if (maxData>1) { // work upwards
			while (maxTick<maxData)
				maxTick*=10;
		} else if (maxData<1) { // work downwards
			while ((maxTick/10)>maxData)
				maxTick/=10;
		}
	}

	/**
	* Calculate the optimum major tick distance.
	*/
	public void calcMajorTick() {
		majorTick=10;
		majorTickCount=(int)Math.round(log10(maxTick/minTick))+1;
	}

	/**
	* Calculate the optimum minor tick distance.
	*/
	public void calcMinorTick() {
		minorTick=1;
		minorTickCount=8; // 1 and 10 are major ticks
	}

	/**
	* Handles axis flags.
	*/
	public void handleAxisFlags() {
		if (minAxisFlag==AT_MINOR_TICK) {
			if (minAxis+minTick<=minData) {
				while ((minAxis+minTick)<=minData) {
					minAxis+=minTick;
				}
				majorTickCount--;
				minTick*=10;
			}
		}

		if (maxAxisFlag==AT_MINOR_TICK) {
			if (maxAxis-(maxTick/10)>=maxData) {
				majorTickCount--;
				maxTick/=10;
				while ((maxAxis-maxTick)>=maxData) {
					maxAxis-=maxTick;
				}
			}
		}
	}

	/**
	*	Transform a value onto the log axis
	*/
	public double transform(double value) {
		return log10(value);
	}

	/**
	*	Untransform a value from the log axis
	*/
	public double untransform(double value) {
		return Math.pow(10.0, value);
	}

	/**
	*	Returns the number of minor tick marks within each major one
	*	The last major tick may have less minors that the other.
	*/
	public int getMinorTickCount(int majorTickIndex) {
		if (!isCalibrated)
			calibrate();

		if (majorTickIndex==majorTickCount-1)
			return (int)((maxAxis-maxTick)/maxTick);
		else if (majorTickIndex==-1) {
			return (int)((minTick-minAxis)/(minTick/10));
		} else
			return minorTickCount;
	}


	/**
	*	Returns the value of the majorTickIndex'th major tick
	*/
	public double getMajorTickValue(int majorTickIndex) {
		if (!isCalibrated)
			calibrate();

		if (majorTickIndex==majorTickCount-1)
			return maxTick;
		else
			return Math.pow(10, majorTickIndex)*minTick;
	}

	/**
	*	Returns the value of the minorTickIndex'th minor tick
	*/
	public double getMinorTickValue(int minorTickIndex, int majorTickIndex) {
		if (!isCalibrated)
			calibrate();

		if (majorTickIndex==-1)
			return minTick-((minorTickIndex+1)*(minTick/10));
		else
			return (minorTickIndex+2)*getMajorTickValue(majorTickIndex);

		//	This last line is the equivalent of:
		//	return ((minorTickNo+1)*getMajorTick(majorTickNo))+getMajorTick(majorTickNo);
	}
}

