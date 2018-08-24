/*
 * CustomAxis.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 * @author Guy Baele
 */

public class CustomAxis extends Axis.AbstractAxis {

    /**
     * Empty constructor
     */
    public CustomAxis() {
    }

    /**
     * Axis flag constructor
     */
    public CustomAxis(int minAxisFlag, int maxAxisFlag) {
        setAxisFlags(minAxisFlag, maxAxisFlag);
    }

    /**
     * Transform a value
     */
    public double transform(double value) {
        return value;    // a linear transform !
    }

    /**
     * Untransform a value
     */
    public double untransform(double value) {
        return value;    // a linear transform !
	}

    /**
     *	Set the range of the data
     */
    public void setRange(double minValue, double maxValue) {
        if (!Double.isNaN(minValue)) {
            this.minData = minValue;
            this.minAxisFlag = (int)minValue;
        }
        if (!Double.isNaN(maxValue)) {
            this.maxData = maxValue;
            this.maxAxisFlag = (int)maxValue;
        }

        isCalibrated = false;
    }

    public void calibrate() {

        //do nothing
        majorTick=1;
        minorTick=1;

        minTick=minAxisFlag;
        maxTick=maxAxisFlag;

        majorTickCount = (int)((maxTick-minTick+1)/majorTick);
        minorTickCount = 0;

        minAxis=minTick-0.5;
        maxAxis=maxTick+0.5;

    }

}

