/*
 * DynamicalVariable.java
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

package dr.evomodel.epidemiology;

import java.util.ArrayList;
import java.util.List;
import java.util.Collections;

/**
 * @author Trevor Bedford
 */
public class DynamicalVariable {

    private String name = "variable";
    private List<Double> times = new ArrayList<Double>();
    private List<Double> values = new ArrayList<Double>();
    private double currentTime = 0.0;
    private double currentValue = 0.0;

    private List<Double> storedTimes = new ArrayList<Double>();
    private List<Double> storedValues = new ArrayList<Double>();
    private double storedTime = 0.0;
    private double storedValue = 0.0;

    // initialize variable
    public DynamicalVariable(String n, double t0, double v0) {
        name = n;
        currentTime = t0;
        currentValue = v0;
        add(currentTime, currentValue);
    }

    // add a time / value pair
    public void add(double t, double v) {
        Double tD = new Double(t);
        Double vD = new Double(v);
        times.add(tD);
        values.add(vD);
    }

    public void pushCurrentState() {
        add(currentTime,currentValue);
    }

    public void modCurrentTime(double dt) {
        currentTime += dt;
    }

    public void modCurrentValue(double dv) {
        currentValue += dv;
    }

    // reset arrays
    public void reset(double t0, double v0) {
        currentTime = t0;
        currentValue = v0;
        times.clear();
        values.clear();
        add(currentTime, currentValue);
    }

    // copy values to stored state
    public void store() {

        storedTimes.clear();
        for (Double t : times) {
            Double tD = new Double(t.doubleValue());
            storedTimes.add(tD);
        }
        storedValues.clear();
        for (Double v : values) {
            Double vD = new Double(v.doubleValue());
            storedValues.add(vD);
        }
        storedTime = currentTime;
        storedValue = currentValue;

    }

    // copy values from stored state
    public void restore() {

        times.clear();
        for (Double t : storedTimes) {
            Double tD = new Double(t.doubleValue());
            times.add(tD);
        }
        values.clear();
        for (Double v : storedValues) {
            Double vD = new Double(v.doubleValue());
            values.add(vD);
        }
        currentTime = storedTime;
        currentValue = storedValue;

    }

    public int size() {
        return times.size();
    }

    public String getName() {
        return name;
    }

    public double getTimeAtIndex(int i) {
        return times.get(i).doubleValue();
    }

    public double getValueAtIndex(int i) {
        return values.get(i).doubleValue();
    }

    // pull value closest to given point in time
    // uses linear interpolation when the exact time point isn't present
    public double getValue(double t) {

        double rval = 0.0;
        Double tD = new Double(t);
        int index = Collections.binarySearch(times, tD);

        // found exact match
        if (index >= 0) {
            Double vD = values.get(index);
            rval = vD.doubleValue();
        }

        // need to interpolate between two adjacent values
        else {
            int a = Math.abs(index) - 2;
            int b = Math.abs(index) - 1;
            double t0 = getTimeAtIndex(a);
            double v0 = getValueAtIndex(a);
            double t1 = getTimeAtIndex(b);
            double v1 = getValueAtIndex(b);
            rval = v0 + (t-t0) * ( (v1-v0) / (t1-t0) );
        }

        return rval;
    }

    // return integral between two points in time
    // approximation via trapezoidal rule
    public double getIntegral(double start, double finish) {

        // first index after start and last index before finish
        int a = 0;
        int b = 0;
        double sum = 0.0;

        // if not exact match, need first index after than start
        int index = Collections.binarySearch(times, new Double(start));
        if (index >= 0)
            a = index;
        else
            a = Math.abs(index) - 1;

        // if not exact match, need last index before than finish
        index = Collections.binarySearch(times, new Double(finish));
        if (index >= 0)
            b = index;
        else
            b = Math.abs(index) - 2;

        // from start to a
        sum += 0.5 * (getTimeAtIndex(a) - start) * (getValueAtIndex(a) + getValue(start));

        // between a and b
        for (index = a; index < b; index++) {
            sum += 0.5 * (getTimeAtIndex(index+1) - getTimeAtIndex(index)) * (getValueAtIndex(index+1) + getValueAtIndex(index));
        }

        // b to finish
        sum += 0.5 * (finish - getTimeAtIndex(b)) * (getValue(finish) + getValueAtIndex(b));

        return sum;

    }

    public double getAverage(double start, double finish) {
        return 0.5*(getValue(start) + getValue(finish));
    }

    public void print() {
        System.out.println(name);
        for (int i=0; i<size(); i++) {
            System.out.println(times.get(i) + "\t" + values.get(i));
        }
    }

}