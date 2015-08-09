/*
 * DynamicalForce.java
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

/**
 * @author Trevor Bedford
 */
public class DynamicalForce {

    private String name = "force";
    private double coefficient = 0.0;
    private List<DynamicalVariable> multipliers = new ArrayList<DynamicalVariable>();
    private List<DynamicalVariable> divisors = new ArrayList<DynamicalVariable>();
    private DynamicalVariable increasingVariable;
    private DynamicalVariable decreasingVariable;

    private double storedCoefficient = 0.0;

    public DynamicalForce(String n, double c, DynamicalVariable increasing, DynamicalVariable decreasing) {
        name = n;
        coefficient = c;
        increasingVariable = increasing;
        decreasingVariable = decreasing;
    }

    public String getName() {
        return name;
    }

    public void addMultiplier(DynamicalVariable var) {
        multipliers.add(var);
    }

    public void addDivisor(DynamicalVariable var) {
        divisors.add(var);
    }

    // reset
    public void reset(double c) {
        coefficient = c;
    }

    public void store() {
        storedCoefficient = coefficient;
    }

    public void restore() {
        coefficient = storedCoefficient;
    }

    public double getForce(double t) {
        double force = coefficient;
        for (DynamicalVariable var : multipliers) {
            force *= var.getValue(t);
        }
        for (DynamicalVariable var : divisors) {
            force /= var.getValue(t);
        }
        if (increasingVariable.getValue(t) < 0 || decreasingVariable.getValue(t) < 0) {
            force = 0.0;
        }
        return force;
    }

    public void modCurrentValue(double t, double dt) {
        double dv = getForce(t) * dt;
        increasingVariable.modCurrentValue(dv);
        decreasingVariable.modCurrentValue(-dv);
    }

}