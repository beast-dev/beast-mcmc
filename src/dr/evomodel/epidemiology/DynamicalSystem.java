/*
 * DynamicalSystem.java
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
import java.util.HashMap;

/**
 * @author Trevor Bedford
 */
public class DynamicalSystem {

    private List<DynamicalVariable> variables = new ArrayList<DynamicalVariable>();
    private List<DynamicalForce> forces = new ArrayList<DynamicalForce>();
    private HashMap<String,DynamicalVariable> varMap = new HashMap<String,DynamicalVariable>();
    private HashMap<String,DynamicalForce> forceMap = new HashMap<String,DynamicalForce>();

    private double currentTime = 0.0;
    private double startTime = 0.0;
    private double timeStep = 0.0;

    private double storedCurrentTime = 0.0;
    private double storedStartTime = 0.0;

    public static void main(String[] args) {

        DynamicalSystem syst = new DynamicalSystem(0, 0.001);

        double transmissionRate = 0.027;
        double recoveryRate = 0.00054;

        syst.addVariable("susceptibles", 330);
        syst.addVariable("infecteds", 23325);
        syst.addVariable("recovereds", 4524);
        syst.addVariable("total", 330 + 23325 + 4524);
        syst.addForce("contact", transmissionRate, new String[]{"infecteds","susceptibles"}, new String[]{"total"}, "susceptibles", "infecteds");
        syst.addForce("recovery", recoveryRate, new String[]{"infecteds"}, new String[]{}, "infecteds", "recovereds");

//        while (syst.getTime() < 400) {
//            syst.step();
//        }
//        syst.print(0,400,1);

        double val = syst.getValue("susceptibles", 500);
        System.out.println(val);
        syst.print(0,500,10);

    }

    public DynamicalSystem(double t0, double dt) {
        startTime = t0;
        currentTime = startTime;
        timeStep = dt;
    }

    public double getTime() {
        return currentTime;
    }

    public int size() {
        return variables.size();
    }

    public DynamicalVariable getVar(String n) {
        return varMap.get(n);
    }

   public DynamicalForce getForce(String n) {
        return forceMap.get(n);
    }

    public void resetVar(String n, double v0) {
        DynamicalVariable var = getVar(n);
        var.reset(startTime, v0);
    }

    public void resetForce(String n, double c) {
        DynamicalForce frc = getForce(n);
        frc.reset(c);
    }

    public void resetTime() {
        currentTime = startTime;
    }

    // copy values to stored state
    public void store() {
        storedCurrentTime = currentTime;
        storedStartTime = startTime;
        for (DynamicalVariable var : variables) {
            var.store();
        }
        for (DynamicalForce frc : forces) {
            frc.store();
        }
    }

    // copy values from stored state
    public void restore() {
        currentTime = storedCurrentTime;
        startTime = storedStartTime;
        for (DynamicalVariable var : variables) {
            var.restore();
        }
        for (DynamicalForce frc : forces) {
            frc.restore();
        }
    }

    // get value of indexed variable at time t
    // dynamically extend trace
    public double getValue(int index, double t) {
        while (currentTime < t) {
            step();
        }
        DynamicalVariable var = variables.get(index);
        return var.getValue(t);
    }

    // get value of named variable at time t
    // dynamically extend trace
    public double getValue(String n, double t) {
        while (currentTime < t) {
            step();
        }
        DynamicalVariable var = getVar(n);
        return var.getValue(t);
    }

    // get integral of indexed variable between times start and finish
    // dynamically extend trace
    public double getIntegral(int index, double start, double finish) {
        while (currentTime < finish) {
            step();
        }
        DynamicalVariable var = variables.get(index);
        return var.getIntegral(start, finish);
    }

    // get integral of named variable between times start and finish
    // dynamically extend trace
    public double getIntegral(String n, double start, double finish) {
        while (currentTime < finish) {
            step();
        }
        DynamicalVariable var = getVar(n);
        return var.getIntegral(start, finish);
    }

   // get average of indexed variable between times start and finish
    // dynamically extend trace
    public double getAverage(int index, double start, double finish) {
        while (currentTime < finish) {
            step();
        }
        DynamicalVariable var = variables.get(index);
        return var.getAverage(start, finish);
    }

    // get average of named variable between times start and finish
    // dynamically extend trace
    public double getAverage(String n, double start, double finish) {
        while (currentTime < finish) {
            step();
        }
        DynamicalVariable var = getVar(n);
        return var.getAverage(start, finish);
    }

    public void addVariable(String n, double v0) {
        DynamicalVariable var = new DynamicalVariable(n, startTime, v0);
        varMap.put(n, var);
        variables.add(var);
    }

    public void addForce(String n, double coeff, String[] mult, String[] div, String increasing, String decreasing) {
        DynamicalForce frc = new DynamicalForce(n, coeff, getVar(increasing), getVar(decreasing));
        for (String s : mult) {
            frc.addMultiplier(getVar(s));
        }
        for (String s : div) {
            frc.addDivisor(getVar(s));
        }
        forceMap.put(n, frc);
        forces.add(frc);
    }

    public void step() {
        for (DynamicalForce frc : forces) {
            frc.modCurrentValue(currentTime, timeStep);
        }
        for (DynamicalVariable var : variables) {
            var.modCurrentTime(timeStep);
            var.pushCurrentState();
        }
        currentTime += timeStep;
    }

    public void print(double start, double finish, double step) {

        System.out.print("time");
        for (DynamicalVariable var : variables) {
            System.out.print("\t" + var.getName());
        }
        System.out.println();

        for (double t=start; t<=finish; t += step) {
            System.out.printf("%.3f", t);
            for (DynamicalVariable var : variables) {
                System.out.printf("\t%.3f", var.getValue(t));
            }
            System.out.println();
        }

    }

    public String printValues(String n, double start, double finish, double step) {

        String out = "";
        DynamicalVariable var = getVar(n);
        for (double t=start; t<=finish; t += step) {
            double v = var.getValue(t);
            out += "\t";
            out += Double.toString(v);
        }
        return out;

    }

}