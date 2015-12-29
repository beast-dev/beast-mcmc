/*
 * CoalGenData.java
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

package dr.app.coalgen;

import dr.evolution.coalescent.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.inference.trace.LogFileTraces;

import java.io.File;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class CoalGenData {
    public static final String version = "1.0";

    public static String[] demographicModels = {"Constant Population",
            "Exponential Growth (Growth Rate)",
            "Exponential Growth (Doubling Time)",
            "Logistic Growth (Growth Rate)",
            "Logistic Growth (Doubling Time)",
            "Expansion (Growth Rate)",
            "Expansion (Doubling Time)",
            "Piecewise Constant (Skyline)",
            "Piecewise Linear (Skyline)",
            "Piecewise Linear (GMRF Skyride)"};

    public static String[][] argumentGuesses = {
            {"populationsize", "population", "popsize", "n0", "size", "pop"},
            {"ancestralsize", "ancestralproportion", "proportion", "ancestral", "n1"},
            {"growthrate", "growth", "rate", "r"},
            {"doublingtime", "doubling", "time", "t"},
            {"logisticshape", "halflife", "t50", "time50", "shape"},
            {"populationsize", "population", "popsize", "n0", "size", "pop"},
            {"groupsize"}};

    public static String[] argumentNames = new String[]{
            "Population Size", "Ancestral Proportion", "Growth Rate", "Doubling Time", "Logistic Shape (Half-life)",
            "Population Sizes", "Group Sizes"
    };

    public int[][] argumentIndices = {{0}, {0, 2}, {0, 3}, {0, 2, 4}, {0, 3, 4}, {0, 1, 2}, {0, 1, 3}, {5, 6}, {5, 6}, {5}};

    public File logFile = null;
    public File treesFile = null;
    public LogFileTraces traces = null;
    public TaxonList taxonList = new Taxa();
    public int demographicModel = 0;

    public int popSizeCount = 0;
    public int groupSizeCount = 0;

    public double[] argumentValues = new double[] {1000, 0.1, 0.5, 10.0, 50.0, 1000, 1};
    public boolean[] argumentFromTraces = new boolean[argumentNames.length];
    public int[] argumentTraces = new int[argumentNames.length];
    public int replicateCount = 1000;

    public CoalGenData() {
    }

    public void setDemographicModel(int model) {
        demographicModel = model;
        argumentTraces = new int[argumentIndices[demographicModel].length];
        for (int i = 0; i < argumentTraces.length; i++) {
            argumentTraces[i] = -1;
        }
    }

    private int current = 0;

    public boolean hasNext() {
        if (traces != null) {
            return (current < traces.getStateCount());
        }

        return (current < replicateCount);
    }

    private DemographicFunction demo = null;

    public DemographicFunction nextDemographic() {
        int n = replicateCount;

        if (traces != null) {
            n = traces.getStateCount();
        }

        if (current == n) {
            current = 0;
            demo = null;
            return null;
        }

        if (demo == null) {

            current = 0;

            if (traces != null) {
                traces.setBurnIn(0);
            }

            if (demographicModel == 0) { // Constant Size
                demo = new ConstantPopulation(Units.Type.YEARS);
            } else if (demographicModel == 1 || demographicModel == 2) { // Exponential Growth
                demo = new ExponentialGrowth(Units.Type.YEARS);
            } else if (demographicModel == 3 || demographicModel == 4) { // Logistic Growth
                demo = new LogisticGrowth(Units.Type.YEARS);
            } else if (demographicModel == 5 || demographicModel == 6) { // Expansion Growth
                demo = new Expansion(Units.Type.YEARS);
            } else if (demographicModel == 7) { // Piecewise Constant (skyline)
                demo = new PiecewiseConstantPopulation(Units.Type.YEARS);
            } else if (demographicModel == 8) { // Piecewise Linear (skyline)
                demo = new PiecewiseLinearPopulation(Units.Type.YEARS);
            } else if (demographicModel == 9) { // Piecewise Linear (skyride)
                demo = new PiecewiseLinearPopulation(Units.Type.YEARS);
            }

            guessArguments();

            if (demographicModel >= 7) { // is skyline or skyride
                popSizeCount = getTraceRange(argumentTraces[0]);
                if (demographicModel < 9) { // is skyline not skyride
                    groupSizeCount = getTraceRange(argumentTraces[1]);
                }
            }

        }

        if (traces != null) {
            double value;
            if (demographicModel == 0) { // Constant Size
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((ConstantPopulation) demo).setN0(value);
            } else if (demographicModel == 1) { // Exponential Growth (Growth Rate)
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((ExponentialGrowth) demo).setN0(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][1]], current);
                ((ExponentialGrowth) demo).setGrowthRate(value);
            } else if (demographicModel == 2) { // Exponential Growth (Doubling Time)
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((ExponentialGrowth) demo).setN0(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][1]], current);
                ((ExponentialGrowth) demo).setDoublingTime(value);
            } else if (demographicModel == 3) { // Logistic Growth (Growth Rate)
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((LogisticGrowth) demo).setN0(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][1]], current);
                ((LogisticGrowth) demo).setGrowthRate(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][2]], current);
                ((LogisticGrowth) demo).setTime50(value);
            } else if (demographicModel == 4) { // Logistic Growth (Doubling Time)
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((LogisticGrowth) demo).setN0(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][1]], current);
                ((LogisticGrowth) demo).setDoublingTime(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][2]], current);
                ((LogisticGrowth) demo).setTime50(value);
            } else if (demographicModel == 5) { // Expansion (Growth Rate)
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((Expansion) demo).setN0(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][1]], current);
                ((Expansion) demo).setProportion(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][2]], current);
                ((Expansion) demo).setGrowthRate(value);
            } else if (demographicModel == 6) { // Expansion (Doubling Time)
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][0]], current);
                ((Expansion) demo).setN0(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][1]], current);
                ((Expansion) demo).setProportion(value);
                value = traces.getStateValue(argumentTraces[argumentIndices[demographicModel][2]], current);
                ((Expansion) demo).setDoublingTime(value);
            } else if (demographicModel >= 7) { // Piecewise (skyline/skyride)
            }
        } else {
            double value;
            if (demographicModel == 0) { // Constant Size
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((ConstantPopulation) demo).setN0(value);
            } else if (demographicModel == 1) { // Exponential Growth (Growth Rate)
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((ExponentialGrowth) demo).setN0(value);
                value = argumentValues[argumentIndices[demographicModel][1]];
                ((ExponentialGrowth) demo).setGrowthRate(value);
            } else if (demographicModel == 2) { // Exponential Growth (Doubling Time)
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((ExponentialGrowth) demo).setN0(value);
                value = argumentValues[argumentIndices[demographicModel][1]];
                ((ExponentialGrowth) demo).setDoublingTime(value);
            } else if (demographicModel == 3) { // Logistic Growth (Growth Rate)
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((LogisticGrowth) demo).setN0(value);
                value = argumentValues[argumentIndices[demographicModel][1]];
                ((LogisticGrowth) demo).setGrowthRate(value);
                value = argumentValues[argumentIndices[demographicModel][2]];
                ((LogisticGrowth) demo).setTime50(value);
            } else if (demographicModel == 4) { // Logistic Growth (Doubling Time)
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((LogisticGrowth) demo).setN0(value);
                value = argumentValues[argumentIndices[demographicModel][1]];
                ((LogisticGrowth) demo).setDoublingTime(value);
                value = argumentValues[argumentIndices[demographicModel][2]];
                ((LogisticGrowth) demo).setTime50(value);
            } else if (demographicModel == 5) { // Expansion (Growth Rate)
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((Expansion) demo).setN0(value);
                value = argumentValues[argumentIndices[demographicModel][1]];
                ((Expansion) demo).setProportion(value);
                value = argumentValues[argumentIndices[demographicModel][2]];
                ((Expansion) demo).setGrowthRate(value);
            } else if (demographicModel == 6) { // Expansion (Doubling Time)
                value = argumentValues[argumentIndices[demographicModel][0]];
                ((Expansion) demo).setN0(value);
                value = argumentValues[argumentIndices[demographicModel][1]];
                ((Expansion) demo).setProportion(value);
                value = argumentValues[argumentIndices[demographicModel][2]];
                ((Expansion) demo).setDoublingTime(value);
            } else if (demographicModel >= 7) { // Piecewise Constant (skyline)
                throw new IllegalArgumentException("Fixed epoch model not supported");
            }
        }
        current++;

        return demo;
    }

    public void setupSkyline() {
        guessArguments();
        popSizeCount = getTraceRange(argumentTraces[0]);
        if (demographicModel < 9) { // is skyline not skyride
            groupSizeCount = getTraceRange(argumentTraces[1]);
        }     
    }

    public void getNextSkyline(double[] popSizes, double[] groupSizes) {
        int popSizeIndex = argumentTraces[0];
        for (int i = 0; i < popSizeCount; i++) {
            popSizes[i] = traces.getStateValue(popSizeIndex + i, current);
        }

        if (groupSizes != null) { // is skyline not skyride
            int groupSizeIndex = argumentTraces[1];
            for (int i = 0; i < groupSizeCount; i++) {
                groupSizes[i] = traces.getStateValue(groupSizeIndex + i, current);
            }
        }
        current++;
    }
    
    public void reset() {
        current = 0;
    }

    private void guessArguments() {
        int argumentCount = argumentIndices[demographicModel].length;

        for (int i = 0; i < argumentCount; i++) {

            int index = -1;

            for (int j = 0; j < argumentGuesses[i].length; j++) {
                if (index != -1) break;

                index = findArgument(argumentGuesses[argumentIndices[demographicModel][i]][j]);
            }

            argumentTraces[i] = index;
        }
    }

    private int findArgument(String argument) {
        for (int i = 0; i < traces.getTraceCount(); i++) {
            String item = traces.getTraceName(i).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    private int getTraceRange(int first) {
        int i = 1;
        int k = first;

        String name = traces.getTraceName(first);
        String root = name.substring(0, name.length() - 1);
        while (k < traces.getTraceCount() && traces.getTraceName(k).equals(root + i)) {
            i++;
            k++;
        }

        return i - 1;
    }
}

