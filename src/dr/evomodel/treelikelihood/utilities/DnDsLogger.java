/*
 * DnDsLogger.java
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

package dr.evomodel.treelikelihood.utilities;

import dr.evomodel.substmodel.CodonLabeling;
import dr.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.math.EmpiricalBayesPoissonSmoother;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsLogger implements Loggable {

    public DnDsLogger(String name, Tree tree, TreeTrait[] traits, boolean useSmoothing, boolean useDnMinusDs, boolean counts, boolean synonymous) {
        this.tree = tree;
        this.traits = traits;
        numberSites = getNumberSites();
        this.name = name;
        this.useSmoothing = useSmoothing;
        this.useDnMinusDs = useDnMinusDs;
        this.counts = counts;
        this.synonymous = synonymous;

        for (int i = 0; i < NUM_TRAITS; i++) {
            if (traits[i].getIntent() != TreeTrait.Intent.WHOLE_TREE) {
                throw new IllegalArgumentException("Only whole tree traits are currently supported in DnDsLogger");
            }
        }
    }

//    public LogColumn[] getColumns() {
//        LogColumn[] columns;
////        if (conditionalCounts) {
//            columns = new LogColumn[numberSites * 3];
////        } else {
////            columns = new LogColumn[numberSites];
////        }
//        int columnCount = 0;
////        if (!conditionalCounts){
//            for (int i = 0; i < numberSites; i++) {
//                columns[columnCount] = new SmoothedColumn(name, i);
//                columnCount ++;
//            }
////        } else {
//            for (int j = 0; j < numberSites; j++) {
//                //CN and CS are outputted per site
//                columns[columnCount] = new ConditionalColumn(conditionalNon, j, true);
//                columns[columnCount + 1] = new ConditionalColumn(conditionalSyn, j, false);
//                columnCount = columnCount + 2;
//            }
////        }
//        return columns;
//    }

    public LogColumn[] getColumns() {
        LogColumn[] columns;
        if (!counts){
            columns = new LogColumn[numberSites];
        } else {
            columns = new LogColumn[numberSites * 2];
        }
        int columnCount = 0;
        if (!counts){
            for (int i = 0; i < numberSites; i++) {
                columns[columnCount] = new SmoothedColumn(name, i);
                columnCount ++;
            }
        } else {
            for (int j = 0; j < numberSites; j++) {
                //CN and CS are outputted per site
                if (synonymous){
                    columns[columnCount] = new ConditionalColumn(conditionalSyn, j, false, true);
                } else {
                    columns[columnCount] = new ConditionalColumn(conditionalNon, j, true, true);
                }
                columnCount ++;
            }
            for (int k = 0; k < numberSites; k++) {
                if (synonymous){
                    columns[columnCount] = new ConditionalColumn(unconditionalSyn, k, false, false);
                }  else {
                    columns[columnCount] = new ConditionalColumn(unconditionalNon, k, true, false);
                }
                columnCount ++;

            }

        }
        return columns;
    }

    private class ConditionalColumn extends NumberColumn {
        private final int index;
        boolean nonsynonymous = false;
        boolean conditional = true;

        public ConditionalColumn(String label, int index, boolean N, boolean C) {
            super(label + "[" + (index+1) + "]");
            this.index = index;
            if (N){
                nonsynonymous = true;
            }
            if (!C){
                conditional = false;
            }
        }

        @Override
        public double getDoubleValue() {
            if (index == 0) {
                doSmoothing();
            }
            if (conditional){
                if (nonsynonymous){
                    return cachedValues[CN][index];
                }  else {
                    return cachedValues[CS][index];
                }
            } else {
                if (nonsynonymous){
                    return cachedValues[UN][index];
                }  else {
                    return cachedValues[US][index];
                }
            }
        }
    }

    private class SmoothedColumn extends NumberColumn {

        private final int index;

        public SmoothedColumn(String label, int index) {
            super(label + "[" + (index+1) + "]");
            this.index = index;
        }

        @Override
        public double getDoubleValue() {
            if (index == 0) { // Assumes that columns are accessed IN ORDER
                doSmoothing();
            }
            return doCalculation(index);
        }
    }

    private double doCalculation(int index) {
        double returnValue;
        if (!useDnMinusDs) {
            returnValue = (cachedValues[CN][index] / cachedValues[UN][index]) /
                    (cachedValues[CS][index] / cachedValues[US][index]);

        } else {
            returnValue =  (cachedValues[CN][index] / cachedValues[UN][index]) -
                    (cachedValues[CS][index] / cachedValues[US][index]);
        }
        return returnValue;
    }

    private int getNumberSites() {
        double[] values = (double[]) traits[0].getTrait(tree, tree.getRoot());
        return values.length;
    }

    private void doSmoothing() {

        if (cachedValues == null) {
            cachedValues = new double[NUM_TRAITS][];
        }

        for (int i = 0; i < NUM_TRAITS; i++) {
            if (useSmoothing) {
                cachedValues[i] = EmpiricalBayesPoissonSmoother.smooth((double[]) traits[i].getTrait(tree, tree.getRoot()));
            } else {
                cachedValues[i] = (double[]) traits[i].getTrait(tree, tree.getRoot());
            }
        }
    }

    private final TreeTrait[] traits;
    private final Tree tree;
    private final int numberSites;
    private final String name;
    private final String conditionalNon = "CN";
    private final String conditionalSyn = "CS";
    private final String unconditionalNon = "UN";
    private final String unconditionalSyn = "US";
    private final boolean useSmoothing;
    private final boolean useDnMinusDs;
    private final boolean counts;
    private final boolean synonymous;

    private final static int NUM_TRAITS = 4;
    private final static int CS = 0;
    private final static int US = 1;
    private final static int CN = 2;
    private final static int UN = 3;

    private double[][] cachedValues;

    public static String[] traitNames = new String[] {
            CodonPartitionedRobustCounting.SITE_SPECIFIC_PREFIX + CodonLabeling.SYN.getText(),
            CodonPartitionedRobustCounting.UNCONDITIONED_PREFIX + CodonLabeling.SYN.getText(),
            CodonPartitionedRobustCounting.SITE_SPECIFIC_PREFIX + CodonLabeling.NON_SYN.getText(),
            CodonPartitionedRobustCounting.UNCONDITIONED_PREFIX + CodonLabeling.NON_SYN.getText()
    };
}


