/*
 * ExponentialSkythingLikelihood.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.ExponentialBSPGrowth;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.BayesianSkylineLikelihoodParser;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.math.MathUtils;

import java.util.Date;

/**
 * A likelihood function for a sky**** function where exponential growth or decay parameters are what are piecewise
 * constant
 *
 * @version $Id: BayesianSkylineLikelihood.java,v 1.5 2006/03/06 11:26:49 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class ExponentialSkythingLikelihood extends OldAbstractCoalescentLikelihood {

    // PUBLIC STUFF

    public ExponentialSkythingLikelihood(Tree tree,
                                         Parameter slopeParameter,
                                         Parameter startingPopSize,
                                         int type) {
        super(BayesianSkylineLikelihoodParser.SKYLINE_LIKELIHOOD);

        groupSizeParameter = new Parameter.Default(slopeParameter.getDimension(), 1);
        popSizeParameter = new Parameter.Default(slopeParameter.getDimension()-1);

        this.slopeParameter = slopeParameter;
        this.startingPopSize = startingPopSize;
        int events = tree.getExternalNodeCount() - 1;
        int paramDim = slopeParameter.getDimension();
        this.type = type;

        if (paramDim != events) {
            throw new IllegalArgumentException("There are more groups than coalescent nodes in the tree.");
        }

        this.tree = tree;
        if (tree instanceof TreeModel) {
            addModel((TreeModel)tree);
        }
        addVariable(slopeParameter);

        setupIntervals();

        addStatistic(new GroupHeightStatistic());
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double getLogLikelihood() {

        setupIntervals();

        double logL = 0.0;

        double currentTime = 0.0;

        int groupIndex=0;
        int[] groupSizes = getGroupSizes();
        double[] groupEnds = getGroupHeights();

        int subIndex = 0;


        ExponentialBSPGrowth eg = new ExponentialBSPGrowth(Type.YEARS);

        double startGroupPopSize = 0;
        double endGroupPopSize;

        for (int j = intervalCount - 1; j >=0; j--) {
            if(j==intervalCount - 1){
                endGroupPopSize = startingPopSize.getParameterValue(0);
            } else {
                endGroupPopSize = startGroupPopSize;
            }
            double startTime = currentTime;
            double endTime = currentTime + intervals[j];
            startGroupPopSize = endGroupPopSize*Math.exp(slopeParameter.getParameterValue(j)*(endTime - startTime));


            eg.setupN1(endGroupPopSize, slopeParameter.getParameterValue(j), endTime - startTime);

            if (getIntervalType(j) == CoalescentEventType.COALESCENT) {
                subIndex += 1;
                if (subIndex >= groupSizes[groupIndex]) {
                    groupIndex += 1;
                    subIndex = 0;
                }
            }

            logL += calculateIntervalLikelihood(eg, intervals[j], currentTime, lineageCounts[j], getIntervalType(j));

            // insert zero-length coalescent intervals
            int diff = getCoalescentEvents(j)-1;
            for (int k = 0; k < diff; k++) {
                eg.setup(startGroupPopSize, startGroupPopSize, endTime - startTime);
                logL += calculateIntervalLikelihood(eg, 0.0, currentTime, lineageCounts[j]-k-1,
                        CoalescentEventType.COALESCENT);
                subIndex += 1;
                if (subIndex >= groupSizes[groupIndex]) {
                    groupIndex += 1;
                    subIndex = 0;
                }
            }

            currentTime += intervals[j];
        }

        return logL;
    }

    /**
     * @return the pop size for the given time. If linear model is being used then this pop size is
     * interpolated between the two pop sizes at either end of the grouped interval.
     */
    public final double getPopSize(int groupIndex, double midTime, double[] groupHeights) {

        return popSizeParameter.getParameterValue(groupIndex);

    }

    /* GAL: made public to give BayesianSkylineGibbsOperator access */
    public final int[] getGroupSizes() {
        if (groupSizeParameter.getParameterValue(0) < 2.0) {
            throw new IllegalArgumentException("For linear model first group size must be >= 2.");
        }

        int[] groupSizes = new int[groupSizeParameter.getDimension()];

        for (int i = 0; i < groupSizes.length; i++) {
            double g = groupSizeParameter.getParameterValue(i);
            if (g != Math.round(g)) {
                throw new RuntimeException("Group size " + i + " should be integer but found:" + g);
            }
            groupSizes[i] = (int)Math.round(g);
        }
        return groupSizes;
    }

    private  int getGroupCount() {
        return groupSizeParameter.getDimension();
    }

    private  int getGroupSize(int groupIndex) {
        double g = groupSizeParameter.getParameterValue(groupIndex);
        if (g != Math.round(g)) {
            throw new RuntimeException("Group size " + groupIndex + " should be integer but found:" + g);
        }
        return (int)Math.round(g);
    }

    /* GAL: made public to give BayesianSkylineGibbsOperator access */
    public final double[] getGroupHeights() {
        double[] groupEnds = new double[getGroupCount()];

        double timeEnd = 0.0;
        int groupIndex = 0;
        int subIndex = 0;
        for (int i = 0; i < intervalCount; i++) {

            timeEnd += intervals[i];

            if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                subIndex += 1;
                if (subIndex >= getGroupSize(groupIndex)) {
                    groupEnds[groupIndex] = timeEnd;
                    groupIndex += 1;
                    subIndex = 0;
                }
            }
        }
        groupEnds[getGroupCount()-1] = timeEnd;

        return groupEnds;
    }

    private double getGroupHeight(int groupIndex) {
        return getGroupHeights()[groupIndex];
    }

    final public int getType() {
        return type;
    }

    final public Parameter getPopSizeParameter() {
        return popSizeParameter;
    }

    final public Parameter getGroupSizeParameter() {
        return groupSizeParameter;
    }

    // ****************************************************************
    // Implementing Demographic Reconstructor
    // ****************************************************************

    public String getTitle() {
        final String title = "Bayesian Skything (exponential)\n" +
                "Generated " + (new Date()).toString() + " [seed=" + MathUtils.getSeed() + "]";
        return title;
    }

    // ****************************************************************
    // Inner classes
    // ****************************************************************

    public class GroupHeightStatistic extends Statistic.Abstract {

        public GroupHeightStatistic() {
            super("groupHeight");
        }

        public int getDimension() { return getGroupCount(); }

        public double getStatisticValue(int i) {
            return getGroupHeight(i);
        }

    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /** The demographic model. */
    private final Parameter slopeParameter;

    private final Parameter popSizeParameter;

    private final Parameter groupSizeParameter;

    private final Parameter startingPopSize;

    private final int type;

}
