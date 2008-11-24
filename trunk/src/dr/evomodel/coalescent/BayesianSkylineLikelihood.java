/*
 * BayesianSkylineLikelihood.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.*;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;
import dr.math.MathUtils;

import java.util.logging.Logger;
import java.util.Date;

/**
 * A likelihood function for the generalized skyline plot coalescent. Takes a tree and population size and group size parameters.
 *
 * @version $Id: BayesianSkylineLikelihood.java,v 1.5 2006/03/06 11:26:49 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class BayesianSkylineLikelihood extends OldAbstractCoalescentLikelihood implements DemographicReconstructor {

    // PUBLIC STUFF

    public static final String SKYLINE_LIKELIHOOD = "generalizedSkyLineLikelihood";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String GROUP_SIZES = "groupSizes";

    public static final String TYPE = "type";
    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    public static final int STEPWISE_TYPE = 0;
    public static final int LINEAR_TYPE = 1;
    public static final int EXPONENTIAL_TYPE = 2;

    public BayesianSkylineLikelihood(Tree tree,
                                     Parameter popSizeParameter,
                                     Parameter groupSizeParameter,
                                     int type) {
        super(SKYLINE_LIKELIHOOD);

        this.groupSizeParameter = groupSizeParameter;
        this.popSizeParameter = popSizeParameter;
        int events = tree.getExternalNodeCount() - 1;
        int paramDim1 = popSizeParameter.getDimension();
        int paramDim2 = groupSizeParameter.getDimension();
        this.type = type;

        if (type == EXPONENTIAL_TYPE) {
            if (paramDim1 != (paramDim2+1)) {
                throw new IllegalArgumentException("Dimension of population parameter must be one greater than dimension of group size parameter.");
            }
        } else if (type == LINEAR_TYPE) {
            if (paramDim1 != (paramDim2+1)) {
                throw new IllegalArgumentException("Dimension of population parameter must be one greater than dimension of group size parameter.");
            }
        } else { // STEPWISE_TYPE
            if (paramDim1 != paramDim2) {
                throw new IllegalArgumentException("Dimension of population parameter and group size parameters should be the same.");
            }
        }

        if (paramDim2 > events) {
            throw new IllegalArgumentException("There are more groups than coalescent nodes in the tree.");
        }

        int eventsCovered = 0;
        for (int i = 0; i < getGroupCount(); i++) {
            eventsCovered += getGroupSize(i);
        }

        if (eventsCovered != events) {

            if (eventsCovered == 0 || eventsCovered == paramDim2) {
                double[] uppers = new double[paramDim2];
                double[] lowers = new double[paramDim2];

                // For these special cases we assume that the XML has not specified initial group sizes
                // or has set all to 1 and we set them here automatically...
                int eventsEach = events / paramDim2;
                int eventsExtras = events % paramDim2;
                for (int i = 0; i < paramDim2; i++) {
                    if (i < eventsExtras) {
                        groupSizeParameter.setParameterValue(i, eventsEach + 1);
                    } else {
                        groupSizeParameter.setParameterValue(i, eventsEach);
                    }
                    uppers[i] = Double.MAX_VALUE;
                    lowers[i] = 1.0;
                }

                if (type == EXPONENTIAL_TYPE || type == LINEAR_TYPE) {
                    lowers[0] = 2.0;
                }
                groupSizeParameter.addBounds(new Parameter.DefaultBounds(uppers, lowers));
            } else {
                // ... otherwise assume the user has made a mistake setting initial group sizes.
                throw new IllegalArgumentException("The sum of the initial group sizes does not match the number of coalescent events in the tree.");
            }
        }

        if ((type == EXPONENTIAL_TYPE || type == LINEAR_TYPE) && groupSizeParameter.getParameterValue(0) < 2.0) {
            throw new IllegalArgumentException("For linear or exponential model first group size must be >= 2.");
        }

        this.tree = tree;
        if (tree instanceof TreeModel) {
            addModel((TreeModel)tree);
        }
        addParameter(popSizeParameter);

        addParameter(groupSizeParameter);

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

        if (type == EXPONENTIAL_TYPE) {
            ExponentialBSPGrowth eg = new ExponentialBSPGrowth(Units.Type.YEARS);

            for (int j = 0; j < intervalCount; j++) {
                double startGroupPopSize = popSizeParameter.getParameterValue(groupIndex);
                double endGroupPopSize = popSizeParameter.getParameterValue(groupIndex+1);
                double startTime = currentTime;
                double endTime = currentTime + intervals[j];

                eg.setup(startGroupPopSize, endGroupPopSize, endTime - startTime);

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
        } else {
            ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);

            for (int j = 0; j < intervalCount; j++) {

                // set the population size to the size of the middle of the current interval
                final double ps = getPopSize(groupIndex, currentTime + (intervals[j]/2.0), groupEnds);
                cp.setN0(ps);
                if (getIntervalType(j) == CoalescentEventType.COALESCENT) {
                    subIndex += 1;
                    if (subIndex >= groupSizes[groupIndex]) {
                        groupIndex += 1;
                        subIndex = 0;
                    }
                }

                logL += calculateIntervalLikelihood(cp, intervals[j], currentTime, lineageCounts[j], getIntervalType(j));

                // insert zero-length coalescent intervals
                int diff = getCoalescentEvents(j)-1;
                for (int k = 0; k < diff; k++) {
                    cp.setN0(getPopSize(groupIndex, currentTime, groupEnds));
                    logL += calculateIntervalLikelihood(cp, 0.0, currentTime, lineageCounts[j]-k-1,
                            CoalescentEventType.COALESCENT);
                    subIndex += 1;
                    if (subIndex >= groupSizes[groupIndex]) {
                        groupIndex += 1;
                        subIndex = 0;
                    }
                }

                currentTime += intervals[j];
            }
        }
        return logL;
    }

    /**
     * @return the pop size for the given time. If linear model is being used then this pop size is
     * interpolated between the two pop sizes at either end of the grouped interval.
     */
    public final double getPopSize(int groupIndex, double midTime, double[] groupHeights) {
        if (type == LINEAR_TYPE) {

            double startGroupPopSize = popSizeParameter.getParameterValue(groupIndex);
            double endGroupPopSize = popSizeParameter.getParameterValue(groupIndex+1);

            double startGroupTime = 0.0;
            if (groupIndex > 0) {
                startGroupTime = groupHeights[groupIndex-1];
            }
            double endGroupTime = groupHeights[groupIndex];

            // calculate the gradient
            double m = (endGroupPopSize-startGroupPopSize)/(endGroupTime-startGroupTime);

            // calculate the population size at midTime using linear interpolation
            final double midPopSize = (m * (midTime-startGroupTime)) + startGroupPopSize;

            return midPopSize;
        } else {
            return popSizeParameter.getParameterValue(groupIndex);
        }
    }

    /* GAL: made public to give BayesianSkylineGibbsOperator access */
    public final int[] getGroupSizes() {
        if ((type == EXPONENTIAL_TYPE || type == LINEAR_TYPE) && groupSizeParameter.getParameterValue(0) < 2.0) {
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
        final String title = "Bayesian Skyline (" + getChangeType().name() + ")\n" +
                "Generated " + (new Date()).toString() + " [seed=" + MathUtils.getSeed() + "]";
        return title;
    }

    public ChangeType getChangeType() {
        switch (getType()) {
            case 0: return ChangeType.STEPWISE;
            case 1: return ChangeType.LINEAR;
            case 2: return ChangeType.EXPONENTIAL;
            default: throw new IllegalArgumentException("Unknown change type index");
        }
    }

    public double[] getIntervals() {
        return new double[0];
    }

    public double[] getPopSizes() {
        return new double[0];
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return SKYLINE_LIKELIHOOD; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject)xo.getChild(POPULATION_SIZES);
            Parameter param = (Parameter)cxo.getChild(Parameter.class);

            cxo = (XMLObject)xo.getChild(GROUP_SIZES);
            Parameter param2 = (Parameter)cxo.getChild(Parameter.class);

            cxo = (XMLObject)xo.getChild(CoalescentLikelihood.POPULATION_TREE);
            TreeModel treeModel = (TreeModel)cxo.getChild(TreeModel.class);

            int type = LINEAR_TYPE;
            String typeName = LINEAR;
            if (xo.hasAttribute(LINEAR) &&!xo.getBooleanAttribute(LINEAR)) {
                type = STEPWISE_TYPE;
                typeName = STEPWISE;
            }

            if (xo.hasAttribute(TYPE)) {
                if (xo.getStringAttribute(TYPE).equalsIgnoreCase(STEPWISE)) {
                    type = STEPWISE_TYPE;
                    typeName = STEPWISE;
                } else if (xo.getStringAttribute(TYPE).equalsIgnoreCase(LINEAR)) {
                    type = LINEAR_TYPE;
                    typeName = LINEAR;
                } else if (xo.getStringAttribute(TYPE).equalsIgnoreCase(EXPONENTIAL)) {
                    type = EXPONENTIAL_TYPE;
                    typeName = EXPONENTIAL;
                }
                else throw new XMLParseException("Unknown Bayesian Skyline type: " + xo.getStringAttribute(TYPE));
            }

            Logger.getLogger("dr.evomodel").info("Bayesian skyline plot: " + param.getDimension() + " " + typeName + " control points");

            return new BayesianSkylineLikelihood(treeModel, param, param2, type);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the population size vector.";
        }

        public Class getReturnType() { return Likelihood.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new XORRule(
                        AttributeRule.newBooleanRule(LINEAR),
                        AttributeRule.newStringRule(TYPE)
                ),
                new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(GROUP_SIZES, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(CoalescentLikelihood.POPULATION_TREE, new XMLSyntaxRule[] {
                        new ElementRule(TreeModel.class)
                }),
        };


    };

    /** The demographic model. */
    private final Parameter popSizeParameter;

    private final Parameter groupSizeParameter;

    private final int type;

}
