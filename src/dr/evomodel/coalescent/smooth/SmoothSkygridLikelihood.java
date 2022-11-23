/*
 * SmoothSkygridLikelihood.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.smooth;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.util.*;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Xiang Ji
 * @author Yuwei Bao
 * @author Marc A. Suchard
 */
public class SmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Citable, Reportable {

    private final List<TreeModel> trees;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;
    private final Parameter smoothRate;
    private final SmoothSkygridPopulationSizeInverse populationSizeInverse;
    private final OldSmoothLineageCount lineageCount;

    private final GlobalSigmoidSmoothFunction smoothFunction;

    private final List<BigFastTreeIntervals> intervalsList;

    public SmoothSkygridLikelihood(String name,
                                   List<TreeModel> trees,
                                   Parameter logPopSizeParameter,
                                   Parameter gridPointParameter,
                                   Parameter smoothRate) {
        super(name);
        this.trees = trees;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;
        this.smoothRate = smoothRate;
        this.smoothFunction = new GlobalSigmoidSmoothFunction();
        this.populationSizeInverse = new SmoothSkygridPopulationSizeInverse(logPopSizeParameter, gridPointParameter, smoothFunction, smoothRate);
        this.lineageCount = new OldSmoothLineageCount(trees.get(0), smoothFunction, smoothRate);
        intervalsList = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            intervalsList.add(new BigFastTreeIntervals(trees.get(i)));
        }
    }

    @Override
    public String getReport() {
        return "smoothSkygrid(" + getLogLikelihood() + ")";
    }

    class OldSmoothLineageCount {

        private final Tree tree;
        private final GlobalSigmoidSmoothFunction smoothFunction;

        private final Parameter smoothRate;

        OldSmoothLineageCount(Tree tree, GlobalSigmoidSmoothFunction smoothFunction, Parameter smoothRate) {
            this.tree = tree;
            this.smoothFunction = smoothFunction;
            this.smoothRate = smoothRate;
        }

        double getLineageCount(double time) {
            double sum = 0;
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                sum += smoothFunction.getSmoothValue(time, tree.getNodeHeight(tree.getNode(i)), 0.0, 1.0, smoothRate.getParameterValue(0));
            }
            for (int i = tree.getExternalNodeCount(); i < tree.getNodeCount(); i++) {
                sum += smoothFunction.getSmoothValue(time, tree.getNodeHeight(tree.getNode(i)), 0.0, -1.0, smoothRate.getParameterValue(0));
            }
            return sum;
        }
    }

    class SmoothSkygridPopulationSizeInverse {

        private final Parameter logPopSizeParameter;
        private final Parameter gridPointParameter;
        private final GlobalSigmoidSmoothFunction smoothFunction;
        private final Parameter smoothRate;

        SmoothSkygridPopulationSizeInverse(Parameter logPopSizeParameter,
                                           Parameter gridPointParameter,
                                           GlobalSigmoidSmoothFunction smoothFunction,
                                           Parameter smoothRate) {
            this.logPopSizeParameter = logPopSizeParameter;
            this.gridPointParameter = gridPointParameter;
            this.smoothRate = smoothRate;
            this.smoothFunction = smoothFunction;
        }

        double getPopulationSizeInverse(double time) {
            double sum = 0;
            for(int i = 0; i < gridPointParameter.getDimension(); i++) {
                double increment = smoothFunction.getSmoothValue(time, gridPointParameter.getParameterValue(i),
                        i == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(0)) : 0.0,
                        i == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(1)) : Math.exp(-logPopSizeParameter.getParameterValue(i + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(i)),
                        smoothRate.getParameterValue(0));
                sum += increment;
            }
            return sum;
        }
    }

    @Override
    public Type getUnits() {
        return null;
    }

    @Override
    public void setUnits(Type units) {

    }
    @Override
    protected double calculateLogLikelihood() {
        assert(trees.size() == 1);
        Tree tree = trees.get(0);
        BigFastTreeIntervals intervals = intervalsList.get(0);
        double logPopulationSizeInverse = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
//            logPopulationSizeInverse += Math.log(populationSizeInverse.getPopulationSizeInverse(tree.getNodeHeight(node)));
            logPopulationSizeInverse += Math.log(getSmoothPopulationSizeInverse(tree.getNodeHeight(node), tree.getNodeHeight(tree.getRoot())));
        }

        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());

        final double firstInversePopulationSize = Math.exp(-logPopSizeParameter.getParameterValue(0));
        double singleSigmoidIntegralSums =  ((double) intervals.getLineageCount(0) * (intervals.getLineageCount(0) - 1))
                * smoothFunction.getSingleIntegration(startTime, endTime, intervals.getIntervalTime(0), smoothRate.getParameterValue(0));
        for (int i = 1; i < intervals.getIntervalCount(); i++) {
            final double lineageCountDifference = ((double) intervals.getLineageCount(i) * (intervals.getLineageCount(i) - 1)
                    - intervals.getLineageCount(i - 1) * (intervals.getLineageCount(i - 1) - 1)) / 2;
            final double smoothIntegral = smoothFunction.getSingleIntegration(startTime, endTime, intervals.getIntervalTime(i), smoothRate.getParameterValue(0));
            singleSigmoidIntegralSums += lineageCountDifference * smoothIntegral;
        }
        singleSigmoidIntegralSums *= firstInversePopulationSize;

        int maxGridIndex = getMaxGridIndex(endTime);

        double pairSigmoidIntegralSums = 0;

        for (int j = 0; j < maxGridIndex + 1; j++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
            final double gridTime = gridPointParameter.getParameterValue(j);
            pairSigmoidIntegralSums += intervals.getLineageCount(0) * (intervals.getLineageCount(0) - 1) * (nextPopSizeInverse - currentPopSizeInverse) *
                    smoothFunction.getPairProductIntegration(startTime, endTime, intervals.getIntervalTime(0), gridTime, smoothRate.getParameterValue(0));
        }

        for (int i = 1; i < intervals.getIntervalCount(); i++) {

            final double lineageCountDifference = ((double) intervals.getLineageCount(i) * (intervals.getLineageCount(i) - 1)
                    - intervals.getLineageCount(i - 1) * (intervals.getLineageCount(i - 1) - 1)) / 2;
            final double intervalTime = intervals.getIntervalTime(i);

            for (int j = 0; j < maxGridIndex + 1; j++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
                final double gridTime = gridPointParameter.getParameterValue(j);
                pairSigmoidIntegralSums += lineageCountDifference * (nextPopSizeInverse - currentPopSizeInverse) *
                        smoothFunction.getPairProductIntegration(startTime, endTime, intervalTime, gridTime, smoothRate.getParameterValue(0));
            }
        }

        return logPopulationSizeInverse - singleSigmoidIntegralSums - pairSigmoidIntegralSums;
    }

    private int getMaxGridIndex(double endTime) {
        int maxGridIndex = gridPointParameter.getDimension() - 1;
        while (gridPointParameter.getParameterValue(maxGridIndex) > endTime && maxGridIndex > 0) {
            maxGridIndex--;
        }
        return maxGridIndex;
    }

    private double getSmoothPopulationSizeInverse(double t, double endTime) {
        int maxGridIndex = getMaxGridIndex(endTime);

        double populationSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(0));

        for (int j = 0; j < maxGridIndex + 1; j++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
            final double gridTime = gridPointParameter.getParameterValue(j);
            populationSizeInverse += (nextPopSizeInverse - currentPopSizeInverse) *
                    smoothFunction.getSmoothValue(t, gridTime, 0, 1, smoothRate.getParameterValue(0));
        }

        return populationSizeInverse;
    }

    private double oldCalculateLogLikelihood() {
        assert(trees.size() == 1);
        Tree tree = trees.get(0);
        double logPopulationSizeInverse = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            logPopulationSizeInverse += Math.log(populationSizeInverse.getPopulationSizeInverse(tree.getNodeHeight(node)));
        }
        double integralBit = 0;
        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double stepLocation1 = tree.getNodeHeight(tree.getNode(i));
            final double preStepValue1 = 0;
            final double postStepValue1 = i < tree.getExternalNodeCount() ? 1 : -1;
            for (int j = 0; j < tree.getNodeCount(); j++) {
                final double stepLocation2 = tree.getNodeHeight(tree.getNode(j));
                final double preStepValue2 = j == 0 ? -1 : 0;
                final double postStepValue2 = (j < tree.getExternalNodeCount() ? 1 : -1) + preStepValue2;
                for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                    final double stepLocation3 = gridPointParameter.getParameterValue(k);
                    final double preStepValue3 = k == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(0)) : 0;
                    final double postStepValue3 = k == 0? Math.exp(-logPopSizeParameter.getParameterValue(1)) :
                            Math.exp(-logPopSizeParameter.getParameterValue(k + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(k));
                    final double analytic = -0.5 * smoothFunction.getTripleProductIntegration(startTime, endTime,
                            stepLocation1, preStepValue1, postStepValue1,
                            stepLocation2, preStepValue2, postStepValue2,
                            stepLocation3, preStepValue3, postStepValue3,
                            smoothRate.getParameterValue(0));
                    integralBit += analytic;
                }
            }
        }
        return logPopulationSizeInverse + integralBit;
    }


    public static double getReciprocalPopSizeInInterval(double time, OldSmoothLineageCount lineageCount,
                                                        SmoothSkygridPopulationSizeInverse populationSizeInverse) {
        return 0.5 * lineageCount.getLineageCount(time) * (lineageCount.getLineageCount(time) - 1) * populationSizeInverse.getPopulationSizeInverse(time);
    }

    public OldSmoothLineageCount getLineageCount() {
        return lineageCount;
    }

    public SmoothSkygridPopulationSizeInverse getPopulationSizeInverse() {
        return populationSizeInverse;
    }

    @Override
    public int getNumberOfCoalescentEvents() {
        int sum = 0;
        for (Tree tree : trees) {
            sum += tree.getInternalNodeCount();
        }
        return sum;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        throw new RuntimeException("Not yet implemented.");
    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Differentiable skygrid coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(CommonCitations.GILL_2013_IMPROVING,
                new Citation(
                        new Author[] {
                                new Author( "Y", "Bao"),
                                new Author("MA", "Suchard"),
                                new Author( "X", "Ji"),
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
    }

}
