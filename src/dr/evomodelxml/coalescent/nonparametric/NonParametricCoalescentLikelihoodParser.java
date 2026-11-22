/*
 * GMRFSkyrideLikelihoodParser.java
 *
 * Copyright (c) 2002-2024 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.coalescent.nonparametric;

import dr.evolution.coalescent.TreeIntervalList;
import dr.evomodelxml.coalescent.skyride.GMRFSkyrideLikelihoodParser;
import dr.evomodel.coalescent.intervals.FastTreeIntervals;
import dr.evomodel.epochs.EpochProvider;
import dr.evomodel.coalescent.nonparametric.NonParametricCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.smooth.SmoothSkygridLikelihoodParser;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class NonParametricCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "NonParametricCoalescentLikelihood";
    private static final String POPULATION_TREE = GMRFSkyrideLikelihoodParser.POPULATION_TREE;
    private static final String POPULATION_PARAMETER = GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER;
    //    public static final String PRECISION_PARAMETER = "precisionParameter";

    public static final String PLOIDY = "ploidy";

    private static final String GRID_POINTS = GMRFSkyrideLikelihoodParser.GRID_POINTS;
    private static final String NUM_GRID_POINTS = GMRFSkyrideLikelihoodParser.NUM_GRID_POINTS;
    private static final String CUT_OFF = GMRFSkyrideLikelihoodParser.CUT_OFF;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject xoPopParameter = xo.getChild(POPULATION_PARAMETER);
        Parameter logPopSizes = (Parameter) xoPopParameter.getChild(Parameter.class);
        Transform.ParsedTransform pt = (Transform.ParsedTransform) xoPopParameter.getChild(Transform.ParsedTransform.class);
        if (pt != null) {
            for (int i = 0; i < logPopSizes.getDimension(); i++) {
                logPopSizes.setParameterValue(i, pt.transform.transform(logPopSizes.getParameterValue(i)));
            }
        }
        Parameter gridPoints;
        EpochProvider epochs = (EpochProvider) xo.getChild(EpochProvider.class);
        if (epochs != null) {
            if (xo.hasChildNamed(GRID_POINTS) || xo.hasChildNamed(NUM_GRID_POINTS) || xo.hasChildNamed(CUT_OFF)) {
                throw new XMLParseException("The " + PARSER_NAME + " element has both epochs and grid points: " +
                        "give one or the other");
            }
            if (logPopSizes.getDimension() != epochs.getEpochCount()) {
                throw new XMLParseException("The " + PARSER_NAME + " element has " + epochs.getEpochCount() +
                        " epochs but its population size parameter has dimension " + logPopSizes.getDimension() +
                        ": there must be one population size for each epoch");
            }
            // the transition times between the epochs are the grid points
            gridPoints = epochs.getTransitionTimeParameter();
        } else {
            gridPoints = parseGridPoints(xo);
        }
        int nGridPoints = gridPoints.getDimension();

        //        cxo = xo.getChild(PRECISION_PARAMETER);
        //        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        List<TreeModel> trees = new ArrayList<>();
        XMLObject cxo = xo.getChild(POPULATION_TREE);
        for (int i = 0; i < cxo.getChildCount(); i++) {
            trees.add((TreeModel) cxo.getChild(i));
        }

        List<TreeIntervalList> intervalLists = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            intervalLists.add(new FastTreeIntervals(trees.get(i)));
        }

        Parameter ploidyFactors = parsePloidyFactors(xo, trees, nGridPoints);

        NonParametricCoalescentLikelihood likelihood = new NonParametricCoalescentLikelihood(intervalLists, logPopSizes, gridPoints, ploidyFactors);

        return likelihood;
    }

    public static Parameter parseGridPoints(XMLObject xo) throws XMLParseException {
        Parameter gridPoints;
        if (xo.getChild(GRID_POINTS) != null) {
            gridPoints = (Parameter) xo.getElementFirstChild(GRID_POINTS);
        } else {
            int numGridPoints = (int) ((Parameter) xo.getElementFirstChild(NUM_GRID_POINTS)).getParameterValue(0);
            double cutOff = ((Parameter) xo.getElementFirstChild(CUT_OFF)).getParameterValue(0);
            double[] gridLocations = new double[numGridPoints];

            for (int pt = 0; pt < gridLocations.length; pt++) {
                gridLocations[pt] = (pt + 1) * (cutOff / numGridPoints);
            }
            gridPoints = new Parameter.Default(gridLocations);
        }
        return gridPoints;
    }


    protected Parameter parsePloidyFactors(XMLObject xo, List<TreeModel> trees, int nGridPoints) {
        Parameter ploidyFactors = null;
        if (xo.getChild(PLOIDY) != null) {
            XMLObject cxo = xo.getChild(PLOIDY);
            ploidyFactors = (Parameter) cxo.getChild(Parameter.class);
        } else {
            ploidyFactors = new Parameter.Default(PLOIDY, trees.size());
            for (int i = 0; i < trees.size(); i++) {
                ploidyFactors.setParameterValue(i, 1.0);
            }
//            if (nGridPoints != 0) {
//                ploidyFactors = new Parameter.Default(PLOIDY, nGridPoints);
//                for(int i = 0; i < nGridPoints; i++){
//                    ploidyFactors.setParameterValue(i, 1.0);
//                }
//            } else {
//                ploidyFactors = new Parameter.Default(PLOIDY, trees.size());
//                for (int i = 0; i < trees.size(); i++) {
//                    ploidyFactors.setParameterValue(i, 1.0);
//                }
//            }
        }
        return ploidyFactors;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
            }),

            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
                    new ElementRule(Transform.ParsedTransform.class, true)
            }),

            new XORRule(new XMLSyntaxRule[]{
                    new ElementRule(EpochProvider.class),
                    new ElementRule(GRID_POINTS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new AndRule(
                            new ElementRule(CUT_OFF, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            }),
                            new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
                    )
            }),

    };

    @Override
    public String getParserDescription() {
        return "A likelihood for the multilocus nonparametric coalescent model";
    }

    @Override
    public Class getReturnType() {
        return NonParametricCoalescentLikelihood.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}