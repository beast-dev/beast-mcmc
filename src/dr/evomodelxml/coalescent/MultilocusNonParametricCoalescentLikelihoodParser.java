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

package dr.evomodelxml.coalescent;

import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.MultilocusNonparametricCoalescentLikelihood;
import dr.evomodel.coalescent.NewMultilocusNonparametricCoalescentLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.coalescent.smooth.SkyGlideLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.smooth.SmoothSkygridLikelihoodParser;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class MultilocusNonParametricCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "multiLocusNPCoalescentLikelihood";
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
        Parameter gridPoints = SmoothSkygridLikelihoodParser.getGridPoints(xo); //TODO should I call this parser?
        int nGridPoints = gridPoints.getDimension();

        //        cxo = xo.getChild(PRECISION_PARAMETER);
        //        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        List<TreeModel> trees = new ArrayList<>();
        XMLObject cxo = xo.getChild(POPULATION_TREE);
        for (int i = 0; i < cxo.getChildCount(); i++) {
            trees.add((TreeModel) cxo.getChild(i));
        }

        List<BigFastTreeIntervals> intervalLists = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals(trees.get(i));
            intervalLists.add(treeIntervals);
        }

        Parameter ploidyFactors = parsePloidyFactors(xo, trees, nGridPoints);

        MultilocusNonparametricCoalescentLikelihood likelihood = new MultilocusNonparametricCoalescentLikelihood(intervalLists, logPopSizes, gridPoints, ploidyFactors);

        List<TreeIntervals> newIntervalLists = new ArrayList<>();
        for (int i = 0; i < trees.size(); ++i) {
            TreeIntervals treeIntervals = new TreeIntervals(trees.get(i));
            newIntervalLists.add(treeIntervals);
        }

        NewMultilocusNonparametricCoalescentLikelihood newLikelihood = new NewMultilocusNonparametricCoalescentLikelihood(
                newIntervalLists,logPopSizes, gridPoints, ploidyFactors);


        System.err.println("old: " + likelihood.getLogLikelihood());
        System.err.println("new: " + newLikelihood.getLogLikelihood());
        System.exit(-1);

        return likelihood;
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

            new XORRule(
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
            ),

    };

    @Override
    public String getParserDescription() {
        return "A likelihood for the multilocus nonparametric coalescent model";
    }

    @Override
    public Class getReturnType() {
        return MultilocusNonparametricCoalescentLikelihood.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}