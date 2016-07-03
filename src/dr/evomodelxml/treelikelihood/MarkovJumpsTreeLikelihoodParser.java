/*
 * MarkovJumpsTreeLikelihoodParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.inference.markovjumps.MarkovJumpsRegisterAcceptor;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Map;
import java.util.Set;

/**
 * @author Marc Suchard
 */

public class MarkovJumpsTreeLikelihoodParser extends AncestralStateTreeLikelihoodParser {

    public static final String MARKOV_JUMP_TREE_LIKELIHOOD = "markovJumpsTreeLikelihood";
    public static final String JUMP_TAG = "jumps";
    public static final String JUMP_TAG_NAME = "jumpTagName";
    public static final String COUNTS = MarkovJumpsType.COUNTS.getText();
    public static final String REWARDS = MarkovJumpsType.REWARDS.getText();
    public static final String SCALE_REWARDS = "scaleRewardsByTime";
    public static final String USE_UNIFORMIZATION = "useUniformization";
    public static final String SAVE_HISTORY = "saveCompleteHistory";
    public static final String LOG_HISTORY = "logCompleteHistory";
    public static final String COMPACT_HISTORY = "compactHistory";
    public static final String NUMBER_OF_SIMULANTS = "numberOfSimulants";
    public static final String REPORT_UNCONDITIONED_COLUMNS = "reportUnconditionedValues";


    public String getParserName() {
        return MARKOV_JUMP_TREE_LIKELIHOOD;
    }

    protected BeagleTreeLikelihood createTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                                        BranchModel branchModel,
                                                        GammaSiteRateModel siteRateModel,
                                                        BranchRateModel branchRateModel,
                                                        TipStatesModel tipStatesModel,
                                                        boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                        boolean delayScaling,
                                                        Map<Set<String>, Parameter> partialsRestrictions,
                                                        XMLObject xo) throws XMLParseException {

        DataType dataType = branchModel.getRootSubstitutionModel().getDataType();

        String stateTag = xo.getAttribute(RECONSTRUCTION_TAG_NAME,RECONSTRUCTION_TAG);
        String jumpTag = xo.getAttribute(JUMP_TAG_NAME, JUMP_TAG);

        boolean scaleRewards = xo.getAttribute(SCALE_REWARDS,true);

        boolean useMAP = xo.getAttribute(MAP_RECONSTRUCTION, false);
        boolean useMarginalLogLikelihood = xo.getAttribute(MARGINAL_LIKELIHOOD, true);

        boolean useUniformization = xo.getAttribute(USE_UNIFORMIZATION, false);
        boolean reportUnconditionedColumns = xo.getAttribute(REPORT_UNCONDITIONED_COLUMNS, false);
        int nSimulants = xo.getAttribute(NUMBER_OF_SIMULANTS, 1);

        if (patternList.areUnique()) {
            throw new XMLParseException("Markov Jumps reconstruction cannot be used with compressed (unique) patterns.");
        }

        MarkovJumpsBeagleTreeLikelihood treeLikelihood = new MarkovJumpsBeagleTreeLikelihood(
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                delayScaling,
                partialsRestrictions,
                dataType,
                stateTag,
                useMAP,
                useMarginalLogLikelihood,
                useUniformization,
                reportUnconditionedColumns,
                nSimulants
        );

        int registersFound = parseAllChildren(xo, treeLikelihood, dataType.getStateCount(), jumpTag,
                MarkovJumpsType.COUNTS, false); // For backwards compatibility

        XMLObject cxo = xo.getChild(COUNTS);
        if (cxo != null) {
            registersFound += parseAllChildren(cxo, treeLikelihood, dataType.getStateCount(), jumpTag,
                    MarkovJumpsType.COUNTS, false);
        }

        cxo = xo.getChild(REWARDS);
        if (cxo != null) {
            registersFound += parseAllChildren(cxo, treeLikelihood, dataType.getStateCount(), jumpTag,
                    MarkovJumpsType.REWARDS, scaleRewards);
        }

        if (registersFound == 0) { // Some default values for testing
//            double[] registration = new double[dataType.getStateCount()*dataType.getStateCount()];
//            MarkovJumpsCore.fillRegistrationMatrix(registration,dataType.getStateCount()); // Count all transitions
//            Parameter registerParameter = new Parameter.Default(registration);
//            registerParameter.setId(jumpTag);
//            treeLikelihood.addRegister(registerParameter,
//                                       MarkovJumpsType.COUNTS,
//                                       false);
            // Do nothing, should run the same as AncestralStateBeagleTreeLikelihood
        }

        boolean saveCompleteHistory = xo.getAttribute(SAVE_HISTORY, false);
        if (saveCompleteHistory) {
            Parameter allCounts = new Parameter.Default(dataType.getStateCount() * dataType.getStateCount());
            for (int i = 0; i < dataType.getStateCount(); ++i) {
                for (int j = 0; j < dataType.getStateCount(); ++j) {
                    if (j == i) {
                        allCounts.setParameterValue(i * dataType.getStateCount() + j, 0.0);
                    } else {
                        allCounts.setParameterValue(i * dataType.getStateCount() + j, 1.0);
                    }
                }
            }
            allCounts.setId(MarkovJumpsBeagleTreeLikelihood.TOTAL_COUNTS);
            treeLikelihood.setLogHistories(xo.getAttribute(LOG_HISTORY, false));
            treeLikelihood.setUseCompactHistory(xo.getAttribute(COMPACT_HISTORY, false));
            treeLikelihood.addRegister(allCounts, MarkovJumpsType.HISTORY, false);
        }

        return treeLikelihood;
    }

    public static int parseAllChildren(XMLObject xo,
                                       MarkovJumpsRegisterAcceptor acceptor,
                                       int stateCount,
                                       String jumpTag,
                                       MarkovJumpsType type,
                                       boolean scaleRewards) throws XMLParseException {
        int registersFound = 0;
        for(int i = 0; i < xo.getChildCount(); i++) {
            Object obj = xo.getChild(i);
            if (obj instanceof Parameter) {
                Parameter registerParameter = (Parameter) obj;
                if (type == MarkovJumpsType.COUNTS &&
                        registerParameter.getDimension() != stateCount * stateCount) {
                    if (registerParameter. getDimension() == 1) {
                        // if the dimension hasn't been set then default to counting all jumps
                        registerParameter.setDimension(stateCount * stateCount);
                        for (int j = 0; j < stateCount; j++) {
                            for (int k = 0; k < stateCount; k++) {
                                registerParameter.setParameterValueQuietly((j * stateCount) + k,  (j == k ? 0.0 : 1.0));
                            }
                        }
                    } else {
                        throw new XMLParseException("Markov Jumps register parameter " + registerParameter.getId() + " is of the wrong dimension");
                    }
                }
                if (type == MarkovJumpsType.REWARDS &&
                        registerParameter.getDimension() != stateCount) {
                    if (registerParameter.getDimension() == 1) {
                        // if the dimension hasn't been set then default to getting rewards for all states
                        registerParameter.setDimension(stateCount);
                        for (int j = 0; j < stateCount; j++) {
                            registerParameter.setParameterValueQuietly(j, 1.0);
                        }
                    } else {
                        throw new XMLParseException("Markov Rewards register parameter " + registerParameter.getId() + " is of the wrong dimension");
                    }
                }

                if (registerParameter.getId() == null) {
                    registerParameter.setId(jumpTag+(registersFound+1));
                }
                acceptor.addRegister(registerParameter, type, scaleRewards);
                registersFound++;
            }
        }
        return registersFound;
    }

    public static XMLSyntaxRule[]  rules =
            new XMLSyntaxRule[] {
                    AttributeRule.newBooleanRule(BeagleTreeLikelihoodParser.USE_AMBIGUITIES, true),
                    AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
                    AttributeRule.newStringRule(JUMP_TAG_NAME, true),
                    AttributeRule.newBooleanRule(SCALE_REWARDS,true),
                    AttributeRule.newBooleanRule(USE_UNIFORMIZATION,true),
                    AttributeRule.newBooleanRule(REPORT_UNCONDITIONED_COLUMNS, true),
                    AttributeRule.newIntegerRule(NUMBER_OF_SIMULANTS,true),
                    AttributeRule.newBooleanRule(SAVE_HISTORY, true),
                    AttributeRule.newBooleanRule(LOG_HISTORY, true),
                    AttributeRule.newBooleanRule(COMPACT_HISTORY, true),
                    new ElementRule(PARTIALS_RESTRICTION, new XMLSyntaxRule[] {
                            new ElementRule(TaxonList.class),
                            new ElementRule(Parameter.class),
                    }, true),
                    new ElementRule(PatternList.class),
                    new ElementRule(TreeModel.class),
                    new ElementRule(GammaSiteRateModel.class),
                    new ElementRule(BranchModel.class, true),
                    new ElementRule(SubstitutionModel.class, true),
                    new ElementRule(BranchRateModel.class, true),
                    AttributeRule.newStringRule(BeagleTreeLikelihoodParser.SCALING_SCHEME, true),
                    new ElementRule(Parameter.class,0,Integer.MAX_VALUE), // For backwards compatibility
                    new ElementRule(COUNTS,
                            new XMLSyntaxRule[] {
                                    new ElementRule(Parameter.class,0,Integer.MAX_VALUE)
                            },true),
                    new ElementRule(REWARDS,
                            new XMLSyntaxRule[] {
                                    new ElementRule(Parameter.class,0,Integer.MAX_VALUE)
                            },true),
                    new ElementRule(FrequencyModel.class, true),
            };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}