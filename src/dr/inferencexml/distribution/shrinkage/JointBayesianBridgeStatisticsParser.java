/*
 * JointBayesianBridgeStatisticsParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
package dr.inferencexml.distribution.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.JointBayesianBridgeStatistics;
import dr.xml.*;

import java.util.List;

/**
 * @author Alexander Fisher
 */

public class JointBayesianBridgeStatisticsParser extends AbstractXMLObjectParser {
    private static final String JOINT_BAYESIAN_BRIDGE_STATISTICS = "jointBayesianBridgeStatistics";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<AutoCorrelatedBranchRatesDistribution> incProviderList;

        incProviderList = xo.getAllChildren(AutoCorrelatedBranchRatesDistribution.class);
        BayesianBridgeDistributionModel prior;
        BayesianBridgeDistributionModel previousPrior;
        for (int i = 0; i < incProviderList.size(); i++) {
            if (!(incProviderList.get(i).getPrior() instanceof BayesianBridgeDistributionModel)) {
                throw new XMLParseException("Must supply rates model with a BayesianBridgeDistributionModel prior");
            }

            prior = (BayesianBridgeDistributionModel) incProviderList.get(i).getPrior();
            if (i > 0) {
                previousPrior = (BayesianBridgeDistributionModel) incProviderList.get(i - 1).getPrior();
                if (!(bayesianBridgePriorsEqual(prior, previousPrior))) {
                    throw new XMLParseException("Bayesian bridge prior must be identical across rate models.");
                }
            }

        }

        return new JointBayesianBridgeStatistics(incProviderList);
    }

    private boolean bayesianBridgePriorsEqual(BayesianBridgeDistributionModel prior, BayesianBridgeDistributionModel previousPrior) {
        boolean globalScaleEqual = (prior.getGlobalScale() == previousPrior.getGlobalScale());
        boolean localScaleEqual = (prior.getLocalScale() == previousPrior.getLocalScale());
        boolean exponentEqual = (prior.getExponent() == previousPrior.getExponent());
        boolean slabEqual = (prior.getSlabWidth() == previousPrior.getSlabWidth());

        return (globalScaleEqual && localScaleEqual && exponentEqual && slabEqual);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AutoCorrelatedBranchRatesDistribution.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return JointBayesianBridgeStatistics.class;
    }

    @Override
    public String getParserName() {
        return JOINT_BAYESIAN_BRIDGE_STATISTICS;
    }

}
