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
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.distribution.shrinkage.JointBayesianBridgeStatistics;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider.equivalent;

/**
 * @author Alexander Fisher
 */

public class JointBayesianBridgeStatisticsParser extends AbstractXMLObjectParser {

    private static final String JOINT_BAYESIAN_BRIDGE_STATISTICS = "jointBayesianBridgeStatistics";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<BayesianBridgeStatisticsProvider> providers =
                xo.getAllChildren(BayesianBridgeStatisticsProvider.class);

        BayesianBridgeStatisticsProvider base = providers.get(0);
        for (int i = 1; i < providers.size(); ++i) {
            BayesianBridgeStatisticsProvider next = providers.get(i);

            if (!equivalent(base, next)) {
                    throw new XMLParseException("Bayesian bridge prior must be identical across providers.");
                }
            }

        return new JointBayesianBridgeStatistics(providers);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BayesianBridgeStatisticsProvider.class, 1, Integer.MAX_VALUE),
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
