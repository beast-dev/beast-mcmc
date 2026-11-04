/*
 * EpochBranchBreakingStatisticParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.EpochBranchBreakingStatistic;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.xml.*;

/**
 */
public class EpochBranchBreakingStatisticParser extends AbstractXMLObjectParser {

    public static final String MODE = "mode";
    public static final String MEAN_MAX = "meanMax";
    public static final String PROP_BROKEN = "propBroken";
    public static final String EPOCH_BREAKING_STATISTIC = "EpochBranchBreakingStatistic";

    public String getParserName() {
        return EPOCH_BREAKING_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String mode = MEAN_MAX;
        if ( xo.hasAttribute(MODE) ) {
            mode = xo.getStringAttribute(MODE);
            if ( !mode.equals(MEAN_MAX) && !mode.equals(PROP_BROKEN) ) {
                throw new RuntimeException("Invalid mode for " + EPOCH_BREAKING_STATISTIC);
            }
        }

        EpochBranchModel epoch = (EpochBranchModel)xo.getChild(EpochBranchModel.class);

        return new EpochBranchBreakingStatistic(xo.getId().toString(), epoch, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns information about how an epoch model breaks up the branches of a tree.";
    }

    public Class getReturnType() {
        return EpochBranchBreakingStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(EpochBranchModel.class),
            new StringAttributeRule("mode", "What to return; propBroken for proportion of branches broken by an epoch time; meanMax for the average (across branches) maximum proportion of time spent in any one model. Default meanMax.", true),
    };

}
