/*
 * CountableMixtureBranchRatesParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class TimeDependentBranchRatesParser extends AbstractXMLObjectParser {

    public static final String TIME_DEPENDENT_BRANCH_RATES = "timeDependentBranchRates";
    public static final String INTERCEPT = "intercept";
    public static final String TIME_COEFFICIENT = "branchTimeEffect";
    public static final String RANDOM_EFFECTS = "randomEffects";
    public static final String IN_LOG_SPACE = "inLogSpace";

    public String getParserName() {
        return TIME_DEPENDENT_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter intercept = (Parameter) xo.getElementFirstChild(INTERCEPT);
        Parameter timeCoefficient = (Parameter) xo.getElementFirstChild(TIME_COEFFICIENT);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        List<AbstractBranchRateModel> randomEffects = null;
        if (xo.hasChildNamed(RANDOM_EFFECTS)) {
            XMLObject cxo = xo.getChild(RANDOM_EFFECTS);
            randomEffects = new ArrayList<AbstractBranchRateModel>();
            for (int i = 0; i < cxo.getChildCount(); ++i) {
                randomEffects.add((AbstractBranchRateModel)cxo.getChild(i));
            }
        }

        boolean inLogSpace = xo.getAttribute(IN_LOG_SPACE, true);
        if (!inLogSpace){
            System.err.print("log space = false for Time Dependent Branch Rates not implemented yet");
            System.exit(-1);
        }

        Logger.getLogger("dr.evomodel").info("Using a time dependent branch rates model.");

        return new TimeDependentBranchRates(treeModel, intercept, timeCoefficient, randomEffects, inLogSpace);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a clock consisting of a mixture of fixed effect that is a function of the branch (midpoint) time and random effects.";
    }

    public Class getReturnType() {
        return CountableMixtureBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(INTERCEPT, Parameter.class, "The intercept for the branch time-dependent rate function", false),
            new ElementRule(TIME_COEFFICIENT, Parameter.class, "The coefficient for the branch time-dependent rate function", false),

            new ElementRule(RANDOM_EFFECTS,
                    new XMLSyntaxRule[] {
                            new ElementRule(AbstractBranchRateModel.class, 0, Integer.MAX_VALUE),
                    },
                    "Possible random effects", true),
            AttributeRule.newBooleanRule(IN_LOG_SPACE, true),
    };
}
