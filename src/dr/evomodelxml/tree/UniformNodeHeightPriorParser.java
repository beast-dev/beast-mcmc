/*
 * UniformNodeHeightPriorParser.java
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

package dr.evomodelxml.tree;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.UniformNodeHeightPrior;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class UniformNodeHeightPriorParser extends AbstractXMLObjectParser {

    public static final String UNIFORM_ROOT_PRIOR = "uniformRootPrior";
    public static final String UNIFORM_NODE_HEIGHT_PRIOR = "uniformNodeHeightPrior";
    public static final String MAX_ROOT_HEIGHT = "maxRootHeight";
    public static final String ANALYTIC = "analytic";
    public static final String MC_SAMPLE = "mcSampleSize";
    public static final String MARGINAL = "marginal";
    public static final String LEADING_TERM = "approximate";

    public String getParserName() {
        return UNIFORM_NODE_HEIGHT_PRIOR;
    }

    public String[] getParserNames() {
        return new String[] {UNIFORM_ROOT_PRIOR, UNIFORM_NODE_HEIGHT_PRIOR};
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("\nConstructing a uniform node height prior:");

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        if (xo.hasAttribute(MAX_ROOT_HEIGHT)) {
            // the Nicholls & Gray variant
            double maxRootHeight = xo.getDoubleAttribute(MAX_ROOT_HEIGHT);
            Logger.getLogger("dr.evomodel").info("\tUsing joint variant with a max root height = "+maxRootHeight+"\n");
            return new UniformNodeHeightPrior(treeModel, maxRootHeight);
        } else {
             // the Bloomquist & Suchard variant or Welch, Rambaut & Suchard variant
            boolean useAnalytic = xo.getAttribute(ANALYTIC,true);
            boolean marginal = xo.getAttribute(MARGINAL,true);
            boolean leadingTerm = xo.getAttribute(LEADING_TERM,false);
            Logger.getLogger("dr.evomodel").info("\tUsing conditional variant with "+(useAnalytic ? "analytic" : "Monte Carlo integrated")+" expressions");
            if (useAnalytic) {
                Logger.getLogger("dr.evomodel").info("\t\tSubvariant: "+(marginal ? "marginal" : "conditional"));
                Logger.getLogger("dr.evomodel").info("\t\tApproximation: "+leadingTerm);
            }
            Logger.getLogger("dr.evomodel").info("\tPlease reference:");
            Logger.getLogger("dr.evomodel").info("\t\t (1) Welch, Rambaut and Suchard (in preparation) and");
            Logger.getLogger("dr.evomodel").info("\t\t (2) Bloomquist and Suchard (in press) Systematic Biology\n");
            if (!useAnalytic) {
//                    if( treeModel.getExternalNodeCount() > MAX_ANALYTIC_TIPS)
//                        throw new XMLParseException("Analytic evaluation of UniformNodeHeight is unreliable for > "+MAX_ANALYTIC_TIPS+" taxa");
                int mcSampleSize = xo.getAttribute(MC_SAMPLE, UniformNodeHeightPrior.DEFAULT_MC_SAMPLE);
                return new UniformNodeHeightPrior(treeModel,useAnalytic,mcSampleSize);
            }

            return new UniformNodeHeightPrior(treeModel, useAnalytic, marginal,leadingTerm);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the demographic function.";
    }

    public Class getReturnType() {
        return UniformNodeHeightPrior.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(ANALYTIC, true),
            AttributeRule.newDoubleRule(MAX_ROOT_HEIGHT, true),
            AttributeRule.newIntegerRule(MC_SAMPLE,true),
            AttributeRule.newBooleanRule(MARGINAL,true),
            AttributeRule.newBooleanRule(LEADING_TERM,true),
            new ElementRule(TreeModel.class)
    };

}
