/*
 * DnDsLoggerParser.java
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

import dr.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.evomodel.treelikelihood.utilities.DnDsLogger;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.xml.*;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "dNdSLogger";
    public static final String USE_SMOOTHING = "smooth";
    public static final String USE_DNMINUSDS = "dn-ds";
    public static final String COUNTS = "counts";
    public static final String SYNONYMOUS = "synonymous";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String[] names = DnDsLogger.traitNames;
        TreeTrait[] foundTraits = new TreeTrait[names.length];

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object obj = xo.getChild(i);
            if (obj instanceof CodonPartitionedRobustCounting) {
                CodonPartitionedRobustCounting thisCount = (CodonPartitionedRobustCounting) obj;
                for (int j = 0; j < names.length; j++) {
                    TreeTrait trait = thisCount.getTreeTrait(names[j]);
                    if (trait != null) {
                        foundTraits[j] = trait;
                    }
                }
            }
        }

        for (int i = 0; i < foundTraits.length; i++) {
            if (foundTraits[i] == null) {
                throw new XMLParseException("Unable to find trait '" + names[i] + "'");
            }
        }

        Tree tree = (Tree) xo.getChild(Tree.class);

        // Use AttributeRules for options here

        boolean useSmoothing = xo.getAttribute(USE_SMOOTHING, true);
        boolean useDnMinusDs = xo.getAttribute(USE_DNMINUSDS, false);
        boolean conditionalCounts = xo.getAttribute(COUNTS, false);
        boolean synonymous = xo.getAttribute(SYNONYMOUS, false);

        return new DnDsLogger(xo.getId(), tree, foundTraits, useSmoothing, useDnMinusDs, conditionalCounts, synonymous);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static XMLSyntaxRule[] rules = {
            new ElementRule(CodonPartitionedRobustCounting.class, 2, 2),
            new ElementRule(Tree.class),
            AttributeRule.newBooleanRule(USE_SMOOTHING, true),
            AttributeRule.newBooleanRule(COUNTS, true),
            AttributeRule.newBooleanRule(SYNONYMOUS, true),
            AttributeRule.newBooleanRule(USE_DNMINUSDS, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DnDsLogger.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
