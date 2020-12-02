/*
 * SubtreeLeapOperatorParser.java
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

package dr.evomodelxml.bigfasttree;

import dr.evomodel.bigfasttree.CladeAwareSubtreeLeap;
import dr.evomodel.bigfasttree.CladeNodeModel;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;


/**
 */
public class CladeAwareSubtreeLeapOperatorParser extends AbstractXMLObjectParser {

    public static final String CLADE_AWARE_SUBTREE_LEAP = "cladeAwareSubtreeLeap";

    public static final String SIZE = "size";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String DISTANCE_KERNEL = "distanceKernel";
    public static final String SLIDE_ONLY = "slideOnly";

    public String getParserName() {
        return CLADE_AWARE_SUBTREE_LEAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        CladeNodeModel cladeModel =  (CladeNodeModel) xo.getChild(CladeNodeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);


        final double size = xo.getAttribute(SIZE, Double.NaN);
        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

        CladeAwareSubtreeLeap.DistanceKernelType distanceKernel = CladeAwareSubtreeLeap.DistanceKernelType.NORMAL;
//        if (xo.hasAttribute(DISTANCE_KERNEL)) {
//            try {
//                distanceKernel = CladeAwareSubtreeLeap.DistanceKernelType.valueOf(xo.getStringAttribute(DISTANCE_KERNEL).trim().toUpperCase());
//            } catch (IllegalArgumentException iae) {
//                throw new XMLParseException("Unrecognised distanceKernel attribute: " + xo.getStringAttribute(DISTANCE_KERNEL));
//            }
//        }

        if (size <= 0.0) {
            throw new XMLParseException("The SubTreeLeap size attribute must be positive and non-zero.");
        }

        if (targetAcceptance <= 0.0 || targetAcceptance >= 1.0) {
            throw new XMLParseException("Target acceptance probability has to lie in (0, 1)");
        }

        final boolean slideOnly = xo.getAttribute(SLIDE_ONLY, false);

        return new CladeAwareSubtreeLeap(cladeModel, weight, size,distanceKernel, mode, targetAcceptance,slideOnly);

    }

    public String getParserDescription() {
        return "An operator that moves a subtree a certain patristic distance.";
    }

    public Class getReturnType() {
        return SubtreeLeapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SIZE, false),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newStringRule(DISTANCE_KERNEL, true),
            AttributeRule.newBooleanRule(SLIDE_ONLY, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(CladeNodeModel.class)

    };

}