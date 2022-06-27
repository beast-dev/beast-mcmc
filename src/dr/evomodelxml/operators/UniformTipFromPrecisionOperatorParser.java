/*
 * UniformTipFromPrecisionOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.evomodel.operators.UniformTipFromPrecisionOperator;
import dr.xml.*;

/**
 */
public class UniformTipFromPrecisionOperatorParser extends AbstractXMLObjectParser {
    public final static String UTFP = "uniformTipFromPrecisionOperator";
    public static final String LOWER = "lower";
    public static final String UPPER = "upper";

    public String getParserName() {
        return UTFP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        if( parameter.getDimension() == 0 ) {
            throw new XMLParseException("parameter with 0 dimension.");
        }

        Double lower = null;
        Double upper = null;

        if (xo.hasAttribute(LOWER)) {
            lower = xo.getDoubleAttribute(LOWER);
        }

        if (xo.hasAttribute(UPPER)) {
            upper = xo.getDoubleAttribute(UPPER);
        }

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Taxon taxon = (Taxon) xo.getChild(Taxon.class);

        return new UniformTipFromPrecisionOperator(parameter, weight, taxon, tree, lower, upper);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "An operator that picks new tip age values uniformly at random (constrained by the tree).";
    }

    public Class getReturnType() {
        return UniformTipFromPrecisionOperator.class;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(LOWER, true),
            AttributeRule.newDoubleRule(UPPER, true),
            new ElementRule(Parameter.class),
            new ElementRule(Taxon.class),
            new ElementRule(TreeModel.class)
    };
}
