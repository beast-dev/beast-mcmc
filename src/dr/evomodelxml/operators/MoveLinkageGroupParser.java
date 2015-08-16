/*
 * MoveLinkageGroupParser.java
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

import dr.evolution.MetagenomeData;
import dr.evomodel.operators.MoveLinkageGroup;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling (koadman)
 */

public class MoveLinkageGroupParser extends AbstractXMLObjectParser {
    public static final String MOVE_LINKAGE_GROUP = "moveLinkageGroup";

	public String getParserDescription() {
		return "Operator to reassign metagenomic reads from one linkage group to another";
	}

	public Class getReturnType() {
		return MoveLinkageGroup.class;
	}

	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final HiddenLinkageModel hlm = (HiddenLinkageModel) xo.getChild(HiddenLinkageModel.class);
        return new MoveLinkageGroup(hlm, weight);
	}

	public String getParserName() {
		return MOVE_LINKAGE_GROUP;
	}

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(MCMCOperator.WEIGHT, "Weight of the move", true),
            new ElementRule(HiddenLinkageModel.class)
    };
}
