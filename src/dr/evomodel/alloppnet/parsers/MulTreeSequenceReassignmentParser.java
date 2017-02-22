/*
 * MulTreeSequenceReassignmentParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.parsers;


import dr.evomodel.alloppnet.operators.MulTreeSequenceReassignment;
import dr.evomodel.alloppnet.speciation.MulSpeciesBindings;
import dr.evomodel.alloppnet.speciation.MulSpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */


public class MulTreeSequenceReassignmentParser extends AbstractXMLObjectParser {
	public static final String MULTREE_SEQUENCE_REASSIGNMENT = "mulTreeSequenceReassignment";
	
	
	public String getParserName() {
		return MULTREE_SEQUENCE_REASSIGNMENT;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		MulSpeciesBindings mulspb = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);
		MulSpeciesTreeModel multree = (MulSpeciesTreeModel) xo.getChild(MulSpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new MulTreeSequenceReassignment(multree, mulspb, weight);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(MulSpeciesBindings.class),
                new ElementRule(MulSpeciesTreeModel.class)
        };
	}

	@Override
	public String getParserDescription() {
		return "Operator which reassigns sequences within an allopolyploid species.";
	}

	@Override
	public Class getReturnType() {
		return MulTreeSequenceReassignment.class;
	}

}
