/*
 * AlloppSpeciesBindingsApSpInfoParser.java
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


import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


/**
 * 
 * Parses a possibly-allopolyploid species made of individual organisms. 
 * 
 * @author Graham Jones
 *         Date: 18/04/2011
 */



/*
 * 
 * Parses a diploid or allopolyploid species, recording ploidy level
 * and running through all individual organisms belonging to the species. 
 * 
 * Part of parsing a AlloppSpeciesBindings.
 * 
 */


public class AlloppSpeciesBindingsApSpInfoParser extends
		AbstractXMLObjectParser {
	public static final String APSP = "apsp";
	public static final String PLOIDYLEVEL = "ploidylevel";

	public String getParserName() {
		return APSP;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		AlloppSpeciesBindings.Individual[] individuals = 
			new AlloppSpeciesBindings.Individual[xo.getChildCount()];
        for (int ni = 0; ni < individuals.length; ++ni) {
        	individuals[ni] = (AlloppSpeciesBindings.Individual) xo.getChild(ni);
        }
        return new AlloppSpeciesBindings.ApSpInfo(xo.getId(), xo.getIntegerAttribute(PLOIDYLEVEL), individuals);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {
			AttributeRule.newDoubleRule(PLOIDYLEVEL),
			new ElementRule(AlloppSpeciesBindings.Individual.class, 1, Integer.MAX_VALUE)
			};
	}

	@Override
	public String getParserDescription() {
		return "A diploid or allopolyploid species made of individuals";
	}

	@Override
	public Class getReturnType() {
		return AlloppSpeciesBindings.ApSpInfo.class;
	}

}
