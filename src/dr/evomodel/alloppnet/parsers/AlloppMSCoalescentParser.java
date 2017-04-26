/*
 * AlloppMSCoalescentParser.java
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

import dr.evomodel.alloppnet.speciation.AlloppMSCoalescent;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */

public class AlloppMSCoalescentParser extends AbstractXMLObjectParser {

    public static final String ALLOPPMSCOALESCENT = "apspCoalescent";


    public String getParserName() {
	    return ALLOPPMSCOALESCENT;
	}



	
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
        final AlloppSpeciesNetworkModel apspnetwork = 
        	       (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);
        return new AlloppMSCoalescent(apsp, apspnetwork);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class),
        };
    }

	

	@Override
	public String getParserDescription() {
		return "Likelihood of a set of gene trees embedded in a allopolyploid species network.";
	}

	@Override
	public Class getReturnType() {
		return AlloppMSCoalescent.class;
	}

}
