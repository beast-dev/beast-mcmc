/*
 * AlloppNetworkPriorParser.java
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

import dr.evomodel.alloppnet.speciation.AlloppNetworkPrior;
import dr.evomodel.alloppnet.speciation.AlloppNetworkPriorModel;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.xml.*;


/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


public class AlloppNetworkPriorParser extends AbstractXMLObjectParser {

	public static final String APSPNETWORKPRIOR = "apspNetworkPrior";
	public static final String MODEL = "model";
	public static final String APSPNETWORK = "apspNetwork";


	public String getParserName() {
		return APSPNETWORKPRIOR;
	}
	
	
/*
	<apspNetworkPrior id="apspnetwork.prior" units="substitutions">
    <eventRate>
      <parameter id="apspnetwork.prior.eventRate" value="0.1" lower="0.0" upper="Infinity"/>
    </eventRate>
    <apspNetwork>
      <apspNetwork idref="apspnetwork"/>
		</apspNetwork>
	</apspNetworkPrior>
*/
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final XMLObject mxo = xo.getChild(MODEL);
		final AlloppNetworkPriorModel anpm = (AlloppNetworkPriorModel) mxo.getChild(AlloppNetworkPriorModel.class);
		final XMLObject asnmxo = xo.getChild(APSPNETWORK);
		final AlloppSpeciesNetworkModel asnm = (AlloppSpeciesNetworkModel) asnmxo.getChild(AlloppSpeciesNetworkModel.class);
		return new AlloppNetworkPrior(anpm, asnm);	
	}
	
	private  XMLSyntaxRule[] modelRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(AlloppNetworkPriorModel.class)
		};

	}

	private  XMLSyntaxRule[] asnmRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(AlloppSpeciesNetworkModel.class)
		};

	}
	
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(MODEL, modelRules()),
				new ElementRule(APSPNETWORK, asnmRules()),
				
		};
	}

	@Override
	public String getParserDescription() {
		return "Prior for an allopolyploid species network.";
	}

	@Override
	public Class getReturnType() {
		return AlloppNetworkPrior.class;
	}

}
