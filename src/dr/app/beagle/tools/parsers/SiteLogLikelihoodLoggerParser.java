/*
 * SiteLogLikelihoodLoggerParser.java
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

package dr.app.beagle.tools.parsers;

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.tools.SiteLogLikelihoodLogger;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class SiteLogLikelihoodLoggerParser extends AbstractXMLObjectParser {

	public static final String SITE_LOGLIKELIHOOD_LOGGER = "siteLogLikelihood";

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		SiteLogLikelihoodLogger siteLogLikelihoodLogger;
		BeagleTreeLikelihood beagleTreeLikelihood = null;

		for (int i = 0; i < xo.getChildCount(); i++) {
			beagleTreeLikelihood = (BeagleTreeLikelihood) xo.getChild(i);
		}

		siteLogLikelihoodLogger = new SiteLogLikelihoodLogger(
				beagleTreeLikelihood);

		return siteLogLikelihoodLogger;
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] { new ElementRule(BeagleTreeLikelihood.class) };
	}// END: getSyntaxRules

	@Override
	public String getParserName() {
		return SITE_LOGLIKELIHOOD_LOGGER;
	}// END: getParserName

	@Override
	public String getParserDescription() {
		return "Beagle site logLikelihood";
	}// END: getParserDescription

	@Override
	public Class<SiteLogLikelihoodLogger> getReturnType() {
		return SiteLogLikelihoodLogger.class;
	}// getReturnType

}// END: class
