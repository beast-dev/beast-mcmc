/*
 * PolarCoordinatesParser.java
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

package dr.evoxml;

import dr.geo.math.SphericalPolarCoordinates;
import dr.xml.*;

/**
 * @author Alexei Drummond
 *
 * @version $Id: PolarCoordinatesParser.java,v 1.1 2005/04/11 11:51:33 alexei Exp $
 */
public class PolarCoordinatesParser extends AbstractXMLObjectParser {

		public String getParserName() { return "latLong"; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			double latitude = xo.getDoubleAttribute("latitude");
			double longitude = xo.getDoubleAttribute("longitude");

			return new SphericalPolarCoordinates(latitude, longitude);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A latitude/longitude pair representing a point on the surface of the Earth.";
		}

		public Class getReturnType() { return SphericalPolarCoordinates.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newDoubleRule("longitude"),
			AttributeRule.newDoubleRule("latitude")
		};
	}

