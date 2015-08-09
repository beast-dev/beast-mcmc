/*
 * XMLUnits.java
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

package dr.evoxml.util;

import dr.evolution.util.Units;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Alexei Drummond
 * @version $Id: XMLUnits.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public interface XMLUnits extends Units {

    final static String GENERATIONS = "generations";
    final static String DAYS = "days";
    final static String MONTHS = "months";
    final static String YEARS = "years";
    public final static String SUBSTITUTIONS = "substitutions";
    // Mutations has been replaced with substitutions...
    final static String MUTATIONS = "mutations";
    final static String UNKNOWN = "unknown";

    public final static String UNITS = "units";

    XMLSyntaxRule UNITS_RULE = new StringAttributeRule("units", "the units", UNIT_NAMES, false);
    XMLSyntaxRule[] SYNTAX_RULES = {UNITS_RULE};

    class Utils {

        public static Units.Type getUnitsAttr(XMLObject
                xo) throws XMLParseException {

            Units.Type units = dr.evolution.util.Units.Type.GENERATIONS;
            if (xo.hasAttribute(UNITS)) {
                String unitsAttr = (String) xo.getAttribute(UNITS);
                if (unitsAttr.equals(YEARS)) {
                    units = dr.evolution.util.Units.Type.YEARS;
                } else if (unitsAttr.equals(MONTHS)) {
                    units = dr.evolution.util.Units.Type.MONTHS;
                } else if (unitsAttr.equals(DAYS)) {
                    units = dr.evolution.util.Units.Type.DAYS;
                } else if (unitsAttr.equals(SUBSTITUTIONS) || unitsAttr.equals(MUTATIONS)) {
                    units = dr.evolution.util.Units.Type.SUBSTITUTIONS;
                }
            }
            return units;
        }

    }
}
