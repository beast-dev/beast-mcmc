/*
 * CachedReport.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.xml.unittest;

import dr.xml.*;

public class CachedReport implements Reportable {
    private String report;
    private final Reportable reportable;

    CachedReport(Reportable reportable) {
        this.reportable = reportable;
        this.report = null;
    }

    @Override
    public String getReport() {
        if (report == null) {
            report = reportable.getReport();
        }
        return report;
    }

    private static final String CACHED_REPORT = "cachedReport";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Class getReturnType() {
            return CachedReport.class;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Reportable reportable = (Reportable) xo.getChild(Reportable.class);
            return new CachedReport(reportable);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Reportable.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Reportable object that caches the report rather than re-computing it.";
        }

        @Override
        public String getParserName() {
            return CACHED_REPORT;
        }


    };
}
