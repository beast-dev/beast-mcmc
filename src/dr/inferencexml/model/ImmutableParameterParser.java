/*
 * ImmutableParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * Created by Guy Baele on 18/12/15.
 */
public class ImmutableParameterParser extends AbstractXMLObjectParser {

    public static final String IMMUTABLE_PARAMETER = "immutableParameter";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Statistic statistic = (Statistic) xo.getChild(Statistic.class);

        Parameter.Abstract immutableParameter = new Parameter.Abstract() {
            public void setParameterValueNotifyChangedAll(int dim, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void setParameterValueQuietly(int dim, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void storeValues() {
                //do nothing
            }
            public void restoreValues() {
                //do nothing
            }
            public void acceptValues() {
                //do nothing
            }

            public int getDimension() {
                return statistic.getDimension();
            }

            public void setParameterValue(int dim, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public double getParameterValue(int dim) {
                return statistic.getStatisticValue(dim);
            }
            public String getParameterName() {
                if (getId() == null)
                    return "immutable." + statistic.getStatisticName();
                return getId();
            }
            public void adoptValues(Parameter source) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void addDimension(int index, double value) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public double removeDimension(int index) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public void addBounds(Bounds<Double> bounds) {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
            public Bounds<Double> getBounds() {
                throw new RuntimeException("Forbidden call to ImmutableParameter.");
            }
        };

        return immutableParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Statistic.class),
    };

    public String getParserDescription() {
        return "An immutable parameter generated from a statistic.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return IMMUTABLE_PARAMETER;
    }

}
