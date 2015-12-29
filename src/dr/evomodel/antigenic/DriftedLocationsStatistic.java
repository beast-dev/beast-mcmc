/*
 * DriftedLocationsStatistic.java
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

package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.xml.*;

/**
 * @author Trevor Bedford
 * @author Marc A. Suchard
 */
public class DriftedLocationsStatistic extends Statistic.Abstract implements VariableListener {

    public static final String DRIFTED_LOCATIONS_STATISTIC = "driftedLocationsStatistic";

    public DriftedLocationsStatistic(MatrixParameter locationsParameter, Parameter offsetsParameter, Parameter locationDriftParameter) {
        this.locationsParameter = locationsParameter;
        locationsParameter.addParameterListener(this);
        this.offsetsParameter = offsetsParameter;
        offsetsParameter.addParameterListener(this);
        this.locationDriftParameter = locationDriftParameter;
        locationDriftParameter.addParameterListener(this);
    }

    public int getDimension() {
        return locationsParameter.getDimension();
    }

    // strain index
    public int getColumnIndex(int dim) {
        return dim / locationsParameter.getRowDimension();
    }

    // dimension index
    public int getRowIndex(int dim) {
        int x = getColumnIndex(dim);
        return dim - x * locationsParameter.getRowDimension();
    }

    public double getStatisticValue(int dim) {

        double val;

        // x is location count, y is location dimension
        int x = getColumnIndex(dim);
        int y = getRowIndex(dim);
        Parameter loc = locationsParameter.getParameter(x);

        if (y == 0) {
            val = loc.getParameterValue(y) + locationDriftParameter.getParameterValue(0) * offsetsParameter.getParameterValue(x);
        }
        else {
            val = loc.getParameterValue(y);
        }


        return val;

    }

    public String getDimensionName(int dim) {
        return locationsParameter.getDimensionName(dim);
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String LOCATIONS = "locations";
        public final static String OFFSETS = "offsets";
        public final static String LOCATION_DRIFT = "locationDrift";

        public String getParserName() {
            return DRIFTED_LOCATIONS_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter locations = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);
            Parameter offset = (Parameter) xo.getElementFirstChild(OFFSETS);
            Parameter locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
            return new DriftedLocationsStatistic(locations, offset, locationDrift);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return DriftedLocationsStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LOCATIONS, MatrixParameter.class),
            new ElementRule(OFFSETS, Parameter.class),
            new ElementRule(LOCATION_DRIFT, Parameter.class)
        };
    };

    private MatrixParameter locationsParameter;
    private Parameter offsetsParameter;
    private Parameter locationDriftParameter;
}
