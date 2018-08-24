/*
 * OriginDestinationDesignMatrix.java
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

package dr.inference.model;

import dr.evolution.datatype.DataType;
import dr.geo.GreatCircleDistances;
import dr.geo.math.SphericalPolarCoordinates;
import dr.stats.DiscreteStatistics;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 */
public class OriginDestinationDesignMatrix extends DesignMatrix {
    public static final String DESIGN_MATRIX = "originDestinationDesignMatrix";
    public static final String LAT_TRAIT = "latitudeTrait";
    public static final String LONG_TRAIT = "longitudeTrait";
//    public static final String ADD_INTERCEPT = "addIntercept";
//    public static final String FORM = "form";
//    public static final String ROW_DIMENSION = "rowDimension";
//    public static final String COL_DIMENSION = "colDimension";
//    public static final String CHECK_IDENTIFABILITY = "checkIdentifiability";
//    public static final String STANDARDIZE = "standardize";
    public static final String DYNAMIC_STANDARDIZATION = "dynamicStandardization";

//    public static final String INTERCEPT = "intercept";

    final private DesignMatrix baseMatrix;
    final DataType dataType;
//    final int dim;


    final int[] originIndex;
    final int[] destinationIndex;

    abstract class Covariate {
        final protected Parameter parameter;
//        final int index;

//        protected Covariate(Parameter parameter) {
//            this.parameter = parameter;
//        }

        protected Covariate(Parameter parameter) {
            this.parameter = parameter;
        }

        abstract double getValue(int row);
    }

    class Origin extends Covariate {

        protected Origin(Parameter parameter) {
            super(parameter);
        }

        @Override
        double getValue(int row) {
            return parameter.getParameterValue(originIndex[row]);
        }
    }

    class Destination extends Covariate {

        protected Destination(Parameter parameter) {
            super(parameter);
        }

        @Override
        double getValue(int row) {
            return parameter.getParameterValue(destinationIndex[row]);
        }
    }

    class Distance extends Covariate {

        final private Parameter latParameter;
        final private Parameter longParameter;

        protected Distance(Parameter latParameter, Parameter longParameter) {
            super(latParameter);
            this.latParameter = latParameter;
            this.longParameter = longParameter;
        }

        @Override
        double getValue(int row) {

            double originLat = latParameter.getParameterValue(originIndex[row]);
            double originLong = longParameter.getParameterValue(originIndex[row]);

            double destinationLat = latParameter.getParameterValue(destinationIndex[row]);
            double destinationLong = longParameter.getParameterValue(destinationIndex[row]);

            SphericalPolarCoordinates origin = new SphericalPolarCoordinates(originLat, originLong);
            SphericalPolarCoordinates destination = new SphericalPolarCoordinates(destinationLat, destinationLong);
            return origin.distance(destination);
        }
    }

    final private List<Covariate> covariates = new ArrayList<Covariate>();

    final int rowDimension;

    public OriginDestinationDesignMatrix(String name,
                                         DesignMatrix baseMatrix,
                                         DataType dataType,
                                         boolean dynamicStandardization,
                                         Parameter latParameter,
                                         Parameter longParameter) {
        super(name, dynamicStandardization);

        this.dataType = dataType;
        this.baseMatrix = baseMatrix;

        System.err.println("Iterating over " + baseMatrix.getParameters().size() + " parameters");

        System.err.println("unique count: " + baseMatrix.getUniqueParameterCount());
//        System.exit(-1);

//        for (int i = 0; i < baseMatrix.getP)

        for (int i = 0; i < baseMatrix.getUniqueParameterCount(); ++i) {
            Parameter param = baseMatrix.getUniqueParameter(i);
//            System.err.println("Name: " + param.getId());
            super.addParameter(param);
            if (param != latParameter && param != longParameter) {
                covariates.add(new Origin(param));
                covariates.add(new Destination(param));
            }
        }

        if (latParameter != null && longParameter != null) {
            covariates.add(new Distance(latParameter, longParameter));
        }

        this.rowDimension = dataType.getStateCount() * (dataType.getStateCount() - 1);

        originIndex = new int[rowDimension];
        destinationIndex = new int[rowDimension];

        int index = 0;
        for (int origin = 0; origin < dataType.getStateCount(); ++origin) {
            for (int destination = origin + 1; destination < dataType.getStateCount(); ++destination) {
                originIndex[index] = origin;
                destinationIndex[index] = destination;
                ++index;
            }
        }

        for (int destination = 0; destination < dataType.getStateCount(); ++destination) {
            for (int origin = destination + 1; origin < dataType.getStateCount(); ++origin) {
                originIndex[index] = origin;
                destinationIndex[index] = destination;
                ++index;
            }
        }

        System.err.println("#row = " + baseMatrix.getRowDimension());
        System.err.println("#col = " + baseMatrix.getColumnDimension());

        System.err.println("#row = " + getRowDimension());
        System.err.println("#col = " + getColumnDimension());
        System.exit(-1);
    }

    public int getColumnDimension() {
        return covariates.size();
    }

    public int getRowDimension() {
        return rowDimension;
    }

    public void addParameter(Parameter param) {
        throw new RuntimeException("Should not be called");
    }

    public void removeParameter(Parameter param) {
        throw new RuntimeException("Should not be called");
    }

    protected double getRawParameterValue(int row, int col) {

        double value = covariates.get(col).getValue(row);
        return value;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


        private Parameter getMatchingParameter(DesignMatrix designMatrix, String name) {
            for (int i = 0; i < designMatrix.getUniqueParameterCount(); ++i) {
                Parameter param = designMatrix.getUniqueParameter(i);
                if (param.getId().compareTo(name) == 0) {
                    return param;
                }
            }
            return null;
        }

        public String getParserName() {
            return DESIGN_MATRIX;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean dynamicStandardization = xo.getAttribute(DYNAMIC_STANDARDIZATION, false);
            String name = (xo.hasId() ? xo.getId() : DESIGN_MATRIX);

            final DataType dataType = (DataType) xo.getChild(DataType.class);
            final DesignMatrix designMatrix = (DesignMatrix) xo.getChild(DesignMatrix.class);

            Parameter latParameter = null;
            Parameter longParameter = null;

            if (xo.hasAttribute(LAT_TRAIT)) {
                if (xo.hasAttribute(LONG_TRAIT)) {
                    String latName = xo.getStringAttribute(LAT_TRAIT);
                    String longName = xo.getStringAttribute(LONG_TRAIT);

                    latParameter = getMatchingParameter(designMatrix, latName);
                    longParameter = getMatchingParameter(designMatrix, longName);

                    if (latParameter == null) {
                        throw new XMLParseException("Unable to find trait named `" + latName + "`");
                    }

                    if (longParameter == null) {
                        throw new XMLParseException("Unable to find trait named `" + longName + "`");
                    }
                    
                } else {
                    throw new XMLParseException("Must provide both latitude and longitude trait names");
                }
            }

            OriginDestinationDesignMatrix matrix = new OriginDestinationDesignMatrix(name,
                    designMatrix, dataType, dynamicStandardization,
                    latParameter, longParameter);

            return matrix;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
//                AttributeRule.newBooleanRule(STANDARDIZE, true),
                new ElementRule(DesignMatrix.class),
                new ElementRule(DataType.class),
                AttributeRule.newStringRule(LAT_TRAIT, true),
                AttributeRule.newStringRule(LONG_TRAIT, true),
        };

        public Class getReturnType() {
            return OriginDestinationDesignMatrix.class;
        }
    };

//    private final boolean dynamicStandardization;
//    private boolean standardizationKnown = false;
//
//    private double[] standardizationMean = null;
//    private double[] standardizationStDev = null;
//    private double[] storedStandardizationMean = null;
//    private double[] storedStandardizationStDev = null;
}
