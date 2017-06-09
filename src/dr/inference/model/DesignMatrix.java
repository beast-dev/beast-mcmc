/*
 * DesignMatrix.java
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

package dr.inference.model;

import dr.stats.DiscreteStatistics;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc Suchard
 */
public class DesignMatrix extends MatrixParameter {
    public static final String DESIGN_MATRIX = "designMatrix";
    public static final String ADD_INTERCEPT = "addIntercept";
    public static final String FORM = "form";
    public static final String ROW_DIMENSION = "rowDimension";
    public static final String COL_DIMENSION = "colDimension";
    public static final String CHECK_IDENTIFABILITY = "checkIdentifiability";
    public static final String STANDARDIZE = "standardize";
    public static final String DYNAMIC_STANDARDIZATION = "dynamicStandardization";

    public static final String INTERCEPT = "intercept";

    public DesignMatrix(String name, boolean dynamicStandardization) {
        super(name);
        this.dynamicStandardization = dynamicStandardization;
        init();
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.variableChangedEvent(variable, index, type);
        standardizationKnown = false;
    }

    protected double getRawParameterValue(int row, int col) {
        return super.getParameterValue(row, col);
    }

    public double getParameterValue(int row, int col) {
        double value = getRawParameterValue(row, col);
        if (dynamicStandardization) {
            if (!standardizationKnown) {
                computeStandarization();
                standardizationKnown = true;
            }
            value = (value - standardizationMean[col]) / standardizationStDev[col];
        }
        return value;
    }

//    public double getParameterValue(int index) {
//        throw new RuntimeException("Univariate value from a design matrix");
//    }

    public void addParameter(Parameter param) {
        super.addParameter(param);
        clearCache();     // Changed size
    }

    public void removeParameter(Parameter param) {
        super.removeParameter(param);
        clearCache();     // Changed size
    }

    private void clearCache() {
        standardizationMean = null;
        standardizationStDev = null;

        storedStandardizationMean = null;
        storedStandardizationStDev = null;
    }

    private void computeStandarization() {
        if (standardizationMean == null) {
            standardizationMean = new double[getColumnDimension()];
        }
        if (standardizationStDev == null) {
            standardizationStDev = new double[getColumnDimension()];
        }
        for (int col = 0; col < getColumnDimension(); col++) {
            if ((getParameter(col).getId()).toLowerCase().indexOf(INTERCEPT) >= 0) {
                standardizationMean[col] = 0.0;
                standardizationStDev[col] = 1.0;
            } else {
                double[] vector = getParameter(col).getParameterValues();
                standardizationMean[col] = DiscreteStatistics.mean(vector);
                standardizationStDev[col] = Math.sqrt(DiscreteStatistics.variance(vector, standardizationMean[col]));
            }
        }
    }

    protected void storeValues() {
        super.storeValues();

        if (dynamicStandardization) {
            if (storedStandardizationMean == null) {
                storedStandardizationMean = new double[standardizationMean.length];
            }
            System.arraycopy(standardizationMean, 0, storedStandardizationMean, 0, standardizationMean.length);

            if (storedStandardizationStDev == null) {
                storedStandardizationStDev = new double[standardizationStDev.length];
            }
            System.arraycopy(standardizationStDev, 0, storedStandardizationStDev, 0, standardizationStDev.length);
        }
    }

    protected void restoreValues() {
        super.restoreValues();

        if (dynamicStandardization) {
            double[] tmp = standardizationMean;
            standardizationMean = storedStandardizationMean;
            storedStandardizationMean = tmp;

            tmp = standardizationStDev;
            standardizationStDev = storedStandardizationStDev;
            storedStandardizationStDev = tmp;
        }
    }

    public DesignMatrix(String name, Parameter[] parameters, boolean dynamicStandardization) {
        super(name, parameters);
        this.dynamicStandardization = dynamicStandardization;
        init();
    }

    private void init() {
        standardizationKnown = false;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public static void standardize(double[] vector) {
        double mean = DiscreteStatistics.mean(vector);
        double stDev = Math.sqrt(DiscreteStatistics.variance(vector, mean));
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = (vector[i] - mean) / stDev;
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


        public String getParserName() {
            return DESIGN_MATRIX;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean dynamicStandardization = xo.getAttribute(DYNAMIC_STANDARDIZATION, false);
            String name = (xo.hasId() ? xo.getId() : DESIGN_MATRIX);

            DesignMatrix designMatrix = new DesignMatrix(name, dynamicStandardization);
            boolean addIntercept = xo.getAttribute(ADD_INTERCEPT, false);
            boolean standardize = xo.getAttribute(STANDARDIZE, false);

            int dim = 0;

            if (xo.hasAttribute(FORM)) {
                String type = xo.getStringAttribute(FORM);
                if (type.compareTo("J") == 0) {
                    int rowDim = xo.getAttribute(ROW_DIMENSION, 1);
                    int colDim = xo.getAttribute(COL_DIMENSION, 1);
                    for (int i = 0; i < colDim; i++) {
                        Parameter parameter = new Parameter.Default(rowDim);
                        designMatrix.addParameter(parameter);
                    }
                } else
                    throw new XMLParseException("Unknown designMatrix form.");
            } else {

                for (int i = 0; i < xo.getChildCount(); i++) {
                    Parameter parameter = (Parameter) xo.getChild(i);
                    designMatrix.addParameter(parameter);
                    if (i == 0)
                        dim = parameter.getDimension();
                    else if (dim != parameter.getDimension())
                        throw new XMLParseException("Parameter " + (i+1) +" has dimension "+ parameter.getDimension()+ " and not "+dim+". "+
                                "All parameters must have the same dimension to construct a rectangular design matrix");
                }
            }

            if (standardize) {
                // Standardize all covariates except intercept
                for (int j = 0; j < designMatrix.getColumnDimension(); ++j) {
                    Parameter columnParameter = designMatrix.getParameter(j);
                    double[] column = columnParameter.getParameterValues();
                    standardize(column);
                    for (int i = 0; i < column.length; ++i) {
                        columnParameter.setParameterValueQuietly(i, column[i]);
                    }
                    columnParameter.setParameterValueNotifyChangedAll(0, columnParameter.getParameterValue(0));
                }
            }

            if (addIntercept) {
                Parameter intercept = new Parameter.Default(dim);
                intercept.setId(INTERCEPT);
                designMatrix.addParameter(intercept);
            }

            return designMatrix;
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
                AttributeRule.newBooleanRule(ADD_INTERCEPT, true),
                AttributeRule.newBooleanRule(CHECK_IDENTIFABILITY, true),
                new ElementRule(Parameter.class, 0, Integer.MAX_VALUE), // TODO or have the following                            
                AttributeRule.newStringRule(FORM, true),     // TODO Should have to include both FORM and DIMENSION at the same time
                AttributeRule.newIntegerRule(COL_DIMENSION, true),
                AttributeRule.newIntegerRule(ROW_DIMENSION, true),
                AttributeRule.newBooleanRule(STANDARDIZE, true),
        };

        public Class getReturnType() {
            return DesignMatrix.class;
        }
    };

    private final boolean dynamicStandardization;
    private boolean standardizationKnown = false;

    private double[] standardizationMean = null;
    private double[] standardizationStDev = null;
    private double[] storedStandardizationMean = null;
    private double[] storedStandardizationStDev = null;
}
