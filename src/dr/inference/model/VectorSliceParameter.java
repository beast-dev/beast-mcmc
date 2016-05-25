/*
 * VectorSliceParameter.java
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

import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class VectorSliceParameter extends CompoundParameter {

    public final static String VECTOR_SLICE_PARAMETER = "vectorSlice";
    public static final String SLICE_DIMENSION = "sliceDimension";

    private final int sliceDimension;
    private Bounds<Double> bounds;

    public Bounds<Double> getBounds() {

        if (bounds == null) {
            bounds = new sliceBounds();
        }
        return bounds;
    }

    //TODO test add bounds function

    public VectorSliceParameter(String name, int sliceDimension) {
        super(name);
        this.sliceDimension = sliceDimension;
    }

    public int getDimension() {
        return getParameterCount();
    }

    public double getParameterValue(int dim) {
        Parameter parameter = getParameter(dim);
        return parameter.getParameterValue(sliceDimension);
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void setParameterValue(int dim, double value) {
        Parameter parameter = getParameter(dim);
        parameter.setParameterValue(sliceDimension, value);
    }

    public void setParameterValueQuietly(int dim, double value) {
        Parameter parameter = getParameter(dim);
        parameter.setParameterValueQuietly(sliceDimension, value);
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        Parameter parameter = getParameter(dim);
        parameter.setParameterValueNotifyChangedAll(sliceDimension, value);
    }

    public String getDimensionName(int dim) {

        return getParameter(dim).getVariableName() + Integer.toString(sliceDimension + 1);
    }

    private class sliceBounds implements Bounds<Double>{


        @Override
        public Double getUpperLimit(int dimension) {
            return getParameter(dimension).getBounds().getUpperLimit(sliceDimension);
        }

        @Override
        public Double getLowerLimit(int dimension) {
            return getParameter(dimension).getBounds().getLowerLimit(sliceDimension);
        }

        @Override
        public int getBoundsDimension() {
            return getDimension();
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return VECTOR_SLICE_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            int sliceDimension = xo.getIntegerAttribute(SLICE_DIMENSION);
            VectorSliceParameter vectorSlice = new VectorSliceParameter(xo.getId(), sliceDimension - 1);

            for (int i = 0; i < xo.getChildCount(); i++) {
                Parameter parameter = (Parameter) xo.getChild(i);
                if (parameter instanceof MatrixParameter) {
                    MatrixParameter mp = (MatrixParameter) parameter;
                    for (int j = 0; j < mp.getParameterCount(); ++j) {
                        checkAndAdd(vectorSlice, mp.getParameter(j), sliceDimension);
                    }
                } else {
                    checkAndAdd(vectorSlice, parameter, sliceDimension);
                }
            }
            return vectorSlice;
        }

        private void checkAndAdd(final VectorSliceParameter vectorSlice, final Parameter parameter,
                                 final int sliceDimension) throws XMLParseException {
            vectorSlice.addParameter(parameter);
            if (sliceDimension < 1 || sliceDimension > parameter.getDimension()) {
                throw new XMLParseException("Slice dimension " + sliceDimension + " is invalid for a parameter" +
                        " with dimension = " + parameter.getDimension());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A vector parameter constructed from a slice of component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
                AttributeRule.newIntegerRule(SLICE_DIMENSION),
        };

        public Class getReturnType() {
            return VectorSliceParameter.class;
        }
    };


}