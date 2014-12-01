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
                vectorSlice.addParameter(parameter);
                if (sliceDimension < 1 || sliceDimension > parameter.getDimension()) {
                    throw new XMLParseException("Slice dimension " + sliceDimension + " is invalid for a parameter" +
                    " with dimension = " + parameter.getDimension());
                }
            }
            return vectorSlice;
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