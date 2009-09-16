package dr.inference.model;

import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class MaskedParameter extends Parameter.Abstract implements VariableListener {

    public static final String MASKED_PARAMETER = "maskedParameter";
    public static final String MASKING = "mask";
    public static final String COMPLEMENT = "complement";

    public MaskedParameter(Parameter parameter) {
        this.parameter = parameter;
        this.map = new int[parameter.getDimension()];
        for(int i=0; i<map.length; i++)
            map[i] = i;
        length = map.length;
    }

    public void addMask(Parameter maskParameter, boolean ones) {
        if (maskParameter.getDimension() != parameter.getDimension())
            throw new RuntimeException("Masking parameter '"+maskParameter.getId()+"' dimension must equal base parameter '"+
                    parameter.getId() +"' dimension");
        this.maskParameter = maskParameter;
        maskParameter.addParameterListener(this);
        if (ones)
            equalValue = 1;
        else
            equalValue = 0;
        updateMask();
    }

    private void updateMask() {
        int index = 0;
        for(int i=0; i<maskParameter.getDimension(); i++) {
            // TODO Add a threshold attribute for continuous value masking
            final int maskValue = (int) maskParameter.getParameterValue(i);
            if (maskValue == equalValue) {
                map[index] = i;
                index++;
            }
        }
        length = index;
    }

    public int getDimension() {
        if (length == 0)
            throw new RuntimeException("Zero-dimensional parameter!");
            // TODO Need non-fatal mechanism to check for zero-dimensional parameters
        return length;
    }

    protected void storeValues() {
        parameter.storeParameterValues();
    }

    protected void restoreValues() {
        parameter.restoreParameterValues();
    }

    protected void acceptValues() {
       parameter.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        parameter.adoptParameterValues(source);
    }

    public double getParameterValue(int dim) {
        return parameter.getParameterValue(map[dim]);
    }

    public void setParameterValue(int dim, double value) {
        parameter.setParameterValue(map[dim],value);
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameter.setParameterValueQuietly(map[dim],value);
    }

    public String getParameterName() {
        if (getId() == null)
            return "masked" + parameter.getParameterName();
        return getId();
    }

    public void addBounds(Bounds<Double> bounds) {
        parameter.addBounds(bounds);
    }

    public Bounds<Double> getBounds() {
        return parameter.getBounds();
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        if (variable == maskParameter) {
            updateMask();
        }
        else {
            System.err.println("Called by "+variable.getId());
            throw new RuntimeException("Not yet implemented.");
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            XMLObject cxo = xo.getChild(MASKING);
            Parameter mask = (Parameter) cxo.getChild(Parameter.class);

            if (mask.getDimension() != parameter.getDimension())
                throw new XMLParseException("dim(" + parameter.getId() + ") != dim(" + mask.getId() + ")");

            MaskedParameter maskedParameter = new MaskedParameter(parameter);
            boolean ones = ! xo.getAttribute(COMPLEMENT,false);
            maskedParameter.addMask(mask, ones);

            return maskedParameter;
        }

       public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class),
                new ElementRule(MASKING,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class)
                        }),
                AttributeRule.newBooleanRule(COMPLEMENT,true),
        };

        public String getParserDescription() {
            return "A masked parameter.";
        }

        public Class getReturnType() {
            return Parameter.class;
        }

        public String getParserName() {
            return MASKED_PARAMETER;
        }
    };

    private final Parameter parameter;
    private Parameter maskParameter;
    private final int[] map;
    private int length;
    private int equalValue;
}
