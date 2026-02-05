package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.SIRCompartmentalModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class SIRCompartmentalModelParser extends AbstractXMLObjectParser {

    public static final String SIR_COMPARTMENTAL_MODEL = "sirCompartmentalModel";
    public static final String TRANSMISSION_RATE = "transmissionRate";
    public static final String RECOVERY_RATE = "recoveryRate";
    public static final String SAMPLING_PROPORTION = "samplingProportion";
    public static final String NUM_S = "numS";
    public static final String NUM_I = "numI";
    public static final String NUM_R = "numR";
    public static final String ORIGIN = "origin";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";

    public String getParserName() {
        return SIR_COMPARTMENTAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter transmissionRate = (Parameter) xo.getChild(TRANSMISSION_RATE).getChild(Parameter.class);
        Parameter recoveryRate = (Parameter) xo.getChild(RECOVERY_RATE).getChild(Parameter.class);
        Parameter samplingRate = (Parameter) xo.getChild(SAMPLING_PROPORTION).getChild(Parameter.class);
        Parameter numS = (Parameter) xo.getChild(NUM_S).getChild(Parameter.class);
        Parameter numI = (Parameter) xo.getChild(NUM_I).getChild(Parameter.class);
        Parameter numR = (Parameter) xo.getChild(NUM_R).getChild(Parameter.class);
        Parameter origin = (Parameter) xo.getChild(ORIGIN).getChild(Parameter.class);
        final Parameter numGridPoints = xo.hasChildNamed(NUM_GRID_POINTS) ? (Parameter) xo.getElementFirstChild(NUM_GRID_POINTS): new Parameter.Default(1.0);
        final Parameter cutOff = xo.hasChildNamed(CUT_OFF) ? (Parameter) xo.getElementFirstChild(CUT_OFF) : new Parameter.Default(Double.POSITIVE_INFINITY);

        SIRCompartmentalModel sirModel = new SIRCompartmentalModel(transmissionRate, recoveryRate, samplingRate,
                numS, numI, numR, origin, (int)(numGridPoints.getParameterValue(0)),
                cutOff.getParameterValue(0));

        return sirModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an SIR model";
    }

    public Class getReturnType() {
        return SIRCompartmentalModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TRANSMISSION_RATE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(RECOVERY_RATE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(SAMPLING_PROPORTION,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_S,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_I,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_R,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(ORIGIN,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
    };

}
