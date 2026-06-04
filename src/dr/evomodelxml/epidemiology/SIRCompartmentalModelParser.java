package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.SIRCompartmentalModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;

public class SIRCompartmentalModelParser extends AbstractXMLObjectParser {

    public static final String SIR_COMPARTMENTAL_MODEL = "sirCompartmentalModel";
    public static final String TRANSMISSION_RATE = "transmissionRate";
    public static final String RECOVERY_RATE = "recoveryRate";
    public static final String SAMPLING_PROPORTION = "samplingProportion";
    public static final String RESUSCEPTIBILITY_RATE = "resusRate";
    public static final String NUM_S = "numS";
    public static final String NUM_I = "numI";
    public static final String NUM_R = "numR";
    public static final String ORIGIN = "origin";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String ORIGIN_TIME_NUM_S = "originTimeNumS";

    public String getParserName() {
        return SIR_COMPARTMENTAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ArrayList<Parameter> rateParams = new ArrayList<>();

        rateParams.add((Parameter) xo.getChild(TRANSMISSION_RATE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(RECOVERY_RATE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(SAMPLING_PROPORTION).getChild(Parameter.class));

        // Change this later to make it optional?
        rateParams.add((Parameter) xo.getChild(RESUSCEPTIBILITY_RATE).getChild(Parameter.class));

        final Parameter numGridPoints = (Parameter) xo.getChild(NUM_GRID_POINTS).getChild(Parameter.class);

        final Parameter cutOff = (Parameter) xo.getChild(CUT_OFF).getChild(Parameter.class);

        ArrayList<Parameter> compartmentCounts = new ArrayList<>();

        Parameter numSParam = (Parameter) xo.getChild(NUM_S).getChild(Parameter.class);

        if (numSParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numS parameter must have dimension equal to numGridPoints");
        }

        compartmentCounts.add(numSParam);

        Parameter numIParam = (Parameter) xo.getChild(NUM_I).getChild(Parameter.class);

        if (numIParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numI parameter must have dimension equal to numGridPoints");
        }

        compartmentCounts.add(numIParam);

        Parameter numRParam = (Parameter) xo.getChild(NUM_R).getChild(Parameter.class);

        if (numRParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numR parameter must have dimension equal to numGridPoints");
        }

        compartmentCounts.add(numRParam);

        Parameter origin = (Parameter) xo.getChild(ORIGIN).getChild(Parameter.class);

        Parameter originTimeNumS = (Parameter) xo.getChild(ORIGIN_TIME_NUM_S).getChild(Parameter.class);

        SIRCompartmentalModel sirModel = new SIRCompartmentalModel(rateParams, compartmentCounts,
                origin, originTimeNumS,3, (int)(numGridPoints.getParameterValue(0)),
                cutOff.getParameterValue(0));

        return sirModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an SIRS model";
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
            new ElementRule(RESUSCEPTIBILITY_RATE,
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
            new ElementRule(ORIGIN_TIME_NUM_S,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
    };

}
