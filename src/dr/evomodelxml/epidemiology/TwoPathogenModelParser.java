package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.TwoPathogenModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;

public class TwoPathogenModelParser extends AbstractXMLObjectParser {

    public static final String TWO_PATHOGEN_MODEL = "twoPathogenModel";
    public static final String TRANSMISSION_RATE_ONE = "transmissionRateOne";
    public static final String TRANSMISSION_RATE_TWO = "transmissionRateTwo";
    public static final String MOVE_TO_R_RATE_ONE = "moveToRRateOne";
    public static final String MOVE_TO_R_RATE_TWO = "moveToRRateTwo";
    public static final String MOVE_TO_C_RATE_ONE = "moveToCRateOne";
    public static final String MOVE_TO_C_RATE_TWO = "moveToCRateTwo";
    public static final String SAMPLING_PROPORTION_ONE = "samplingProportionOne";
    public static final String SAMPLING_PROPORTION_TWO = "samplingProportionTwo";
    public static final String RESUSCEPTIBILITY_RATE_ONE = "resusRateOne";
    public static final String RESUSCEPTIBILITY_RATE_TWO = "resusRateTwo";
    public static final String INFECTION_RATE_MODULATION_I = "infectionRateModulationI";
    public static final String INFECTION_RATE_MODULATION_C = "infectionRateModulationC";
    public static final String RECOVERY_RATE_MODULATION = "recoveryRateModulation";
    public static final String NUM_SS = "numSS";
    public static final String NUM_SI = "numSI";
    public static final String NUM_SC = "numSC";
    public static final String NUM_SR = "numSR";

    public static final String NUM_IS = "numIS";
    public static final String NUM_II = "numII";
    public static final String NUM_IC = "numIC";
    public static final String NUM_IR = "numIR";

    public static final String NUM_CS = "numCS";
    public static final String NUM_CI = "numCI";
    public static final String NUM_CC = "numCC";
    public static final String NUM_CR = "numCR";

    public static final String NUM_RS = "numRS";
    public static final String NUM_RI = "numRI";
    public static final String NUM_RC = "numRC";
    public static final String NUM_RR = "numRR";

    public static final String ORIGIN_ONE = "originOne";
    public static final String ORIGIN_TWO = "originTwo";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String ORIGIN_TIME_NUM_SS = "originTimeNumSS";


    public String getParserName() {
        return TWO_PATHOGEN_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ArrayList<Parameter> rateParams = new ArrayList<>();

        rateParams.add((Parameter) xo.getChild(TRANSMISSION_RATE_ONE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(MOVE_TO_C_RATE_ONE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(MOVE_TO_R_RATE_ONE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(SAMPLING_PROPORTION_ONE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(RESUSCEPTIBILITY_RATE_ONE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(TRANSMISSION_RATE_TWO).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(MOVE_TO_C_RATE_TWO).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(MOVE_TO_R_RATE_TWO).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(SAMPLING_PROPORTION_TWO).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(RESUSCEPTIBILITY_RATE_TWO).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(INFECTION_RATE_MODULATION_I).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(INFECTION_RATE_MODULATION_C).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(RECOVERY_RATE_MODULATION).getChild(Parameter.class));

        final Parameter numGridPoints = (Parameter) xo.getChild(NUM_GRID_POINTS).getChild(Parameter.class);

        final Parameter cutOff = (Parameter) xo.getChild(CUT_OFF).getChild(Parameter.class);

        ArrayList<Parameter> compartmentCounts = new ArrayList<>();

        Parameter numSSParam = (Parameter) xo.getChild(NUM_SS).getChild(Parameter.class);
        if (numSSParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numSS parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numSSParam);

        Parameter numSIParam = (Parameter) xo.getChild(NUM_SI).getChild(Parameter.class);
        if (numSIParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numSI parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numSIParam);

        Parameter numSCParam = (Parameter) xo.getChild(NUM_SC).getChild(Parameter.class);
        if (numSCParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numSC parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numSCParam);

        Parameter numSRParam = (Parameter) xo.getChild(NUM_SR).getChild(Parameter.class);
        if (numSRParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numSR parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numSRParam);


        Parameter numISParam = (Parameter) xo.getChild(NUM_IS).getChild(Parameter.class);
        if (numISParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numIS parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numISParam);

        Parameter numIIParam = (Parameter) xo.getChild(NUM_II).getChild(Parameter.class);
        if (numIIParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numII parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numIIParam);

        Parameter numICParam = (Parameter) xo.getChild(NUM_IC).getChild(Parameter.class);
        if (numICParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numIC parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numICParam);

        Parameter numIRParam = (Parameter) xo.getChild(NUM_IR).getChild(Parameter.class);
        if (numIRParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numIR parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numIRParam);


        Parameter numCSParam = (Parameter) xo.getChild(NUM_CS).getChild(Parameter.class);
        if (numCSParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numCS parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numCSParam);

        Parameter numCIParam = (Parameter) xo.getChild(NUM_CI).getChild(Parameter.class);
        if (numCIParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numCI parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numCIParam);

        Parameter numCCParam = (Parameter) xo.getChild(NUM_CC).getChild(Parameter.class);
        if (numCCParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numCC parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numCCParam);

        Parameter numCRParam = (Parameter) xo.getChild(NUM_CR).getChild(Parameter.class);
        if (numCRParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numCR parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numCRParam);


        Parameter numRSParam = (Parameter) xo.getChild(NUM_RS).getChild(Parameter.class);
        if (numRSParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numRS parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numRSParam);

        Parameter numRIParam = (Parameter) xo.getChild(NUM_RI).getChild(Parameter.class);
        if (numRIParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numRI parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numRIParam);

        Parameter numRCParam = (Parameter) xo.getChild(NUM_RC).getChild(Parameter.class);
        if (numRCParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numRC parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numRCParam);

        Parameter numRRParam = (Parameter) xo.getChild(NUM_RR).getChild(Parameter.class);
        if (numRRParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numRR parameter must have dimension equal to numGridPoints");
        }
        compartmentCounts.add(numRRParam);

        // origin changes here
        Parameter originOne = (Parameter) xo.getChild(ORIGIN_ONE).getChild(Parameter.class);
        Parameter originTwo = (Parameter) xo.getChild(ORIGIN_TWO).getChild(Parameter.class);

        Parameter originTimeNumSS = (Parameter) xo.getChild(ORIGIN_TIME_NUM_SS).getChild(Parameter.class);

        //16 or 56? was originally 3
        TwoPathogenModel twoPathogenModel = new TwoPathogenModel(rateParams, compartmentCounts,
                originOne, originTwo, originTimeNumSS,16, (int)(numGridPoints.getParameterValue(0)),
                cutOff.getParameterValue(0));

        return twoPathogenModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an SIRS model";
    }

    public Class getReturnType() {
        return TwoPathogenModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TRANSMISSION_RATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOVE_TO_C_RATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOVE_TO_R_RATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(SAMPLING_PROPORTION_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(RESUSCEPTIBILITY_RATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(TRANSMISSION_RATE_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOVE_TO_C_RATE_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOVE_TO_R_RATE_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(SAMPLING_PROPORTION_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(RESUSCEPTIBILITY_RATE_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_SS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_SI,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_SC,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_SR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_IS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_II,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_IC,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_IR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_CS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_CI,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_CC,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_CR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_RS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_RI,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_RC,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_RR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(ORIGIN_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(ORIGIN_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(ORIGIN_TIME_NUM_SS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
    };

}
