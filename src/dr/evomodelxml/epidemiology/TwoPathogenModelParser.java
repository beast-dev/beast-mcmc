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
    public static final String INFECTION_RATE_MODULATION_R = "infectionRateModulationR";
    public static final String BIRTH_MORTALITY_RATE = "birthMortalityRate";
    public static final String IMPORTATION_RATE_ONE = "importationRateOne";
    public static final String IMPORTATION_RATE_TWO = "importationRateTwo";
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
    public static final String ORIGIN_TIME_NUM_SI = "originTimeNumSI";
    public static final String ORIGIN_TIME_NUM_IS = "originTimeNumIS";

    public static final String SEASONAL_MODEL = "seasonalModel";
    public static final String SEASONAL_PERIOD = "seasonalPeriod";
    public static final String SEASONAL_AMP_ONE = "seasonalAmplitudeOne";
    public static final String MOST_RECENT_SAMPLING_DATE_ONE = "mostRecentSamplingDateOne";
    public static final String SEASONAL_AMP_TWO = "seasonalAmplitudeTwo";
    public static final String SEASONAL_PEAK_ONE = "seasonalPeakDayOne";
    public static final String SEASONAL_PEAK_TWO = "seasonalPeakDayTwo";
    public static final String MOST_RECENT_SAMPLING_DATE_TWO = "mostRecentSamplingDateTwo";


    public String getParserName() {
        return TWO_PATHOGEN_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean seasonal = xo.getBooleanAttribute(SEASONAL_MODEL, false);

        Parameter seasonalPeriod = null;
        if (xo.getChild(SEASONAL_PERIOD) != null) {
            seasonalPeriod = (Parameter) xo.getChild(SEASONAL_PERIOD).getChild(Parameter.class);
        }
        if (seasonal && seasonalPeriod == null) {
            throw new RuntimeException("seasonal model requires seasonalPeriod");
        }

        Parameter seasonalAmplitudeOne = null;
        if(xo.getChild(SEASONAL_AMP_ONE) != null) {
            seasonalAmplitudeOne = (Parameter) xo.getChild(SEASONAL_AMP_ONE).getChild(Parameter.class);
        }
        if (seasonal && seasonalAmplitudeOne == null) {
            throw new RuntimeException("seasonal model requires seasonalAmplitudeOne");
        }

        Parameter seasonalAmplitudeTwo = null;
        if(xo.getChild(SEASONAL_AMP_TWO) != null) {
            seasonalAmplitudeTwo = (Parameter) xo.getChild(SEASONAL_AMP_TWO).getChild(Parameter.class);
        }
        if (seasonal && seasonalAmplitudeTwo == null) {
            throw new RuntimeException("seasonal model requires seasonalAmplitudeTwo");
        }

        Parameter seasonalPeakDayOne = null;
        if(xo.getChild(SEASONAL_PEAK_ONE) != null) {
            seasonalPeakDayOne = (Parameter) xo.getChild(SEASONAL_PEAK_ONE).getChild(Parameter.class);
        }
        if (seasonal && seasonalPeakDayOne == null) {
            throw new RuntimeException("seasonal model requires seasonalPeakDayOne");
        }

        Parameter seasonalPeakDayTwo = null;
        if(xo.getChild(SEASONAL_PEAK_TWO) != null) {
            seasonalPeakDayTwo = (Parameter) xo.getChild(SEASONAL_PEAK_TWO).getChild(Parameter.class);
        }
        if (seasonal && seasonalPeakDayTwo == null) {
            throw new RuntimeException("seasonal model requires seasonalPeakDayTwo");
        }

        Parameter mostRecentSamplingDateOne = null;
        if(xo.getChild(MOST_RECENT_SAMPLING_DATE_ONE) != null) {
            mostRecentSamplingDateOne = (Parameter) xo.getChild(MOST_RECENT_SAMPLING_DATE_ONE).getChild(Parameter.class);
        }

        Parameter mostRecentSamplingDateTwo = null;
        if(xo.getChild(MOST_RECENT_SAMPLING_DATE_TWO) != null) {
            mostRecentSamplingDateTwo = (Parameter) xo.getChild(MOST_RECENT_SAMPLING_DATE_TWO).getChild(Parameter.class);
        }

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

        rateParams.add((Parameter) xo.getChild(INFECTION_RATE_MODULATION_R).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(BIRTH_MORTALITY_RATE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(IMPORTATION_RATE_ONE).getChild(Parameter.class));

        rateParams.add((Parameter) xo.getChild(IMPORTATION_RATE_TWO).getChild(Parameter.class));

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

        double originTimeNumSI = 1.0;
        if(xo.getChild(ORIGIN_TIME_NUM_SI)!=null){
            originTimeNumSI = ((Parameter) xo.getChild(ORIGIN_TIME_NUM_SI).getChild(Parameter.class)).getParameterValue(0);
        }

        double originTimeNumIS = 1.0;
        if(xo.getChild(ORIGIN_TIME_NUM_IS)!=null){
            originTimeNumIS = ((Parameter) xo.getChild(ORIGIN_TIME_NUM_IS).getChild(Parameter.class)).getParameterValue(0);
        }

        TwoPathogenModel twoPathogenModel = new TwoPathogenModel(rateParams, compartmentCounts,
                originOne, originTwo, originTimeNumSS, originTimeNumSI,originTimeNumIS, 79, (int)(numGridPoints.getParameterValue(0)),
                cutOff.getParameterValue(0),seasonal, seasonalPeriod, seasonalAmplitudeOne, seasonalPeakDayOne,
                mostRecentSamplingDateOne, seasonalAmplitudeTwo, seasonalPeakDayTwo, mostRecentSamplingDateTwo);

        return twoPathogenModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an two-pathogen compartmental model (Shrestha et al., 2011)";
    }

    public Class getReturnType() {
        return TwoPathogenModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SEASONAL_MODEL, true),
            new ElementRule(TRANSMISSION_RATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(SEASONAL_PERIOD,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(SEASONAL_AMP_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(SEASONAL_AMP_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(SEASONAL_PEAK_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(SEASONAL_PEAK_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(MOST_RECENT_SAMPLING_DATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOST_RECENT_SAMPLING_DATE_TWO,
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
            new ElementRule(INFECTION_RATE_MODULATION_I,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(INFECTION_RATE_MODULATION_C,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(INFECTION_RATE_MODULATION_R,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(BIRTH_MORTALITY_RATE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(IMPORTATION_RATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(IMPORTATION_RATE_TWO,
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
            new ElementRule(ORIGIN_TIME_NUM_SI,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(ORIGIN_TIME_NUM_IS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(NUM_GRID_POINTS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(CUT_OFF,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
    };

}
