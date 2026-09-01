package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.TwoPathogenModelBirthRates;
import dr.inference.model.Parameter;
import dr.xml.*;

public class TwoPathogenModelBirthRatesParser extends AbstractXMLObjectParser {

    public static final String TWO_PATHOGEN_MODEL_BIRTH_RATES = "twoPathogenModelBirthRates";
    public static final String PATHOGEN_NUMBER = "pathogenNumber";
    public static final String TRANSMISSION_RATE = "transmissionRate";
    public static final String NUM_SS = "numSS";
    public static final String NUM_SI = "numSI";
    public static final String NUM_SC = "numSC";
    public static final String NUM_SR = "numSR";
    public static final String NUM_IS = "numIS";
    public static final String NUM_CS = "numCS";
    public static final String NUM_RS = "numRS";
    public static final String INFECTION_RATE_MODULATION_I = "infectionRateModulationI";
    public static final String INFECTION_RATE_MODULATION_C = "infectionRateModulationC";
    public static final String INFECTION_RATE_MODULATION_R = "infectionRateModulationR";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String ORIGIN_TIME_NUM_SS = "originTimeNumSS";
    public static final String SEASONAL_MODEL = "seasonalModel";
    public static final String SEASONAL_PERIOD = "seasonalPeriod";
    public static final String SEASONAL_AMP = "seasonalAmplitude";
    public static final String SEASONAL_PEAK = "seasonalPeakDay";
    public static final String MOST_RECENT_SAMPLING_DATE_ONE = "mostRecentSamplingDateOne";
    public static final String MOST_RECENT_SAMPLING_DATE_TWO = "mostRecentSamplingDateTwo";


    public String getParserName() {
        return TWO_PATHOGEN_MODEL_BIRTH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean seasonal = xo.getBooleanAttribute(SEASONAL_MODEL, false);

        Parameter pathogenNumber = (Parameter) xo.getChild(PATHOGEN_NUMBER).getChild(Parameter.class);
        int pathNumInt = (int) pathogenNumber.getParameterValue(0);

        Parameter seasonalPeriod = null;
        if (xo.getChild(SEASONAL_PERIOD) != null) {
            seasonalPeriod = (Parameter) xo.getChild(SEASONAL_PERIOD).getChild(Parameter.class);
        }
        if (seasonal && seasonalPeriod == null) {
            throw new RuntimeException("seasonal model requires seasonalPeriod");
        }

        Parameter seasonalAmplitude = null;
        if(xo.getChild(SEASONAL_AMP) != null) {
            seasonalAmplitude = (Parameter) xo.getChild(SEASONAL_AMP).getChild(Parameter.class);
        }
        if (seasonal && seasonalAmplitude == null) {
            throw new RuntimeException("seasonal model requires seasonalAmplitude");
        }

        Parameter seasonalPeakDay = null;
        if(xo.getChild(SEASONAL_PEAK) != null) {
            seasonalPeakDay = (Parameter) xo.getChild(SEASONAL_PEAK).getChild(Parameter.class);
        }
        if (seasonal && seasonalPeakDay == null) {
            throw new RuntimeException("seasonal model requires seasonalPeakDay");
        }

        Parameter mostRecentSamplingDateOne = (Parameter) xo.getChild(MOST_RECENT_SAMPLING_DATE_ONE).getChild(Parameter.class);

        Parameter mostRecentSamplingDateTwo = (Parameter) xo.getChild(MOST_RECENT_SAMPLING_DATE_TWO).getChild(Parameter.class);

        final Parameter numGridPoints = (Parameter) xo.getChild(NUM_GRID_POINTS).getChild(Parameter.class);

        final Parameter cutOff = (Parameter) xo.getChild(CUT_OFF).getChild(Parameter.class);

        Parameter numSSParam = (Parameter) xo.getChild(NUM_SS).getChild(Parameter.class);
        if (numSSParam.getDimension() != numGridPoints.getParameterValue(0)) {
            throw new RuntimeException("numSS parameter must have dimension equal to numGridPoints");
        }

        Parameter numSIParam = null;
        if(xo.getChild(NUM_SI) != null) {
            numSIParam = (Parameter) xo.getChild(NUM_SI).getChild(Parameter.class);
            if (numSIParam.getDimension() != numGridPoints.getParameterValue(0)) {
                throw new RuntimeException("numSI parameter must have dimension equal to numGridPoints");
            }
        }
        if(pathNumInt == 1 && numSIParam == null){
            throw new RuntimeException("numSI parameter must not be null for pathogen 1");
        }

        Parameter numSCParam = null;
        if(xo.getChild(NUM_SC) != null) {
            numSCParam = (Parameter) xo.getChild(NUM_SC).getChild(Parameter.class);
            if (numSCParam.getDimension() != numGridPoints.getParameterValue(0)) {
                throw new RuntimeException("numSC parameter must have dimension equal to numGridPoints");
            }
        }
        if(pathNumInt == 1 && numSCParam == null){
            throw new RuntimeException("numSC parameter must not be null for pathogen 1");
        }

        Parameter numSRParam = null;
        if(xo.getChild(NUM_SR) != null) {
            numSRParam = (Parameter) xo.getChild(NUM_SR).getChild(Parameter.class);
            if (numSRParam.getDimension() != numGridPoints.getParameterValue(0)) {
                throw new RuntimeException("numSR parameter must have dimension equal to numGridPoints");
            }
        }
        if(pathNumInt == 1 && numSRParam == null){
            throw new RuntimeException("numSR parameter must not be null for pathogen 1");
        }

        Parameter numISParam = null;
        if(xo.getChild(NUM_IS) != null) {
            numISParam = (Parameter) xo.getChild(NUM_IS).getChild(Parameter.class);
            if (numISParam.getDimension() != numGridPoints.getParameterValue(0)) {
                throw new RuntimeException("numIS parameter must have dimension equal to numGridPoints");
            }
        }
        if(pathNumInt == 2 && numISParam == null){
            throw new RuntimeException("numIS parameter must not be null for pathogen 2");
        }

        Parameter numCSParam = null;
        if(xo.getChild(NUM_CS) != null) {
            numCSParam = (Parameter) xo.getChild(NUM_CS).getChild(Parameter.class);
            if (numCSParam.getDimension() != numGridPoints.getParameterValue(0)) {
                throw new RuntimeException("numCS parameter must have dimension equal to numGridPoints");
            }
        }
        if(pathNumInt == 2 && numCSParam == null){
            throw new RuntimeException("numCS parameter must not be null for pathogen 2");
        }

        Parameter numRSParam = null;
        if(xo.getChild(NUM_RS) != null) {
            numRSParam = (Parameter) xo.getChild(NUM_RS).getChild(Parameter.class);
            if (numRSParam.getDimension() != numGridPoints.getParameterValue(0)) {
                throw new RuntimeException("numRS parameter must have dimension equal to numGridPoints");
            }
        }
        if(pathNumInt == 2 && numRSParam == null){
            throw new RuntimeException("numRS parameter must not be null for pathogen 2");
        }


        Parameter originTimeNumSS = (Parameter) xo.getChild(ORIGIN_TIME_NUM_SS).getChild(Parameter.class);

        Parameter tRate = (Parameter) xo.getChild(TRANSMISSION_RATE).getChild(Parameter.class);

        Parameter iRateModI = (Parameter) xo.getChild(INFECTION_RATE_MODULATION_I).getChild(Parameter.class);

        Parameter iRateModC = (Parameter) xo.getChild(INFECTION_RATE_MODULATION_C).getChild(Parameter.class);

        Parameter iRateModR = (Parameter) xo.getChild(INFECTION_RATE_MODULATION_R).getChild(Parameter.class);

        return new TwoPathogenModelBirthRates(pathNumInt, seasonal, tRate, numSSParam,
                numSIParam, numSCParam, numSRParam, numISParam, numCSParam, numRSParam, iRateModI, iRateModC,
                iRateModR, originTimeNumSS, seasonalAmplitude, seasonalPeakDay, seasonalPeriod, cutOff,
                numGridPoints, mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element translates parameters from a two-pathogen compartmental model (Shrestha et al., 2011) to birth-death process birth rates.";
    }

    public Class getReturnType() {
        return TwoPathogenModelBirthRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SEASONAL_MODEL, true),
            new ElementRule(SEASONAL_PERIOD,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(SEASONAL_AMP,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(SEASONAL_PEAK,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(TRANSMISSION_RATE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_SS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(NUM_SI,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(NUM_SC,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(NUM_SR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(NUM_IS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    },true),
            new ElementRule(NUM_CS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    },true),
            new ElementRule(NUM_RS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    },true),
            new ElementRule(ORIGIN_TIME_NUM_SS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_GRID_POINTS,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(CUT_OFF,
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
            new ElementRule(PATHOGEN_NUMBER,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOST_RECENT_SAMPLING_DATE_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(MOST_RECENT_SAMPLING_DATE_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
    };

}
