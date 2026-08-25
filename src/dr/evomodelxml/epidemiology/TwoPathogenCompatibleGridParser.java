package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.TwoPathogenCompatibleGrid;
import dr.inference.model.Parameter;
import dr.xml.*;


public class TwoPathogenCompatibleGridParser extends AbstractXMLObjectParser {

    public static final String TWO_PATHOGEN_COMPATIBLE_GRID = "twoPathogenCompatibleGrid";
    public static final String CUT_OFF = "cutOff";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String MOST_RECENT_SAMPLING_DATE_ONE = "mostRecentSamplingDateOne";
    public static final String MOST_RECENT_SAMPLING_DATE_TWO = "mostRecentSamplingDateTwo";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter cutOff = (Parameter) xo.getChild(CUT_OFF).getChild(Parameter.class);
        Parameter numGridPoints = (Parameter) xo.getChild(NUM_GRID_POINTS).getChild(Parameter.class);
        Parameter mostRecentSamplingDateOne = (Parameter) xo.getChild(MOST_RECENT_SAMPLING_DATE_ONE).getChild(Parameter.class);
        Parameter mostRecentSamplingDateTwo = (Parameter) xo.getChild(MOST_RECENT_SAMPLING_DATE_TWO).getChild(Parameter.class);

        return new TwoPathogenCompatibleGrid(cutOff, numGridPoints,
                mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CUT_OFF,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(NUM_GRID_POINTS,
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

    public String getParserDescription() {
        return "Grid points (starting times of intervals) in episodic birth-death sampling model for pathogen with narrower time scale in two-pathogen model";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return TWO_PATHOGEN_COMPATIBLE_GRID;
    }
}
