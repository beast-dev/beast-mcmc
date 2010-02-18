package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.coalescent.TwoEpochDemographicModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class TwoEpochDemographicModelParser extends AbstractXMLObjectParser {

    public static final String EPOCH_1 = "modernEpoch";
    public static final String EPOCH_2 = "ancientEpoch";
    public static final String TRANSITION_TIME = "transitionTime";

    public static final String TWO_EPOCH_MODEL = "twoEpoch";

    public String getParserName() {
        return TWO_EPOCH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(EPOCH_1);
        DemographicModel demo1 = (DemographicModel) cxo.getChild(DemographicModel.class);

        cxo = xo.getChild(EPOCH_2);
        DemographicModel demo2 = (DemographicModel) cxo.getChild(DemographicModel.class);

        cxo = xo.getChild(TRANSITION_TIME);
        Parameter timeParameter = (Parameter) cxo.getChild(Parameter.class);

        return new TwoEpochDemographicModel(demo1, demo2, timeParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of two epochs.";
    }

    public Class getReturnType() {
        return TwoEpochDemographicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(EPOCH_1,
                    new XMLSyntaxRule[]{new ElementRule(DemographicModel.class)},
                    "The demographic model for the recent epoch."),
            new ElementRule(EPOCH_2,
                    new XMLSyntaxRule[]{new ElementRule(DemographicModel.class)},
                    "The demographic model for the ancient epoch."),
            new ElementRule(TRANSITION_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The time that splits the two epochs."),
            XMLUnits.SYNTAX_RULES[0]
    };
}
