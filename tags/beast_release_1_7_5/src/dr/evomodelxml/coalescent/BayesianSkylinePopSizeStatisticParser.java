package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.coalescent.BayesianSkylinePopSizeStatistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 */
public class BayesianSkylinePopSizeStatisticParser extends AbstractXMLObjectParser {
    public static final String TIME = "time";
    public static final String BAYESIAN_SKYLINE_POP_SIZE_STATISTIC = "generalizedSkylinePopSizeStatistic";


    public String getParserDescription() {
        return "The pop sizes at the given times";
    }

    public Class getReturnType() {
        return BayesianSkylinePopSizeStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return null;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double time = xo.getDoubleAttribute(TIME);

        BayesianSkylineLikelihood bsl =
                (BayesianSkylineLikelihood) xo.getChild(BayesianSkylineLikelihood.class);

        return new BayesianSkylinePopSizeStatistic(time, bsl);
    }

    public String getParserName() {
        return BAYESIAN_SKYLINE_POP_SIZE_STATISTIC;
    }

}
