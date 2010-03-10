package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.CoalescentIntervalStatistic;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.xml.*;

/**
 */
public class CoalescentIntervalStatisticParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_INTERVAL_STATISTIC = "coalescentIntervalStatistic";

    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return CoalescentIntervalStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(GMRFSkyrideLikelihood.class),
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new CoalescentIntervalStatistic(
                (GMRFSkyrideLikelihood)xo.getChild(GMRFSkyrideLikelihood.class));
    }

    public String getParserName() {
        return COALESCENT_INTERVAL_STATISTIC;
    }

}
