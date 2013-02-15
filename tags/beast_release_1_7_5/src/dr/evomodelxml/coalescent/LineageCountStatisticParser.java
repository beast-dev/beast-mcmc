package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.LineageCountStatistic;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.xml.*;

/**
 */
public class LineageCountStatisticParser extends AbstractXMLObjectParser {

    public static final String LINEAGE_COUNT_STATISTIC = "lineageCountStatistic";


    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return LineageCountStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(OldAbstractCoalescentLikelihood.class),
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new LineageCountStatistic(
                (OldAbstractCoalescentLikelihood) xo.getChild(OldAbstractCoalescentLikelihood.class));
    }

    public String getParserName() {
        return LINEAGE_COUNT_STATISTIC;
    }

}
