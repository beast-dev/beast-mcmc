package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.CoalescentIntervalStatistic;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;


public class CoalescentIntervalStatisticParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_INTERVAL_STATISTIC = "coalescentIntervalStatistic";
    public static final String GMRFLIKELIHOOD = "gmrfSkyrideLikelihood";

    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return CoalescentIntervalStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
        		new XORRule(
                new ElementRule(GMRFSkyrideLikelihood.class),
        		new ElementRule(CoalescentLikelihood.class))
        };
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	if (xo.getChild(GMRFSkyrideLikelihood.class) != null) {
    		return new CoalescentIntervalStatistic((GMRFSkyrideLikelihood)xo.getChild(GMRFSkyrideLikelihood.class));
    	} else {
    		return new CoalescentIntervalStatistic((CoalescentLikelihood)xo.getChild(CoalescentLikelihood.class));
    	}
    }

    public String getParserName() {
        return COALESCENT_INTERVAL_STATISTIC;
    }

}
