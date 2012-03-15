package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFIntervalHeightsStatistic;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class GMRFIntervalHeightsStatisticParser extends AbstractXMLObjectParser {

    public static final String GMRF_HEIGHTS_STATISTIC = "gmrfHeightsStatistic";

    public String getParserName() {
        return GMRF_HEIGHTS_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        GMRFSkyrideLikelihood skyrideLikelihood = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);

        return new GMRFIntervalHeightsStatistic(name, skyrideLikelihood);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the heights of each internal node in increasing order (or groups them by a group size parameter)";
    }

    public Class getReturnType() {
        return GMRFIntervalHeightsStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(Statistic.NAME, true),
            new ElementRule(GMRFSkyrideLikelihood.class),
    };

}