package dr.evomodelxml.speciation;

import dr.evomodel.speciation.AlloppNumHybsStatistic;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;

import dr.xml.*;

/**
 * @author Graham Jones
 * Date: 08/10/2012
 */
public class AlloppNumHybsStatisticParser  extends AbstractXMLObjectParser {
    public static final String NUMHYBS_STATISTIC = "alloppNumHybsStatistic";
    public static final String APSPNETWORK = "apspNetwork";

    public String getParserName() {
        return NUMHYBS_STATISTIC;
    }


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final XMLObject asnmxo = xo.getChild(APSPNETWORK);
        AlloppSpeciesNetworkModel aspnet = (AlloppSpeciesNetworkModel) asnmxo.getChild(AlloppSpeciesNetworkModel.class);
        return new AlloppNumHybsStatistic(aspnet);
    }



    private  XMLSyntaxRule[] asnmRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(APSPNETWORK, asnmRules()),

        };
    }



    @Override
    public String getParserDescription() {
        return "Statistic for number of hybridizations in allopolyploid network";
    }

    @Override
    public Class getReturnType() {
        return AlloppNumHybsStatistic.class;
    }
}
