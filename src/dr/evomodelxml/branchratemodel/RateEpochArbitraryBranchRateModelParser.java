package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RateEpochArbitraryBranchRateModel;
import dr.evomodel.branchratemodel.RateEpochBranchRateModel;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Karthik Gangavarapu
 * @author Marc Suchard
 */
public class RateEpochArbitraryBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String RATE_EPOCH_ARBITRARY_BRANCH_RATES = "rateEpochArbitraryBranchRates";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ArbitraryBranchRates arbitraryBranchRates = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);
        RateEpochBranchRateModel epochBranchRateModel = (RateEpochBranchRateModel) xo.getChild(RateEpochBranchRateModel.class);
        return new RateEpochArbitraryBranchRateModel(epochBranchRateModel, arbitraryBranchRates);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(ArbitraryBranchRates.class),
            new ElementRule(RateEpochBranchRateModel.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return RateEpochArbitraryBranchRateModel.class;
    }

    @Override
    public String getParserName() {
        return RATE_EPOCH_ARBITRARY_BRANCH_RATES;
    }
}
