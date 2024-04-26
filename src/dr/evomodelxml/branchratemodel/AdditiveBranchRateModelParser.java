package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.AdditiveBranchRateModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Karthik Gangavarapu
 * @author Marc Suchard
 */
public class AdditiveBranchRateModelParser extends AbstractXMLObjectParser {

    public static String ADDITIVE_BRANCH_RATES = "additiveBranchRates";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<AbstractBranchRateModel> branchRateModels = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            branchRateModels.add((AbstractBranchRateModel) xo.getChild(i));
        }
        return new AdditiveBranchRateModel(branchRateModels);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AbstractBranchRateModel.class, 2, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return AdditiveBranchRateModel.class;
    }

    @Override
    public String getParserName() {
        return ADDITIVE_BRANCH_RATES;
    }
}
