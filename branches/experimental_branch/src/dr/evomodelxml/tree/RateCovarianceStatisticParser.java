package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.RateCovarianceStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class RateCovarianceStatisticParser extends AbstractXMLObjectParser {

    public static final String RATE_COVARIANCE_STATISTIC = "rateCovarianceStatistic";

    public String getParserName() {
        return RATE_COVARIANCE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        return new RateCovarianceStatistic(name, tree, branchRateModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that has as its value the covariance of parent and child branch rates";
    }

    public Class getReturnType() {
        return RateCovarianceStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
    };

}
