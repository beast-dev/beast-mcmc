package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.inference.model.IndexedParameter;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/*
 * @author Filippo Monti
 */
public final class RewardsAwareMixtureBranchRatesParser extends AbstractXMLObjectParser {

    public static final String NAME = "rewardsAwareMixtureBranchRates";
    public static final String ATOMIC = "atomic";
    public static final String CTS = "ctsParameter";
    public static final String INDICATOR = "indicator";
    public static final String INCLUDE_ROOT = "includeRoot";

    @Override
    public String getParserName() { return NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        IndexedParameter atomic = (IndexedParameter) xo.getChild(ATOMIC).getChild(IndexedParameter.class);
        Parameter cts = (Parameter) xo.getChild(CTS).getChild(Parameter.class);
        Parameter indicator = (Parameter) xo.getChild(INDICATOR).getChild(Parameter.class);

        TreeParameterModel.Type includeRoot = xo.getAttribute(INCLUDE_ROOT, false)
                ? TreeParameterModel.Type.WITH_ROOT
                : TreeParameterModel.Type.WITHOUT_ROOT;

        final int nBranches = tree.getNodeCount() - 1;
        if (includeRoot == TreeParameterModel.Type.WITHOUT_ROOT) {
            if (cts.getDimension() != nBranches) throw new XMLParseException("cts dim must be nodeCount-1");
            if (indicator.getDimension() != nBranches) throw new XMLParseException("indicator dim must be nodeCount-1");
            if (atomic.getDimension() != nBranches) throw new XMLParseException("atomic dim must be nodeCount-1");
        }
        final ArbitraryBranchRates.BranchRateTransform transform = ArbitraryBranchRatesParser.parseTransform(xo);

        return new RewardsAwareMixtureBranchRates(tree, cts, indicator, atomic, transform, false, includeRoot);
    }

    @Override
    public String getParserDescription() {
        return "Mixture branch-rate model: per-branch continuous rate or atomic rate selected via an IndexedParameter, controlled by an indicator.";
    }

    @Override
    public Class getReturnType() { return RewardsAwareMixtureBranchRates.class; }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(TreeModel.class),
                new ElementRule(ATOMIC, new XMLSyntaxRule[]{ new ElementRule(IndexedParameter.class) }),
                new ElementRule(CTS, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
                new ElementRule(INDICATOR, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
                AttributeRule.newBooleanRule(INCLUDE_ROOT, true),
                new ElementRule(ArbitraryBranchRatesParser.SCALE, Parameter.class, "optional scale parameter", true),
                new ElementRule(ArbitraryBranchRatesParser.LOCATION, Parameter.class, "optional location parameter", true),
                new ElementRule(ArbitraryBranchRatesParser.RANDOM_INDICATOR, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }, true),
                AttributeRule.newBooleanRule(ArbitraryBranchRatesParser.SHRINKAGE, true),
        };
    }

}