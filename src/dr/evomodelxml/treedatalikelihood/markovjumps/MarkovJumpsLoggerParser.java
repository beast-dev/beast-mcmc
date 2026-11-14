package dr.evomodelxml.treedatalikelihood.markovjumps;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.ProcessAlongTree;
import dr.evomodel.treedatalikelihood.markovjumps.MarkovJumpsLogger;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.markovjumps.CompleteHistoryAddOn.UNCONDITIONAL_COUNT_SUFFIX;
import static dr.evomodel.treedatalikelihood.markovjumps.MarkovJumpRewardAddOn.MARGINAL_RATE_SUFFIX;
import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedDiscreteTraitDelegate.SITE_RATE_SUFFIX;
import static dr.evomodelxml.tree.TreeLoggerParser.*;

public class MarkovJumpsLoggerParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "markovJumpsLogger";
    private static final String REPORT_UNCONDITIONAL = "reportUnconditional";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ProcessAlongTree process = (ProcessAlongTree) xo.getChild(ProcessAlongTree.class);
        String name = xo.getStringAttribute(NAME);
        final TreeTrait trait = process.getTreeTrait(name);

        if (trait == null) {

            String childName = "ProcessAlongTree";
            if (process instanceof Likelihood) {
                childName = ((Likelihood) process).prettyName();
            } else if (process instanceof Model) {
                childName = ((Model) process).getModelName();
            }

            throw new XMLParseException("Trait named '" + name + "' not found for '" + childName + "'");
        }

        if (!trait.getLoggable()) {
            throw new XMLParseException("Trait named '" + name + "' is not loggable");
        }

        String tag = xo.getAttribute(TAG, name);
        boolean reportUnconditional = xo.getAttribute(REPORT_UNCONDITIONAL, false);

        TreeTrait marginalRateTrait = process.getTreeTrait(name + MARGINAL_RATE_SUFFIX);
        TreeTrait unconditionalExpectedJumps = process.getTreeTrait(name + UNCONDITIONAL_COUNT_SUFFIX);
        TreeTrait siteSpecificRate = process.getTreeTrait(name + SITE_RATE_SUFFIX);

        return new MarkovJumpsLogger(tag, process, trait,
                marginalRateTrait, unconditionalExpectedJumps, siteSpecificRate,
                reportUnconditional);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(NAME, false, "The name of the trait"),
            AttributeRule.newStringRule(TAG, true, "The label of the trait to be used in the log"),
            AttributeRule.newBooleanRule(REPORT_UNCONDITIONAL, true),
            new ElementRule(ProcessAlongTree.class, "The trait provider"),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return MarkovJumpsLogger.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
