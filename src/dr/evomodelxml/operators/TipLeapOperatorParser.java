package dr.evomodelxml.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class TipLeapOperatorParser extends AbstractXMLObjectParser {

    public static final String TIP_LEAP = "tipLeap";

    public static final String SIZE = "size";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String DISTANCE_KERNEL = "distanceKernel";
    public static final String TAXA = "taxa";

    public String getParserName() {
        return TIP_LEAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Taxa taxa = (Taxa) xo.getChild(Taxa.class);

        List<NodeRef> tips = new ArrayList<NodeRef>();

        for (Taxon taxon : taxa) {
            boolean found = false;
            for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
                NodeRef tip = treeModel.getExternalNode(i);
                if (treeModel.getNodeTaxon(tip).equals(taxon)) {
                    tips.add(tip);
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new XMLParseException("Error constructing " + TIP_LEAP +": " + taxon.getId() + ", not found in tree with id " + treeModel.getId());
            }
        }

        // size attribute is mandatory
        final double size = xo.getAttribute(SIZE, Double.NaN);
        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);
        final String distanceKernel = xo.getAttribute(DISTANCE_KERNEL, "gaussian");

        if (size <= 0.0) {
            throw new XMLParseException("The TipLeap size attribute must be positive and non-zero.");
        }

        if (targetAcceptance <= 0.0 || targetAcceptance >= 1.0) {
            throw new XMLParseException("Target acceptance probability has to lie in (0, 1)");
        }

        return new SubtreeLeapOperator(treeModel, taxa, weight, size, targetAcceptance, distanceKernel, mode);
    }

    public String getParserDescription() {
        return "An operator that moves a tip a certain patristic distance.";
    }

    public Class getReturnType() {
        return SubtreeLeapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(SIZE, false),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class),
            new ElementRule(Taxa.class)
    };

}
