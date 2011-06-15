package dr.app.beagle.evomodel.treelikelihood;

import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeModel;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Takes a TreeTraitProvider that returns multiple-site traits and provides a new Trait for logging by site
 *
 * @author Marc A. Suchard
 */
public class SplitBySiteTraitLogger extends TreeTraitProvider.Helper implements Citable {

    public static final String TRAIT_LOGGER = "splitTraitBySite";
    public static final String TRAIT_NAME = "traitName";
    public static final String SCALE = "scaleByBranchLength";

    public SplitBySiteTraitLogger(AncestralStateBeagleTreeLikelihood treeLikelihood, String traitName, boolean scale) throws XMLParseException {
        TreeTrait trait = treeLikelihood.getTreeTrait(traitName);
        if (trait == null) {
            throw new XMLParseException("TreeTraitProvider does not provide trait named '" + traitName + ".");
        }

        TreeModel tree = treeLikelihood.getTreeModel();
        int length;
        Object obj = trait.getTrait(tree, tree.getNode(0));

        if (obj instanceof double[]) {
            length = ((double[]) obj).length;
        } else {
            throw new XMLParseException("Unknown trait type to split");
        }

        TreeTrait[] partitionedTraits = new TreeTrait[length];
        for (int i = 0; i < length; i++) {
            if (scale) {
                partitionedTraits[i] = new TreeTrait.PickEntryDAndScale(trait, i);
            } else {
                partitionedTraits[i] = new TreeTrait.PickEntryD(trait, i);
            }
        }
        addTraits(partitionedTraits);

        Logger.getLogger("dr.app.beagle").info("\tConstructing a split logger with " + length + " partitions;  please cite:\n"
                + Citable.Utils.getCitationString(this));

    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRAIT_LOGGER;
        }

        /**
         * @return an object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String traitName = xo.getStringAttribute(TRAIT_NAME);
            AncestralStateBeagleTreeLikelihood tree = (AncestralStateBeagleTreeLikelihood) xo.getChild(AncestralStateBeagleTreeLikelihood.class);
            boolean scale = xo.getAttribute(SCALE, false);

            return new SplitBySiteTraitLogger(tree, traitName, scale);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(AncestralStateBeagleTreeLikelihood.class, "The tree which is to be logged"),
                AttributeRule.newStringRule(TRAIT_NAME),
                AttributeRule.newBooleanRule(SCALE, true),
        };

        public String getParserDescription() {
            return null;
        }

        public String getExample() {
            return null;
        }

        public Class getReturnType() {
            return TreeTraitProvider.class;
        }
    };

    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                CommonCitations.SUCHARD_2012
        );
        return citations;
    }
}
