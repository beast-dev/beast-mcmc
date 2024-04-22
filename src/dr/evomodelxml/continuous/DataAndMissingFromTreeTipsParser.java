package dr.evomodelxml.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Gabe Hassler
 * @author Marc Suchard
 */

public class DataAndMissingFromTreeTipsParser extends AbstractXMLObjectParser{

        public final static String DATA_FROM_TREE_TIPS = "dataAndMissingFromTreeTips";
        public final static String DATA = "data";
        public static final String CONTINUOUS = "continuous";


        public String getParserName() {
            return DATA_FROM_TREE_TIPS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, treeModel, true);

            return returnValue;
        }

        public static final XMLSyntaxRule[] rules = {
                new ElementRule(MutableTreeModel.class),
                AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                })
        };

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Returns the data and missing indices from the tips of a tree.";
        }

        @Override
        public Class getReturnType() {
            return TreeTraitParserUtilities.TraitsAndMissingIndices.class;
        }


}