package dr.evomodelxml.operators;

import dr.evomodel.operators.TreeNodeSlide;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TreeNodeSlideParser extends AbstractXMLObjectParser {
    public static final String TREE_NODE_REHEIGHT = "nodeReHeight";

    public String getParserName() {
        return TREE_NODE_REHEIGHT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SpeciesBindings species = (SpeciesBindings) xo.getChild(SpeciesBindings.class);
        SpeciesTreeModel tree = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
//            final double range = xo.getAttribute("range", 1.0);
//            if( range <= 0 || range > 1.0 ) {
//                throw new XMLParseException("range out of range");
//            }
        //final boolean oo = xo.getAttribute("outgroup", false);
        return new TreeNodeSlide(tree, species /*, range*//*, oo*/, weight);
    }

    public String getParserDescription() {
        return "Specialized Species tree operator, transform tree without breaking embedding of gene trees.";
    }

    public Class getReturnType() {
        return TreeNodeSlide.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
               // AttributeRule.newDoubleRule("range", true),
              //  AttributeRule.newBooleanRule("outgroup", true),

                new ElementRule(SpeciesBindings.class),
                new ElementRule(SpeciesTreeModel.class)
        };
    }

}
