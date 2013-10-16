
package dr.evomodelxml.operators;

import dr.evomodel.operators.MulTreeNodeSlide;
import dr.evomodel.speciation.MulSpeciesBindings;
import dr.evomodel.speciation.MulSpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */
public class MulTreeNodeSlideParser extends AbstractXMLObjectParser {
    public static final String MULTREE_NODE_REHEIGHT = "mulTreeNodeReHeight";

    public String getParserName() {
        return MULTREE_NODE_REHEIGHT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	MulSpeciesBindings species = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);
    	MulSpeciesTreeModel tree = (MulSpeciesTreeModel) xo.getChild(MulSpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
//            final double range = xo.getAttribute("range", 1.0);
//            if( range <= 0 || range > 1.0 ) {
//                throw new XMLParseException("range out of range");
//            }
        //final boolean oo = xo.getAttribute("outgroup", false);
        return new MulTreeNodeSlide(tree, species /*, range*//*, oo*/, weight);
    }

    public String getParserDescription() {
        return "Specialized Species tree operator, transform tree without breaking embedding of gene trees.";
    }

    public Class getReturnType() {
        return MulTreeNodeSlide.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
               // AttributeRule.newDoubleRule("range", true),
              //  AttributeRule.newBooleanRule("outgroup", true),

                new ElementRule(MulSpeciesBindings.class),
                new ElementRule(MulSpeciesTreeModel.class)
        };
    }

}
