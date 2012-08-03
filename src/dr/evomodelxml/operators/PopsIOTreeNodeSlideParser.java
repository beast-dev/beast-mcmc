package dr.evomodelxml.operators;

import dr.evomodel.operators.PopsIOTreeNodeSlide;
import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * User: Graham Jones
 * Date: 11/05/12
 */
public class PopsIOTreeNodeSlideParser  extends AbstractXMLObjectParser {
    public static final String PIOTREE_NODESLIDE = "pioTreeNodeSlide";



    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        PopsIOSpeciesBindings piosb = (PopsIOSpeciesBindings) xo.getChild(PopsIOSpeciesBindings.class);
        PopsIOSpeciesTreeModel piostm = (PopsIOSpeciesTreeModel) xo.getChild(PopsIOSpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new PopsIOTreeNodeSlide(piostm, piosb, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(PopsIOSpeciesBindings.class),
                new ElementRule(PopsIOSpeciesTreeModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator for species tree: changes tree without breaking embedding of gene trees.";
    }

    @Override
    public Class getReturnType() {
        return PopsIOTreeNodeSlide.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getParserName() {
        return PIOTREE_NODESLIDE;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
