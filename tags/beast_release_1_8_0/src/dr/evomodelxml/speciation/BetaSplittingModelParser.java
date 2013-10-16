package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.BetaSplittingModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a SpeciationModel. Recognises YuleModel.
 */
public class BetaSplittingModelParser extends AbstractXMLObjectParser {

    public static final String BETA_SPLITTING_MODEL = "betaSplittingModel";
    public static final String PHI = "phi";
    public static final String TREE = "branchingTree";


    public String getParserName() {
        return BETA_SPLITTING_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(PHI);
        Parameter phiParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TREE);
        Tree tree = (Tree) cxo.getChild(Tree.class);

        return new BetaSplittingModel(phiParameter, tree);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "The beta-splitting family of tree branching models (Aldous, 1996;2001).";
    }

    public Class getReturnType() {
        return BetaSplittingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PHI,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, "A parameter that ranges from -infinity (comb-tree) to +infinity (balanced tree)"),
            new ElementRule(TREE,
                    new XMLSyntaxRule[]{new ElementRule(Tree.class)})
    };
}
