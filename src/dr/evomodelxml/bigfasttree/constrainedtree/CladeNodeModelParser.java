package dr.evomodelxml.bigfasttree.constrainedtree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.constrainedtree.CladeNodeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.xml.*;

public class CladeNodeModelParser extends AbstractXMLObjectParser {

    public static final String CLADE_MODEL = "cladeNodeModel";
    public static final String CLADE_TREE = "cladeTree";

    public String getParserName() {
        return CLADE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree cladeTree = (Tree) xo.getElementFirstChild(CLADE_TREE);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        try {
            return new CladeNodeModel(CLADE_MODEL, cladeTree, treeModel);
        } catch (TreeUtils.MissingTaxonException e) {
            e.printStackTrace();
            throw new XMLParseException("tree and clade much match");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a clade structure that provides access to nodes in a semifixed topology";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(CLADE_TREE, new XMLSyntaxRule[] {
                    new ElementRule(Tree.class)
            }),
            new ElementRule(TreeModel.class),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}