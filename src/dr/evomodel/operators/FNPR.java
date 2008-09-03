/**
 *
 */
package dr.evomodel.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * This is an implementation of the Subtree Prune and Regraft (SPR) operator for
 * trees. It assumes explicitely bifurcating rooted trees.
 *
 * @author Sebastian Hoehna
 * @version 1.0
 */
public class FNPR extends AbstractTreeOperator {

    private TreeModel tree = null;

    public static final String FNPR = "FixedNodeheightSubtreePruneRegraft";

    /**
     *
     */
    public FNPR(TreeModel tree, double weight) {
        this.tree = tree;
        setWeight(weight);
        // distances = new int[tree.getNodeCount()];
    }

    /*
    * (non-Javadoc)
    *
    * @see dr.inference.operators.SimpleMCMCOperator#doOperation()
    */
    @Override
    public double doOperation() throws OperatorFailedException {
        NodeRef iGrandfather, iBrother;
        double heightFather;
        final int tipCount = tree.getExternalNodeCount();

        final int nNodes = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i;

        int MAX_TRIES = 1000;

        for (int tries = 0; tries < MAX_TRIES; ++tries) {
           // get a random node whose father is not the root - otherwise
           // the operation is not possible
           do {
              i = tree.getNode(MathUtils.nextInt(nNodes));
           } while (root == i || tree.getParent(i) == root);

           // int childIndex = (MathUtils.nextDouble() >= 0.5 ? 1 : 0);
           // int otherChildIndex = 1 - childIndex;
           // NodeRef iOtherChild = tree.getChild(i, otherChildIndex);

           NodeRef iFather = tree.getParent(i);
           iGrandfather = tree.getParent(iFather);
           iBrother = getOtherChild(tree, iFather, i);
           heightFather = tree.getNodeHeight(iFather);

           // NodeRef newChild = getRandomNode(possibleChilds, iFather);
           NodeRef newChild = tree.getNode(MathUtils.nextInt(nNodes));

           if (tree.getNodeHeight(newChild) < heightFather
                 && root != newChild
                 && tree.getNodeHeight(tree.getParent(newChild)) > heightFather
                 && newChild != iFather
                 && tree.getParent(newChild) != iFather) {
              NodeRef newGrandfather = tree.getParent(newChild);

              tree.beginTreeEdit();

              // prune
              tree.removeChild(iFather, iBrother);
              tree.removeChild(iGrandfather, iFather);
              tree.addChild(iGrandfather, iBrother);

              // reattach
              tree.removeChild(newGrandfather, newChild);
              tree.addChild(iFather, newChild);
              tree.addChild(newGrandfather, iFather);

              // ****************************************************

              try {
                 tree.endTreeEdit();
              } catch (MutableTree.InvalidTreeException ite) {
                 throw new OperatorFailedException(ite.toString());
              }

              tree.pushTreeChangedEvent(i);

              assert tree.getExternalNodeCount() == tipCount;
              
              return 0.0;
           }
        }

        throw new OperatorFailedException("Couldn't find valid SPR move on this tree!");
     }

    /*
    * (non-Javadoc)
    *
    * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
    */
    @Override
    public String getOperatorName() {
        return FNPR;
    }

    public double getTargetAcceptanceProbability() {

        return 0.0234;
    }

    public double getMaximumAcceptanceLevel() {

        return 0.04;
    }

    public double getMaximumGoodAcceptanceLevel() {

        return 0.03;
    }

    public double getMinimumAcceptanceLevel() {

        return 0.005;
    }

    public double getMinimumGoodAcceptanceLevel() {

        return 0.01;
    }

    /*
    * (non-Javadoc)
    *
    * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
    */
    public String getPerformanceSuggestion() {
        return "";
    }

    public static XMLObjectParser FNPR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FNPR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            double weight = xo.getDoubleAttribute("weight");

            return new FNPR(treeModel, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription() {
            return "This element represents a FNPR operator. "
                    + "This operator swaps a random subtree with its uncle.";
        }

        public Class getReturnType() {
            return FNPR.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(TreeModel.class)};

    };

}
