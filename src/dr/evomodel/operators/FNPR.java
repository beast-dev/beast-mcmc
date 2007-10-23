/**
 * 
 */
package dr.evomodel.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * This is an implementation of the Subtree Prune and Regraft (SPR) operator for
 * trees. It assumes explicitely bifurcating rooted trees.
 * 
 * @author Sebastian Hoehna
 * @version 1.0
 * 
 */
public class FNPR extends SimpleMCMCOperator {

   private TreeModel _tree = null;

   public static final String FNPR = "FixedNodeheightSubtreePruneRegraft";

   int mode = CoercableMCMCOperator.DEFAULT;

   /**
    * 
    */
   public FNPR(TreeModel tree, int weight) {
      _tree = tree;
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
      NodeRef iFather, iGrandfather, iBrother;
      double heightFather;
      int tipCount = _tree.getExternalNodeCount();

      final int nNodes = _tree.getNodeCount();
      final NodeRef root = _tree.getRoot();

      NodeRef i;

      int MAX_TRIES = 1000;

      for (int tries = 0; tries < MAX_TRIES; ++tries) {
         // get a random node where you are not the root
         do {
            i = _tree.getNode(MathUtils.nextInt(nNodes));
            // i and his father should not be the root because then you
            // can
            // not
            // find a new node which is older
         } while (root == i || _tree.getParent(i) == root);

         // int childIndex = (MathUtils.nextDouble() >= 0.5 ? 1 : 0);
         // int otherChildIndex = 1 - childIndex;
         // NodeRef iOtherChild = _tree.getChild(i, otherChildIndex);

         iFather = _tree.getParent(i);
         iGrandfather = _tree.getParent(iFather);
         iBrother = getOtherChild(_tree, iFather, i);
         heightFather = _tree.getNodeHeight(iFather);

         // NodeRef newChild = getRandomNode(possibleChilds, iFather);
         NodeRef newChild = _tree.getNode(MathUtils.nextInt(nNodes));

         if (_tree.getNodeHeight(newChild) < heightFather
               && root != newChild
               && _tree.getNodeHeight(_tree.getParent(newChild)) > heightFather
               && newChild != iFather
               && _tree.getParent(newChild) != iFather) {
            NodeRef newGrandfather = _tree.getParent(newChild);

            _tree.beginTreeEdit();

            // prune
            _tree.removeChild(iFather, iBrother);
            _tree.removeChild(iGrandfather, iFather);
            _tree.addChild(iGrandfather, iBrother);

            // reattach
            _tree.removeChild(newGrandfather, newChild);
            _tree.addChild(iFather, newChild);
            _tree.addChild(newGrandfather, iFather);

            // ****************************************************

            try {
               _tree.endTreeEdit();
            } catch (MutableTree.InvalidTreeException ite) {
               throw new OperatorFailedException(ite.toString());
            }

            // TODO
            _tree.pushTreeChangedEvent(i);
            // _tree.pushTreeChangedEvent(newGrandfather);
            // _tree.pushTreeChangedEvent(newGrandfather);

            if (_tree.getExternalNodeCount() != tipCount) {
               throw new RuntimeException(
                     "Lost some tips in SPR operation");
            }

            // return hastingsRatio;
            return 0.0;
         }
      }

      throw new OperatorFailedException(
            "Couldn't find valid SPR move on this tree!");
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
      // TODO Auto-generated method stub
      return null;
   }

   public static XMLObjectParser FNPR_PARSER = new AbstractXMLObjectParser() {

      public String getParserName() {
         return FNPR;
      }

      public Object parseXMLObject(XMLObject xo) throws XMLParseException {

         TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
         int weight = xo.getIntegerAttribute("weight");

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

      private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule("weight"),
            new ElementRule(TreeModel.class) };

   };

}
