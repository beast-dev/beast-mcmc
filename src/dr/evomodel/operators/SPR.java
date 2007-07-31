/**
 * 
 */
package dr.evomodel.operators;

import java.util.HashMap;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
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
public class SPR extends SimpleMCMCOperator {

   private boolean _withoutNNI = false;

   private boolean _justNNI = false;

   private TreeModel _tree = null;

   public static final String SPR = "SubtreePruneRegraft";

   /**
    * 
    */
   public SPR(TreeModel tree, int weight) {
      _tree = tree;
      setWeight(weight);
   }

   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.operators.SimpleMCMCOperator#doOperation()
    */
   @Override
   public double doOperation() throws OperatorFailedException {
      int tipCount = _tree.getExternalNodeCount();

      final int nNodes = _tree.getNodeCount();
      final NodeRef root = _tree.getRoot();

      NodeRef i;

      // get a random node where you are not the root
      do {
         i = _tree.getNode(MathUtils.nextInt(nNodes));
      } while (root == i || _tree.getParent(i) == root);

      // prune the subtree (starting from my father)
      // get parent node
      final NodeRef iFather = _tree.getParent(i);
      // get your brother
      NodeRef iBrother = getOtherChild(_tree, iFather, i);
      double heightBrother = _tree.getNodeHeight(iBrother);
      NodeRef iUncle = getOtherChild(_tree, _tree.getParent(iFather), iFather);

      _tree.beginTreeEdit();

      final NodeRef iGrandfather = _tree.getParent(iFather);
      double heightGrandfather = _tree.getNodeHeight(iGrandfather);
      // set my brother as the new son of my gradnfather
      _tree.removeChild(iGrandfather, iFather);
      _tree.removeChild(iFather, iBrother);
      _tree.addChild(iGrandfather, iBrother);

      // find a random node which is NOT on the subtree we pruned before
      // above this node we will stick the subtree
      NodeRef j;
      do {
         j = getRandomNodeOnOtherSubtree(iFather);
      } while (j == iBrother);
      NodeRef jFather = _tree.getParent(j);

      // special variant if we just want neighbors which aren't NNI neighbors
      // if (_withoutNNI) {
      // if (jFather != root) {
      // NodeRef jUncle = getOtherChild(_tree, _tree.getParent(jFather),
      // jFather);
      // while (j == iUncle || i == jUncle) {
      // j = getRandomNodeOnOtherSubtree(iFather);
      // }
      // jFather = _tree.getParent(j);
      // }
      // }

      if (_justNNI) {
         while (j != iUncle) {
            j = getRandomNodeOnOtherSubtree(iFather);
         }
         jFather = _tree.getParent(j);
      }

      // change the height of the father to be between the height of the
      // new child and the father of him
      // because father becomes the new father of the new child
      // and necessary for hasting ratio
      double ran = Math.random();
      // check if the random value equals 0 or 1 because then the new
      // height would be the same as
      // the uncle's or grandparent's height
      if (ran == 0.0) {
         ran = 0.1;
      } else if (ran == 1.0) {
         ran = 0.9;
      }
      // now calculate the new height for father between the height of the
      // uncle and the grandparent
      double heightNewBrother = _tree.getNodeHeight(j);
      double heightNewGrandfather = _tree.getNodeHeight(jFather);
      double factor = ((MathUtils.nextDouble() * (heightNewGrandfather - heightNewBrother)) + heightNewBrother)
            / _tree.getNodeHeight(iFather);
      // set the new height for the father and his childs
      scaleSubtree(_tree, iFather, factor);

      // set the hastings ratio
       double hastingsRatio = Math.log(Math.min(1,
       (heightNewGrandfather - heightNewBrother)
       / (heightGrandfather - heightBrother)));

      _tree.removeChild(jFather, j);
      _tree.addChild(jFather, iFather);
      _tree.addChild(iFather, j);

      try {
         _tree.endTreeEdit();
      } catch (MutableTree.InvalidTreeException ite) {
         throw new OperatorFailedException(ite.toString());
      }

      // TODO
      _tree.pushTreeChangedEvent(iFather);
      _tree.pushTreeChangedEvent(jFather);

      if (_tree.getExternalNodeCount() != tipCount) {
         throw new RuntimeException("Lost some tips in SPR operation");
      }

      return hastingsRatio;
   }

   private NodeRef getRandomNodeOnOtherSubtree(NodeRef subtreeRoot) {
      HashMap<NodeRef, Boolean> subtreeNodes = new HashMap<NodeRef, Boolean>();

      fillHashMap(subtreeNodes, subtreeRoot);

      int nNodes = _tree.getNodeCount();
      NodeRef root = _tree.getRoot();
      NodeRef i;
      do {
         i = _tree.getNode(MathUtils.nextInt(nNodes));
      } while (root == i || subtreeNodes.containsKey(i));

      return i;
   }

   private void fillHashMap(HashMap<NodeRef, Boolean> nodes, NodeRef root) {
      nodes.put(root, true);
      int childCount = _tree.getChildCount(root);
      NodeRef child;
      for (int i = 0; i < childCount; i++) {
         child = _tree.getChild(root, i);
         if (child == null) {
            child = _tree.getChild(root, 1);
         }
         fillHashMap(nodes, child);
      }
   }

   
   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
    */
   @Override
   public String getOperatorName() {
      return SPR;
   }

   public double getMinimumAcceptanceLevel() {
      return 0.01;
   }

   public double getMinimumGoodAcceptanceLevel() {
      return 0.005;
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

   public static XMLObjectParser SPR_PARSER = new AbstractXMLObjectParser() {

      public String getParserName() {
         return SPR;
      }

      public Object parseXMLObject(XMLObject xo) throws XMLParseException {

         TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
         int weight = xo.getIntegerAttribute("weight");

         return new SPR(treeModel, weight);
      }

      // ************************************************************************
      // AbstractXMLObjectParser implementation
      // ************************************************************************

      public String getParserDescription() {
         return "This element represents a SPR operator. "
               + "This operator swaps a random subtree with its uncle.";
      }

      public Class getReturnType() {
         return SPR.class;
      }

      public XMLSyntaxRule[] getSyntaxRules() {
         return rules;
      }

      private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule("weight"),
            new ElementRule(TreeModel.class) };

   };

}
