/**
 * 
 */
package dr.evomodel.operators;

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
 * This is an operator similar to the p-Edge-Contract-Refine operator proposed
 * by Ganeshkumar Ganapathy, Vijaya Ramachandran and Tandy Warnow.
 * 
 * @author Sebastian Hoehna
 * @version 1.0
 * 
 */
public class ECR extends SimpleMCMCOperator {

   private TreeModel _tree = null;

   public static final String ECR = "EdgeContractRefine";

   // private int _edges = 3;

   /**
    * 
    */
   public ECR(TreeModel tree, int weight) {
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

      // get a random node where the parent and the grandparent aren't the
      // root
      NodeRef i, iFather, iGrandfather;
      final NodeRef root = _tree.getRoot();
      final int nNodes = _tree.getNodeCount();

      if (nNodes < 8) {
         throw new OperatorFailedException(
               "Couldn't find valid ECR move on this tree because tree is to small!!");
      }

      do {
         i = _tree.getNode(MathUtils.nextInt(nNodes));
      } while (_tree.getChildCount(i) != 2 || root == i
            || _tree.getParent(i) == root);

      iFather = _tree.getParent(i);
      iGrandfather = _tree.getParent(iFather);
      NodeRef[] nodes = { _tree.getChild(i, 0), _tree.getChild(i, 1),
            getOtherChild(_tree, iFather, i),
            getOtherChild(_tree, iGrandfather, iFather) };

      int[] positions = MathUtils.shuffled(nodes.length);

      _tree.beginTreeEdit();

      _tree.removeChild(i, nodes[0]);
      _tree.removeChild(i, nodes[1]);
      _tree.removeChild(iFather, nodes[2]);
      _tree.removeChild(iGrandfather, nodes[3]);

      double nodeHeightGrandfather = _tree.getNodeHeight(iGrandfather);
      double minNodeHeight = Math.max(_tree
            .getNodeHeight(nodes[positions[0]]), _tree
            .getNodeHeight(nodes[positions[1]]));
      double oldMinNodeHeight = Math.max(_tree
            .getNodeHeight(nodes[0]), _tree
            .getNodeHeight(nodes[1]));
      double newNodeHeight = (MathUtils.nextDouble()
            * (nodeHeightGrandfather - minNodeHeight) + minNodeHeight);
      
      _tree.setNodeHeight(i, newNodeHeight);
      _tree.addChild(i, nodes[positions[0]]);
      _tree.addChild(i, nodes[positions[1]]);

      double hastingsRatio = (nodeHeightGrandfather - minNodeHeight)/(nodeHeightGrandfather - oldMinNodeHeight);
      
      minNodeHeight = Math.max(newNodeHeight, _tree
            .getNodeHeight(nodes[positions[2]]));
      newNodeHeight = (MathUtils.nextDouble()
            * (nodeHeightGrandfather - minNodeHeight) + minNodeHeight);
      _tree.setNodeHeight(iFather, newNodeHeight);
      _tree.addChild(iFather, nodes[positions[2]]);
      
      oldMinNodeHeight = Math.max(_tree
            .getNodeHeight(i), _tree
            .getNodeHeight(nodes[2]));
      hastingsRatio *= (nodeHeightGrandfather - minNodeHeight)/(nodeHeightGrandfather - oldMinNodeHeight);
      
      _tree.addChild(iGrandfather, nodes[positions[3]]);

      try {
         _tree.endTreeEdit();
      } catch (MutableTree.InvalidTreeException ite) {
         throw new OperatorFailedException(ite.toString());
      }

      // TODO
      _tree.pushTreeChangedEvent(iFather);
      _tree.pushTreeChangedEvent(iGrandfather);

      return Math.log(Math.min(1, hastingsRatio));
   }

   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
    */
   @Override
   public String getOperatorName() {
      return ECR;
   }

   public double getMinimumAcceptanceLevel() {
      return 0.025;
   }

   public double getMinimumGoodAcceptanceLevel() {
      return 0.05;
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

   public static XMLObjectParser ECR_PARSER = new AbstractXMLObjectParser() {

      public String getParserName() {
         return ECR;
      }

      public Object parseXMLObject(XMLObject xo) throws XMLParseException {

         TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
         int weight = xo.getIntegerAttribute("weight");

         return new ECR(treeModel, weight);
      }

      // ************************************************************************
      // AbstractXMLObjectParser implementation
      // ************************************************************************

      public String getParserDescription() {
         // TODO
         return "This element represents a ECR operator. "
               + "This operator does something similar to NNI.";
      }

      public Class getReturnType() {
         return ECR.class;
      }

      public XMLSyntaxRule[] getSyntaxRules() {
         return rules;
      }

      private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule("weight"),
            new ElementRule(TreeModel.class) };

   };

}
