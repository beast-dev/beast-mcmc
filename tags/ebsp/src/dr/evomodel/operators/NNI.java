/**
 * 
 */
package dr.evomodel.operators;

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
 * Implements the Nearest Neighbor Interchange (NNI) operation. This particular
 * implementation assumes explicitely bifurcating trees. It works similar to the
 * Narrow Exchange but with manipulating the height of a node if necessary.
 * 
 * @author Sebastian Hoehna
 * @version 1.0
 * 
 */
public class NNI extends SimpleMCMCOperator {

   private TreeModel          _tree = null;

   public static final String NNI   = "NearestNeighborInterchange";

   /**
    * 
    */
   public NNI(TreeModel tree, int weight) {
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
      double hastingsRatio = 0;

      int tipCount = _tree.getExternalNodeCount();

      final int nNodes = _tree.getNodeCount();
      final NodeRef root = _tree.getRoot();

      NodeRef i;

      // get a random node where neither you or your father is the root
      do {
         i = _tree.getNode(MathUtils.nextInt(nNodes));
      } while (root == i || _tree.getParent(i) == root);

      // get parent node
      final NodeRef iParent = _tree.getParent(i);
      // get parent of parent -> grant parent :)
      final NodeRef iGrandParent = _tree.getParent(iParent);
      // get left child of grant parent -> uncle
      NodeRef iUncle = _tree.getChild(iGrandParent, 0);
      // check if uncle == father
      if (iUncle == iParent) {
         // if so take right child -> sibling of father
         iUncle = _tree.getChild(iGrandParent, 1);
      }

      // check if height of me smaller than height of my grant father
      // (isn't it obvious?!?)
      assert _tree.getNodeHeight(i) < _tree.getNodeHeight(iGrandParent);

      // change the height of my father to be randomly between my uncle's
      // heights and my grandfather's height
      // this is necessary for the hastings ratio to do also if the uncle
      // is
      // younger anyway

      double heightGrandfather = _tree.getNodeHeight(iGrandParent);
      double heightUncle = _tree.getNodeHeight(iUncle);
      double minHeightFather = Math.max(heightUncle, _tree.getNodeHeight(getOtherChild(_tree, iParent, i)));
      double heightI = _tree.getNodeHeight(i);
      double minHeightReverse = Math.max(heightI, _tree.getNodeHeight(getOtherChild(_tree, iParent, i)));

      double ran;
      do {
         ran = Math.random();
      } while (ran == 0.0 || ran == 1.0);

      // now calculate the new height for father between the height of the
      // uncle and the grandparent
      double newHeightFather = minHeightFather
            + (ran * (heightGrandfather - minHeightFather));
      // set the new height for the father
      _tree.setNodeHeight(iParent, newHeightFather);

      // double prForward = 1 / (heightGrandfather - minHeightFather);
      // double prBackward = 1 / (heightGrandfather - minHeightReverse);
      // hastings ratio = backward Prob / forward Prob
      hastingsRatio = Math.log(Math.min(1,
            (heightGrandfather - minHeightFather)
                  / (heightGrandfather - minHeightReverse)));
      // now change the nodes
      exchangeNodes(_tree, i, iUncle, iParent, iGrandParent);

      _tree.pushTreeChangedEvent(iParent);
      _tree.pushTreeChangedEvent(iGrandParent);
      // throw new OperatorFailedException(
      // "Couldn't find valid narrow move on this tree!!");
      // maybe instead of a new try the height of the father can be
      // increment to be then between the uncle and the grant father
      // System.out.println("tries = " + tries);

      if (_tree.getExternalNodeCount() != tipCount) {
         throw new RuntimeException("Lost some tips in NNI operation");
      }

      return hastingsRatio;
   }

   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
    */
   @Override
   public String getOperatorName() {
      return NNI;
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

   public static XMLObjectParser NNI_PARSER = new AbstractXMLObjectParser() {

                                               public String getParserName() {
                                                  return NNI;
                                               }

                                               public Object parseXMLObject(
                                                     XMLObject xo)
                                                     throws XMLParseException {

                                                  TreeModel treeModel = (TreeModel) xo
                                                        .getChild(TreeModel.class);
                                                  int weight = xo
                                                        .getIntegerAttribute("weight");

                                                  return new NNI(treeModel,
                                                        weight);
                                               }

                                               // ************************************************************************
                                               // AbstractXMLObjectParser
                                                // implementation
                                               // ************************************************************************

                                               public String getParserDescription() {
                                                  return "This element represents a NNI operator. "
                                                        + "This operator swaps a random subtree with its uncle.";
                                               }

                                               public Class getReturnType() {
                                                  return NNI.class;
                                               }

                                               public XMLSyntaxRule[] getSyntaxRules() {
                                                  return rules;
                                               }

                                               private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
         AttributeRule.newIntegerRule("weight"),
         new ElementRule(TreeModel.class)                                   };

                                            };

}
