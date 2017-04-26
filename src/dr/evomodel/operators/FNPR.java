/*
 * FNPR.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

/**
 *
 */
package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.FNPRParser;
import dr.math.MathUtils;

/**
 * This is an implementation of the Subtree Prune and Regraft (SPR) operator for
 * trees. It assumes explicitely bifurcating rooted trees.
 *
 * @author Sebastian Hoehna
 * @version 1.0
 */
public class FNPR extends AbstractTreeOperator {

    private TreeModel tree = null;
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
    public double doOperation() {
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

              tree.endTreeEdit();

              tree.pushTreeChangedEvent(i);

              assert tree.getExternalNodeCount() == tipCount;

              return 0.0;
           }
        }

        throw new RuntimeException("Couldn't find valid SPR move on this tree!");
     }

    /*
    * (non-Javadoc)
    *
    * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
    */
    @Override
    public String getOperatorName() {
        return FNPRParser.FNPR;
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

}
