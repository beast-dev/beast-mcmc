/*
 * TreeUniform.java
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

package dr.evomodel.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.TreeUniformParser;
import dr.math.MathUtils;

/**
 *  Uniform height change of several related tree nodes at once.
 *
 *  Currently supported - a parent-son pair and parent-two sons triplet.
 *
 *  @author  Joseph Heled
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class TreeUniform extends AbstractTreeOperator {

    private final int nodesToMove;

    //private static final int MAX_TRIES = 100;
    private final TreeModel tree;

    public TreeUniform(int n, TreeModel tree, double weight) {
        assert n == 2 || n == 3;

        this.tree = tree;
        this.nodesToMove = n;
        setWeight(weight);
    }

    public double doOperation() {
        if( tree.getInternalNodeCount() < 2 ) {
            throw new RuntimeException("no node found");
        }

        tree.beginTreeEdit();

        switch( nodesToMove ) {
            case 2: move2(); break;
            case 3: move3(); break;
        }

        tree.endTreeEdit();

        // AR not sure whether this check is needed...
        try {
            tree.checkTreeIsValid();
        } catch( MutableTree.InvalidTreeException ite ) {
            throw new RuntimeException(ite.toString());
        }

        return 0;
    }


    // insure two elements are sorted in reverse order
    private void sorted(double[] h, int i, int j) {
        if( h[i] < h[j] ) {
           exch(h, i, j);
        }
    }

    private void exch(double[] h, int i, int j) {
        final double t = h[i];
        h[i] = h[j];
        h[j] = t;
    }

    private void move3() /*throws OperatorFailedException*/ {
        final int nInternalNodes = tree.getInternalNodeCount();

        final NodeRef root = tree.getRoot();

        //  find an internal node with 2 internal children
        NodeRef i = root;
        int nTries = 0;

        while( root == i ) {
            i = tree.getInternalNode(MathUtils.nextInt(nInternalNodes));
        }
        boolean ise0 = tree.isExternal(tree.getChild(i, 0));
        boolean ise1 = tree.isExternal(tree.getChild(i, 1));

        if( ise0 || ise1 ) {
            if( ise0 != ise1 ) {
                doMove2(tree.getChild(i, ise0 ? 1 : 0) );
            } else {
                doMove1(i);
            }
            return;
        }

        final NodeRef iParent = tree.getParent(i);

        final NodeRef child0 = tree.getChild(i, 0);
        final NodeRef child1 = tree.getChild(i, 1);

        final double hMin0 = getMaxChildHeight(child0);
        final double hMin1 = getMaxChildHeight(child1);
        // upper limit for i
        final double hMax = tree.getNodeHeight(iParent);

        final double mx = Math.max(hMin0, hMin1);
        final double mn = Math.min(hMin0, hMin1);

        final double d1 = hMax - mx;

        // factor of d1^2 dropped
        final double
                a0 = Math.abs(hMin0 - hMin1) / 2,
                a1 = d1 / 3;

        final double th = a0 / (a0 + a1);

        double[] h = new double[3];

        final int il = (hMin0 < hMin1) ? 0 : 1;

        for(int k = 0; k < 3; ++k) {
            h[k] = MathUtils.uniform(mx, hMax);
        }

        if( MathUtils.nextDouble() < th ) {
            h[1+il] = MathUtils.uniform(mn,mx);
            sorted(h, 0, 2-il);
        } else {
            int iMax =  h[0] > h[1] ? 0 : 1;
            iMax = h[iMax] < h[2] ? 2 : iMax;
            exch(h, 0, iMax);
        }

        assert hMax > h[0] && h[0] > h[1] && h[0] > h[2] && h[1] > hMin0 && h[2] > hMin1;

        NodeRef[] nodes = {i, child0, child1};
        for(int k = 0; k < nodes.length; ++k) {
            tree.setNodeHeight(nodes[k], h[k]);
            tree.pushTreeChangedEvent(nodes[k]);
        }
    }

    private void move2() {
        final int nInternalNodes = tree.getInternalNodeCount();

        final NodeRef root = tree.getRoot();

        NodeRef i = root;

        while( root == i ) {
            i = tree.getInternalNode(MathUtils.nextInt(nInternalNodes)) ;
        }

        doMove2(i);
     //   return true;
    }

   private void doMove2(NodeRef i) {
       final NodeRef iParent = tree.getParent(i);

       if( iParent == tree.getRoot() ) {
           doMove1(i);
           return;
       }

        final NodeRef iGrandParent = tree.getParent(iParent);

        // lower limit for node (max height of children)
        final double hMin = getMaxChildHeight(i);
        // upper limit for parent (it's parent height)
        final double hMax = tree.getNodeHeight(iGrandParent);
        // lower bound for parent
        final double hLim0 = Math.max(hMin, tree.getNodeHeight(getOtherChild(iParent, i)));

        final double d1 = hMax - hLim0;
        // factor of d1 dropped
        final double a0 = (hLim0 - hMin),
                     a1 = d1/2;

        final double th = a0 / (a0 + a1);
        double[] h = new double[2];

        if( MathUtils.nextDouble() < th ) {
           h[0] = MathUtils.uniform(hMin, hLim0);
           h[1] = MathUtils.uniform(hLim0, hMax);
        } else {
            for(int k = 0; k < 2; ++k) {
              h[k] = MathUtils.uniform(hLim0, hMax);
            }
            // sorted(h, 1, 0) ??
            if( h[0] > h[1] ) {
                double t = h[0];
                h[0] = h[1];
                h[1] = t;
            }
        }

        tree.setNodeHeight(iParent, h[1]);
        tree.setNodeHeight(i, h[0]);

        tree.pushTreeChangedEvent(iParent);
        tree.pushTreeChangedEvent(i);
    }

    private void doMove1(NodeRef i) {
        final NodeRef iParent = tree.getParent(i);
        final double hMax = tree.getNodeHeight(iParent);

        // lower limit for node (max height of children)
        final double hMin = getMaxChildHeight(i);
        double h = MathUtils.uniform(hMin, hMax);
        tree.setNodeHeight(i, h);

        tree.pushTreeChangedEvent(i);
    }

    private NodeRef getOtherChild(NodeRef n, NodeRef c) {
        for(int nc = 0; nc < tree.getChildCount(n); ++nc ) {
            final NodeRef child = tree.getChild(n, nc);
            if( child != c ) {
               return child;
           }
        }
        assert true;
        return null;
    }

    private double getMaxChildHeight(NodeRef n) {
        double hMax = -1;
        for(int nc = 0; nc < tree.getChildCount(n); ++nc ) {
            hMax = Math.max(hMax, tree.getNodeHeight(tree.getChild(n, nc) ));
        }
        return hMax;
    }


    public String getPerformanceSuggestion() {
        return "";
    }

    public String getOperatorName() {
        return TreeUniformParser.TREE_UNIFORM + "(" + nodesToMove + "," + tree.getId() + ")";
    }

}