/*
 * TreeNodeSlide.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodelxml.operators.TreeNodeSlideParser;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import jebl.util.FixedBitSet;

import java.util.Arrays;

/**
 * An operator on a species tree based on the ideas of Mau et all.
 *
 * <a href="http://citeseer.ist.psu.edu/rd/27056960%2C5592%2C1%2C0.25%2CDownload/http://citeseer.ist.psu.edu/cache/papers/cs/5768/ftp:zSzzSzftp.stat.wisc.eduzSzpubzSznewtonzSzlastfinal.pdf/mau98bayesian.pdf">
 * Bayesian Phylogenetic Inference via Markov Chain Monte Carlo Methods</a>
 *
 *  @author Joseph Heled
 *         Date: 29/05/2008
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class TreeNodeSlide extends SimpleMCMCOperator {

    private final SpeciesTreeModel tree;
    private final SpeciesBindings species;

    private final int[] preOrderIndexBefore;
    private final int[] preOrderIndexAfter;
	
    private final boolean verbose = false;

   // private double range = 1.0;
  //  private boolean outgroupOnly = false;

    //private boolean verbose = false;

    public TreeNodeSlide(SpeciesTreeModel tree, SpeciesBindings species /*, boolean outgroupOnly*/, double weight) {
        this.tree = tree;
        this.species = species;
    //    this.outgroupOnly = outgroupOnly;

        preOrderIndexBefore = new int[tree.getNodeCount()];
        Arrays.fill(preOrderIndexBefore,  -1);

        preOrderIndexAfter= new int[tree.getNodeCount()];
        Arrays.fill(preOrderIndexAfter,  -1);

        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return "none";
    }

    public String getOperatorName() {
        return TreeNodeSlideParser.TREE_NODE_REHEIGHT + "(" + tree.getId() + "," + species.getId() + ")";
    }

    public double doOperation() {
        operateOneNode(0.0);
        return 0;
    }

    private void ptree(SpeciesTreeModel tree, NodeRef[] order, int[] preOrderIndex) {
        System.err.println(TreeUtils.uniqueNewick(tree, tree.getRoot()));
        System.err.println(TreeUtils.newick(tree));
        ptreenode(tree, tree.getRoot(), preOrderIndex, "");
        for(NodeRef anOrder : order) {
            System.err.print(anOrder.getNumber() + " ");
        }
        System.err.println();
    }

    private void ptreenode(SpeciesTreeModel tree, NodeRef node, int[] preOrderIndex, String indent) {
       System.err.print(indent + node.getNumber() + " " + tree.getNodeHeight(node) +
       " " + ((node != tree.getRoot()) ? tree.getParent(node).getNumber() : -1) + " ");

       if( ! tree.isExternal(node) ) {
          System.err.print("(" + tree.getChild(node, 0).getNumber() + " " + tree.getChild(node, 1).getNumber()
                  + ")");
           int k = tree.getExternalNodeCount() + 2*preOrderIndex[node.getNumber()];
           System.err.println(" [" + tree.sppSplitPopulations.getParameterValue(k) + "," +
                   tree.sppSplitPopulations.getParameterValue(k+1) + "] ");
           ptreenode(tree, tree.getChild(node, 0), preOrderIndex, indent + "  ");
           ptreenode(tree, tree.getChild(node, 1), preOrderIndex, indent + "  ");
       } else {
           System.err.println( tree.getNodeTaxon(node).toString() );
       }
    }

    public void operateOneNode(final double factor) {

//            #print "operate: tree", ut.treerep(t)
     //   if( verbose)  System.out.println("  Mau at start: " + tree.getSimpleTree());

        final int count = tree.getExternalNodeCount();  assert count == species.nSpecies();

        NodeRef[] order = new NodeRef[2 * count - 1];
        boolean[] swapped = new boolean[count-1];
        mauCanonical(tree, order, swapped);

        // internal node to change
        // count-1 - number of internal nodes
        int which = MathUtils.nextInt(count - 1);
//        if( outgroupOnly ) {
//           if( order[1] == tree.getRoot() ) {
//               which = 0;
//           } else if( order[2*count - 3] == tree.getRoot() ) {
//               which = count - 2;
//           }
//        }
//        if( verbose)  {
//            System.out.print("order:");
//            for(int k = 0 ; k < 2*count; k += 2) {
//                System.out.print(tree.getNodeTaxon(order[k]).getId() + ((k == 2*which) ? "*" : " "));
//            }
//            System.out.println();
//        }

        FixedBitSet left = new FixedBitSet(count);
        FixedBitSet right = new FixedBitSet(count);

        for(int k = 0; k < 2*which+1; k += 2) {
           left.set(tree.speciesIndex(order[k]));
        }

        for(int k = 2*(which+1); k < 2*count; k += 2) {
           right.set(tree.speciesIndex(order[k]));
        }

        double newHeight;

        if( factor > 0 ) {
          newHeight = tree.getNodeHeight(order[2*which+1]) * factor;
        } else {
          final double limit = species.speciationUpperBound(left, right);
          newHeight = MathUtils.nextDouble() * limit;
        }

        if( false ) {
            double totChange = 0;
            double tot = 0;
            double totChangeUp = 0.0;
            int topChange = 0;

            for(int which1 = 0; which1 < count - 1; which1++) {
                FixedBitSet left1 = new FixedBitSet(count);
                FixedBitSet right1 = new FixedBitSet(count);

                final NodeRef node = order[2 * which1 + 1];
                final double h = tree.getNodeHeight(node);

                for(int k = 0; k < 2 * which1 + 1; k += 2) {
                    final NodeRef ref = order[k];
                    left1.set(tree.speciesIndex(ref));
                }

                for(int k = 2 * (which1 + 1); k < 2 * count; k += 2) {
                    right1.set(tree.speciesIndex(order[k]));
                }

                final double limit = species.speciationUpperBound(left1, right1);

                final NodeRef parent = tree.getParent(node);
                final double ph = parent != null ? tree.getNodeHeight(parent) : limit;
                final double h1 = tree.getNodeHeight(tree.getChild(node, 0));
                final double h2 = tree.getNodeHeight(tree.getChild(node, 1));
                final double chDown = Math.max(h1, h2);

                final double chup = Math.max(limit - ph, 0);
                totChangeUp += chup;
                double change = chDown + chup;
                totChange += change;
                if( change > limit ) {
                    assert change <= limit;
                }

                tot += limit;

                if( which == which1 && (newHeight < chDown || newHeight > limit - chup) ) {
                    topChange = 1;
                }
//               final double h1 = tree.getNodeHeight(order[2 * which1]);
//               final double h2 = tree.getNodeHeight(order[2 * which1 + 2]);
//               double lowerLimit = 0;
//               double upperLimit = limit;
//               if( h > h1 ) {
//                   assert tree.getParent(order[2 * which]) == node;
//                   lowerLimit = h1;
//               } else {
//                   assert tree.getParent() == node;
//                   upperLimit = Math.min(upperLimit, h1);
//               }
//               if( h > h2 && h2 > lowerLimit ) {
//                   assert tree.getParent(order[2 * which+2]) == node;
//                   lowerLimit = h2;
//               }
//
            }
            final double p = totChange / tot;
            System.err.println("top% " + p + " " + totChangeUp / tot + (topChange > 0 ? "++" : ""));
        }

        tree.beginTreeEdit();

        tree.setPreorderIndices(preOrderIndexBefore);

        final NodeRef node = order[2 * which + 1];

        if( false ) {
        System.err.println("Before " + node.getNumber() + " " + newHeight);
        ptree(tree, order, preOrderIndexBefore);

        if( newHeight > tree.getNodeHeight(node) ) {
            NodeRef parent = tree.getParent(node);
            while( parent != null && tree.getNodeHeight(parent) < newHeight)  {
                // swap
                final int pi = getIndexOf(parent, order);
                final int side = pi < which ? 0 : 1;
                // 0 for right, 1 for left
                //final int side = (1+getSide(which, parent, order, swapped))/2;

                swapPops(node, swapped[which] ? 1-side : side,
                        parent, swapped[pi] ? side : 1-side, preOrderIndexBefore);
                parent = tree.getParent(parent);
            }
        } else {
          // NodeRef child = tree.getChild(node, 0);
            if( which > 0 ) {
                NodeRef nb = order[2*which - 1];
                if( tree.getNodeHeight(node) > tree.getNodeHeight(nb) ) {
                    while( nb != node && tree.getNodeHeight(nb) < newHeight ) {
                        nb = tree.getParent(nb);
                    }

                    if( nb != node ) {
                        final int side = swapped[which] ? 1 : 0;
                        swapPops(node, side, nb, 1-side, preOrderIndexBefore);

                        NodeRef pnb = tree.getParent(nb);
                        while( pnb != node ) {
                            swapPops(nb, 1-side, pnb, 1-side, preOrderIndexBefore);
                            nb = pnb;
                            pnb = tree.getParent(nb);
                        }
                    }
                }
            }
            if( which < count - 2 ) {
                NodeRef nb = order[2*which + 3];
                if( tree.getNodeHeight(node) > tree.getNodeHeight(nb) ) {
                    while( nb != node && tree.getNodeHeight(nb) < newHeight ) {
                        nb = tree.getParent(nb);
                    }
                    if( nb != node ) {
                        final int side = swapped[which] ? 0 : 1;
                        swapPops(node, side, nb, 1-side, preOrderIndexBefore);

                        NodeRef pnb = tree.getParent(nb);
                        while( pnb != node ) {
                            swapPops(nb, 1-side, pnb, 1-side, preOrderIndexBefore);
                            nb = pnb;
                            pnb = tree.getParent(nb);
                        }
                    }
                }
            }
        }
        }
        tree.setNodeHeight(node, newHeight);

        mauReconstruct(tree, order, swapped);

        // restore pre-order of pops -
        {
            tree.setPreorderIndices(preOrderIndexAfter);

            double[] splitPopValues = null;

            for(int k = 0; k < preOrderIndexBefore.length; ++k) {
                final int b = preOrderIndexBefore[k];
                if( b >= 0 ) {
                    final int a = preOrderIndexAfter[k];
                    if( a != b ) {
                        //if( verbose)  System.out.println("pops: " + a + " <- " + b);

                        final Parameter p1 = tree.sppSplitPopulations;
                        if( splitPopValues == null ) {
                            splitPopValues = p1.getParameterValues();
                        }

                        if( tree.constPopulation()  ) {
                           p1.setParameterValue(count + a, splitPopValues[count + b]);
                        } else {
                          for(int i = 0; i < 2; ++i) {
                              p1.setParameterValue(count + 2*a + i, splitPopValues[count + 2*b + i]);
                          }
                        }
                    }
                }
            }
        }

       // System.err.println("After");
       //ptree(tree, order, preOrderIndexAfter);
        
       // if( verbose)  System.out.println("  Mau after: " + tree.getSimpleTree());
       // }

        tree.endTreeEdit();
    }

    private int getSide(int which, NodeRef parent, NodeRef[] order, boolean[] swapped) {
        for(int k = 1; k < order.length; k+=2) {
            if( order[k] == parent ) {
                //if( swapped[(k-1)/2] ) s = -s;
                return k < 2*which+1 ? -1 : +1;
            }
        }
        assert false;
        return 0;
    }

     private int getIndexOf(NodeRef parent, NodeRef[] order) {
        for(int k = 1; k < order.length; k+=2) {
            if( order[k] == parent ) {
               return (k-1)/2;
            }
        }
        assert false;
        return -1;
    }

    private void swapPops(NodeRef node, int i, NodeRef parent, int i1, int[] preOrderIndexBefore) {
        int count = tree.getExternalNodeCount();
        int l1 = count + 2*preOrderIndexBefore[node.getNumber()] + i;
        int l2 = count + 2*preOrderIndexBefore[parent.getNumber()] + i1;
        final double p1 = tree.sppSplitPopulations.getParameterValue(l1);
        final double p2 = tree.sppSplitPopulations.getParameterValue(l2);
        tree.sppSplitPopulations.setParameterValue(l1, p2);
        tree.sppSplitPopulations.setParameterValue(l2, p1);
    }

//    static private void treeMixup(Tree tree, NodeRef node) {
//        final double h = tree.getNodeHeight(node);
//        if( ! tree.isRoot(node) )  {
//            assert tree.getBranchLength(node) ==  (tree.getNodeHeight(tree.getParent(node)) - h);
//        }
//        if( ! tree.isExternal(node) ) {
//           for(int k = 0; k < tree.getChildCount(node); ++k) {
//               assert  h >  tree.getNodeHeight(tree.getChild(node, k));
//               final NodeRef child = tree.getChild(node, k);
//               treeMixup(tree, child);
//           }
//        }
//    }

//    static private void setPreorderIndices(SpeciesTreeModel tree, int[] indices) {
//        setPreorderIndices(tree, tree.getRoot(), 0, indices);
//    }
//
//    static private int setPreorderIndices(SpeciesTreeModel tree, NodeRef node, int loc, int[] indices) {
//        if( ! tree.isExternal(node) ) {
//            int l = setPreorderIndices(tree, tree.getChild(node, 0), loc, indices);
//            indices[node.getNumber()] = l;
//            loc = setPreorderIndices(tree, tree.getChild(node, 1), l+1, indices);
//        }
//        return loc;
//    }

    /**
     * Obtain an ordering of tree tips from randomly swaping the children order in internal nodes.
     *
     * @param tree    tree to create order from
     * @param order   Nodes in their random order (only odd indices are filled)
     * @param wasSwapped  true if internal node was swapped
     */
    static private void mauCanonical(Tree tree, NodeRef[] order, boolean[] wasSwapped) {
        mauCanonicalSub(tree, tree.getRoot(), 0, order, wasSwapped);
    }

    static private int mauCanonicalSub(Tree tree, NodeRef node, int loc, NodeRef[] order, boolean[] wasSwaped) {
        if( tree.isExternal(node) ) {
            order[loc] = node;     assert (loc & 0x1) == 0;
            return loc + 1;
        }

        final boolean swap = MathUtils.nextBoolean();
        //wasSwaped[(loc-1)/2] = swap;
       
        int l = mauCanonicalSub(tree, tree.getChild(node, swap ? 1 : 0), loc, order, wasSwaped);

        order[l] = node;   assert (l & 0x1) == 1;
        wasSwaped[(l-1)/2] = swap;

        l = mauCanonicalSub(tree, tree.getChild(node, swap ? 0 : 1), l+1, order, wasSwaped);
        return l;
    }

    static private void mauReconstruct(SpeciesTreeModel tree, NodeRef[] order, boolean[] swapped) {
        final NodeRef root = mauReconstructSub(tree, 0, swapped.length, order, swapped);
        if( tree.getRoot() != root ) {
            tree.setRoot(root);
        }
    }

    static private NodeRef mauReconstructSub(SpeciesTreeModel tree, int from, int to, NodeRef[] order, boolean[] wasSwaped) {
        if( from == to ) {
            return order[2*from];
        }

        int rootIndex = -1;
        {
            double h = -1;

            for(int i = from; i < to; ++i) {
                final double v = tree.getNodeHeight(order[2 * i + 1]);
                if( h < v ) {
                    h = v;
                    rootIndex = i;
                }
            }
        }

        final NodeRef root = order[2 * rootIndex + 1];
        
        final NodeRef lchild = tree.getChild(root, 0);
        final NodeRef rchild = tree.getChild(root, 1);

        NodeRef lTargetChild = mauReconstructSub(tree, from, rootIndex, order, wasSwaped);
        NodeRef rTargetChild = mauReconstructSub(tree, rootIndex+1, to, order, wasSwaped);

        if( wasSwaped[rootIndex] ) {
            NodeRef z = lTargetChild;
            lTargetChild = rTargetChild;
            rTargetChild = z;
        }

        if( lchild != lTargetChild ) {
            tree.replaceChild(root, lchild, lTargetChild);
        }

        if( rchild != rTargetChild ) {
            tree.replaceChild(root, rchild, rTargetChild);
        }

        return root;
    }

}
