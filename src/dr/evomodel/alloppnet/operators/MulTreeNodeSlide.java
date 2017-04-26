/*
 * MulTreeNodeSlide.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.alloppnet.speciation.MulSpeciesTreeModel;
import dr.evomodel.alloppnet.speciation.MulSpeciesBindings;
import dr.evomodelxml.operators.TreeNodeSlideParser;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import jebl.util.FixedBitSet;

import java.util.Arrays;

/**
 * An operator on a multiply labelled species tree based on the ideas of Mau et all. Very similar to
 * JH's TreeNodeSlide.
 *
 * <a href="http://citeseer.ist.psu.edu/rd/27056960%2C5592%2C1%2C0.25%2CDownload/http://citeseer.ist.psu.edu/cache/papers/cs/5768/ftp:zSzzSzftp.stat.wisc.eduzSzpubzSznewtonzSzlastfinal.pdf/mau98bayesian.pdf">
 * Bayesian Phylogenetic Inference via Markov Chain Monte Carlo Methods</a>
 *
 *  @author Joseph Heled, Graham Jones
 *         Date: 21/12/2011
 */
public class MulTreeNodeSlide extends SimpleMCMCOperator {

    private final MulSpeciesTreeModel multree;
    private final MulSpeciesBindings species;

    private final int[] preOrderIndexBefore;
    private final int[] preOrderIndexAfter;

    public MulTreeNodeSlide(MulSpeciesTreeModel tree, MulSpeciesBindings species /*, boolean outgroupOnly*/, double weight) {
        this.multree = tree;
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
        return TreeNodeSlideParser.TREE_NODE_REHEIGHT + "(" + multree.getId() + "," + species.getId() + ")";
    }

    public double doOperation() {
        operateOneNode(0.0);
        return 0;
    }



    public void operateOneNode(final double factor) {

//            #print "operate: tree", ut.treerep(t)
     //   if( verbose)  System.out.println("  Mau at start: " + tree.getSimpleTree());

        final int count = multree.getExternalNodeCount();  
        assert count == species.nSpSeqs();

        NodeRef[] order = new NodeRef[2 * count - 1];
        boolean[] swapped = new boolean[count-1];
        mauCanonical(multree, order, swapped);

        // internal node to change
        // count-1 - number of internal nodes
        int which = MathUtils.nextInt(count - 1);

        FixedBitSet left = new FixedBitSet(count);
        FixedBitSet right = new FixedBitSet(count);

        for(int k = 0; k < 2*which+1; k += 2) {
           left.set(multree.speciesIndex(order[k]));
        }

        for(int k = 2*(which+1); k < 2*count; k += 2) {
           right.set(multree.speciesIndex(order[k]));
        }

        double newHeight;

        if( factor > 0 ) {
          newHeight = multree.getNodeHeight(order[2*which+1]) * factor;
        } else {
          final double limit = species.speciationUpperBound(left, right);
          newHeight = MathUtils.nextDouble() * limit;
        }

        multree.beginTreeEdit();

        multree.setPreorderIndices(preOrderIndexBefore);

        final NodeRef node = order[2 * which + 1];

        multree.setNodeHeight(node, newHeight);

        mauReconstruct(multree, order, swapped);

        // restore pre-order of pops -
        {
            multree.setPreorderIndices(preOrderIndexAfter);

            double[] splitPopValues = null;

            for(int k = 0; k < preOrderIndexBefore.length; ++k) {
                final int b = preOrderIndexBefore[k];
                if( b >= 0 ) {
                    final int a = preOrderIndexAfter[k];
                    if( a != b ) {
                        //if( verbose)  System.out.println("pops: " + a + " <- " + b);

                        final Parameter p1 = multree.sppSplitPopulations;
                        if( splitPopValues == null ) {
                            splitPopValues = p1.getParameterValues();
                        }

                        if( multree.constPopulation()  ) {
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

        multree.endTreeEdit();
    }


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

    static private void mauReconstruct(MulSpeciesTreeModel tree, NodeRef[] order, boolean[] swapped) {
        final NodeRef root = mauReconstructSub(tree, 0, swapped.length, order, swapped);
        if( tree.getRoot() != root ) {
            tree.setRoot(root);
        }
    }

    static private NodeRef mauReconstructSub(MulSpeciesTreeModel tree, int from, int to, NodeRef[] order, boolean[] wasSwaped) {
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
