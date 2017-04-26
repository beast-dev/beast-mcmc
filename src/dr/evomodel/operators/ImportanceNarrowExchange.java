/*
 * ImportanceNarrowExchange.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Joseph Heled
 * @version 1.0
 */
@SuppressWarnings({"ConstantConditions"})
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class ImportanceNarrowExchange extends AbstractTreeOperator implements TreeLogger.LogUpon {

    private TreeModel tree = null;

    private final double epsilon;

    private int[] nodeCounts;
    private boolean justAccepted;

    private final double[] weights;
    private double totalWeight;

    public ImportanceNarrowExchange(TreeModel tree, PatternList patterns, double epsilon, double weight) throws Exception {
        this.tree = tree;
        setWeight(weight);

        justAccepted = false;
        this.epsilon = epsilon;

        weights = new double[tree.getNodeCount()];
        setTaxaWeights(patterns);
    }

    private void setTaxaWeights(PatternList patterns) throws Exception {
        final DataType type = patterns.getDataType();
        Map<Integer, Integer> counts = new HashMap<Integer, Integer>();

        int[] taxaCounts = new int[patterns.getPatternLength()];

        for(int nPat = 0; nPat < patterns.getPatternCount(); ++nPat) {
            final int[] pattern = patterns.getPattern(nPat);

            counts.clear();
            for( int s : pattern ) {
                if( type.isGapState(s) || type.isAmbiguousState(s) || type.isUnknownState(s) ) {
                    continue;
                }
                if( ! counts.containsKey(s) ) {
                   counts.put(s, 0);
                }
                counts.put(s, counts.get(s)+1);
            }

            if( counts.size() <= 1 ) {
                continue;
            }

            Map.Entry<Integer, Integer> m = null;
            for( Map.Entry<Integer, Integer> e : counts.entrySet()) {
               if( m == null || e.getValue() > m.getValue() ) {
                   m = e;
               }
            }
            assert m != null;

            for(int i = 0; i < pattern.length; ++i) {
                final int s = pattern[i];
                if( ! (type.isGapState(s) || type.isAmbiguousState(s) || type.isUnknownState(s) ) ) {
                    if( s != m.getKey() ) {
                        taxaCounts[i] += patterns.getPatternWeight(nPat);
                    }
                }
            }
        }

        nodeCounts = new int[tree.getNodeCount()];
        Map<Taxon, Integer> taxaWeights = new HashMap<Taxon, Integer>();
        for(int i = 0; i < taxaCounts.length; ++i) {
            taxaWeights.put(patterns.getTaxon(i), taxaCounts[i]);
        }

        for(int i = 0; i < tree.getExternalNodeCount(); ++i) {
            final NodeRef leaf = tree.getExternalNode(i);
            final Taxon nodeTaxon = tree.getNodeTaxon(leaf);

//            assert taxaWeights.containsKey(nodeTaxon) : nodeTaxon;
            if( ! taxaWeights.containsKey(nodeTaxon) ) {
                throw new  Exception("" + nodeTaxon + " in tree " + tree.getId() +
                " not in patterns" + patterns.getId() + ".");
            }
            nodeCounts[leaf.getNumber()] = taxaWeights.get(nodeTaxon) ;
        }
    }


    private int traverseTree(NodeRef n) {
        final int k = n.getNumber();
        if( ! tree.isExternal(n) ) {
            int w = 0;
            for(int nc = 0; nc < tree.getChildCount(n); ++nc ) {
                w += traverseTree(tree.getChild(n, nc));
            }
            nodeCounts[k] = w;
        }

        return nodeCounts[k];
    }

    final private int DEBUG = 0;

    private double nodeWeight(final NodeRef node ) {
        final NodeRef ch0 = tree.getChild(node, 0);
        final NodeRef ch1 = tree.getChild(node, 1);

        if( tree.isExternal(ch0) &&  tree.isExternal(ch1) ) {
            return 0;
        }

        final boolean leftSubtree = tree.getNodeHeight(ch0) < tree.getNodeHeight(ch1);
        final int st0 = nodeCounts[(leftSubtree ? ch0 : ch1).getNumber()];
        final int st1 = nodeCounts[tree.getChild(leftSubtree ? ch1 : ch0, 0).getNumber()];
        final int st2 = nodeCounts[tree.getChild(leftSubtree ? ch1 : ch0, 1).getNumber()];
        final double w = (epsilon + st0)*(epsilon + st1) + (epsilon + st0)*(epsilon + st2) + (epsilon + st1)*(epsilon + st2)
                - 3*epsilon*epsilon;
        return w;
    }

    private int getNode() {
        traverseTree(tree.getRoot());

        totalWeight = 0;

        for(int k = 0; k < tree.getInternalNodeCount(); ++k) {
            final NodeRef node = tree.getInternalNode(k);
            final double w = nodeWeight(node);

            weights[node.getNumber()] = w;
            if( DEBUG > 5 && w > 0 ) {
              System.out.println("" + w + " " + TreeUtils.uniqueNewick(tree, node));
            }
            totalWeight += w;
        }

        double r = MathUtils.nextDouble() * totalWeight;
        for(int k = 0; k < tree.getInternalNodeCount(); ++k) {
            final NodeRef node = tree.getInternalNode(k);
            final int nodeIndex = node.getNumber();
            r -= weights[nodeIndex];
            if( r < 0 ) {
                if( DEBUG > 0 ) {
                    System.out.println("" + weights[nodeIndex] + "/" + totalWeight + " " + TreeUtils.uniqueNewick(tree, node));
                }
                return k;
            }
        }
        //assert false;
        return -1;
    }

    /*
     * (non-Javadoc)
     *
     * @see dr.inference.operators.SimpleMCMCOperator#doOperation()
     */
    @Override
    public double doOperation() {
        int k = getNode();
        if( k < 0 ) {
            throw new RuntimeException("no node found");
        }

        final NodeRef p = tree.getInternalNode(k);
        if( DEBUG > 0 ) {
            System.out.println(TreeUtils.newick(tree));
            System.out.println("" + getAcceptCount() + " - " + getRejectCount());
        }
        assert tree.getChildCount(p) == 2;
        final NodeRef ch0 = tree.getChild(p, 0);
        final NodeRef ch1 = tree.getChild(p, 1);

        final boolean side = tree.getNodeHeight(ch0) < tree.getNodeHeight(ch1);
        final NodeRef iUncle = side ? ch0 : ch1;
        final NodeRef jP = side ? ch1 : ch0;
        final NodeRef j = tree.getChild(jP, MathUtils.nextInt(2));
        exchangeNodes(tree, iUncle, j, p, jP);

        final int jPindex = jP.getNumber();
        nodeCounts[jPindex] += -nodeCounts[j.getNumber()] + nodeCounts[iUncle.getNumber()];         // debug

        // the weights function is symmetric with respect to counts of the three sub-trees involved in the exchange,
        // so the weight of the root node (p) does not change, but the weight of jP has changed

        final double prev = weights[jPindex];
        // The counts below jP are still valid
        final double now = nodeWeight(jP);
        double newTot = totalWeight + (now - prev);

        weights[jPindex] = now;  // debug
        final NodeRef pP = tree.getParent(p);
        if( pP != null ) {
            final int pPindex = pP.getNumber();
            final double prev1 = weights[pPindex];
            final double now1 = nodeWeight(pP);
            weights[pPindex] = now1;    // debug
            newTot += + (now1 - prev1);
        }

        // pr(node before operation) = w/tot. pr(node after operation) = w/newTot
        // log(back/forward) = w/newTot / w/tot = log(tot/newTot)

        double saveTotalWeight = totalWeight;        // debug
        double[] w = new double[weights.length];     // debug

        System.arraycopy(weights, 0, w, 0, w.length);  // debug
        int[] c = new int[nodeCounts.length];          // debug
        System.arraycopy(nodeCounts, 0, c, 0, c.length); // debug
        getNode();                                     // debug

        for(int l = 0; l < c.length; ++l)  {         // debug
            if( c[l] != nodeCounts[l] ) {            // debug
                assert false;                       // debug
            }                                       // debug
        }                                          // debug
        for(int l = 0; l < w.length; ++l)  {                  // debug
            if( Math.abs(weights[l]/w[l] - 1) > 1e-12 ) {     // debug
                assert false;                                 // debug
            }                                                  // debug
        }                                                     // debug
        assert Math.abs(newTot/totalWeight  - 1) < 1e-10;    // debug

        return Math.log(saveTotalWeight / newTot);
    }

    public void reject() {
        super.reject();
        justAccepted = false;
    }
    
    public void accept(double deviation) {
        super.accept(deviation);
        justAccepted = true;
    }

    private final long lFreq = 1000;
    private long lastLog = -lFreq-1;

    public boolean logNow(long state) {

        boolean r = justAccepted;
        if( lastLog + lFreq >= state ) {
            r = false;
        } else if( r ) {
            lastLog = state;
        }
        justAccepted = false;
        return r;       
    }

    /*
     * (non-Javadoc)
     *
     * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
     */
    @Override

    public String getOperatorName() {
        return  "Importance Narrow Exchange" + "(" + tree.getId() + ")";
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
        return "";
    }

}