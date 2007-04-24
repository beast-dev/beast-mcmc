/*
 * CoalescentSimulator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.coalescent;

import dr.evolution.tree.*;
import dr.evolution.util.*;
import dr.xml.*;
import dr.inference.distribution.ParametricDistributionModel;
import dr.math.UnivariateFunction;

import java.util.ArrayList;
import java.util.List;


/**
 * Simulates a set of coalescent intervals given a demographic model.
 *
 * @author Alexei Drummond
 * @version $Id: CoalescentSimulator.java,v 1.43 2005/10/27 10:40:48 rambaut Exp $
 */
public class CoalescentSimulator {

    public static final String COALESCENT_TREE = "coalescentTree";
    public static final String COALESCENT_SIMULATOR = "coalescentSimulator";
    public static final String RESCALE_HEIGHT = "rescaleHeight";
    public static final String ROOT_HEIGHT = "rootHeight";
    public static final String CONSTRAINED_TAXA = "constrainedTaxa";
    public static final String TMRCA_CONSTRAINT = "tmrca";

    /**
     * Simulates a coalescent tree from a set of subtrees.
     */
    public CoalescentSimulator() {
    }

    /**
     * Simulates a coalescent tree from a set of subtrees.
     * @param subtrees an array of tree to be used as subtrees
     * @param model the demographic model to use
     * @param rootHeight an optional root height with which to scale the whole tree
     */
    public SimpleTree simulateTree(Tree[] subtrees, DemographicModel model, double rootHeight) {

        SimpleNode[] roots = new SimpleNode[subtrees.length];
        SimpleTree tree = null;

        dr.evolution.util.Date mostRecent = null;
        for (int i = 0; i < subtrees.length; i++) {
            dr.evolution.util.Date date = Tree.Utils.findMostRecentDate(subtrees[i]);
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < subtrees.length; i++) {
                dr.evolution.util.Date date = Tree.Utils.findMostRecentDate(subtrees[i]);
                if (date != null) {
                    double diff = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    for (int j = 0; j < subtrees[i].getNodeCount(); j++) {
                        NodeRef node = subtrees[i].getNode(j);

/*						if (subtrees[i].isExternal(node)) {
							System.out.print(subtrees[i].getNodeTaxon(node).getId() + " - ");
							System.out.print(subtrees[i].getNodeTaxon(node).getAttribute("date") + " - ");
							System.out.print("Old height: "+Double.toString(subtrees[i].getNodeHeight(node)));
							System.out.println(" New height: "+Double.toString(subtrees[i].getNodeHeight(node) + diff));
						}
*/
                        ((SimpleTree)subtrees[i]).setNodeHeight(node, subtrees[i].getNodeHeight(node) + diff);

                    }
                }
            }
        }

        for (int i =0; i < roots.length; i++) {
            roots[i] = new SimpleNode(subtrees[i], subtrees[i].getRoot());
        }
        // if just one taxonList then finished
        if (roots.length == 1) {
            tree = new SimpleTree(roots[0]);
        } else {
            tree = new SimpleTree(simulator.simulateCoalescent(roots, model.getDemographicFunction()));
        }

        if (rootHeight > 0.0) attemptToScaleTree(tree, rootHeight);

        return tree;
    }


    /**
     * Simulates a coalescent tree, given a taxon list.
     * @param taxa the set of taxa to simulate a coalescent tree between
     * @param model the demographic model to use
     */
    public SimpleTree simulateTree(TaxonList taxa, DemographicModel model) {

        return simulator.simulateTree(taxa, model.getDemographicFunction());
    }

    private static void attemptToScaleTree(MutableTree tree, double rootHeight) {
        double scale = rootHeight / tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef n = tree.getInternalNode(i);
            tree.setNodeHeight(n, tree.getNodeHeight(n) * scale);
        }
        MutableTree.Utils.correctHeightsForTips(tree);
    }

    static class TaxaConstraint {
        final TaxonList taxons;
        final double lower;
        double upper;

        TaxaConstraint(TaxonList taxons, ParametricDistributionModel p) {
            this.taxons = taxons;
            // no constraint
            //upper = -1;

            if( p != null ) {
                final UnivariateFunction univariateFunction = p.getProbabilityDensityFunction();
                lower = univariateFunction.getLowerBound();
                upper = univariateFunction.getUpperBound();
//                if( lower == 0 && upper == Double.POSITIVE_INFINITY ) {
//                     upper = -1;
//                }
            } else {
                lower = 0;
                upper = Double.POSITIVE_INFINITY;
            }
        }

        public boolean realLimits() {
            return lower != 0 || upper !=  Double.POSITIVE_INFINITY;
        }
    }

    static private int sizeOfIntersection(TaxonList tl1, TaxonList tl2) {
        int nIn = 0;
        for(int j = 0; j < tl1.getTaxonCount(); ++j) {
            if( tl2.getTaxonIndex(tl1.getTaxon(j)) >= 0 ) {
                ++ nIn;
            }
        }
        return nIn;
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return COALESCENT_TREE; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoalescentSimulator simulator = new CoalescentSimulator();

            DemographicModel demoModel = (DemographicModel)xo.getChild(DemographicModel.class);
            ArrayList taxonLists = new ArrayList();
            ArrayList subtrees = new ArrayList();

            double rootHeight = -1;
            if (xo.hasAttribute(ROOT_HEIGHT)) {
                rootHeight = xo.getDoubleAttribute(ROOT_HEIGHT);
            }
            if (xo.hasAttribute(RESCALE_HEIGHT)) {
                rootHeight = xo.getDoubleAttribute(RESCALE_HEIGHT);
            }

            // should have one child that is node
            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);

                // AER - swapped the order of these round because Trees are TaxonLists...
                if (child instanceof Tree) {
                    subtrees.add((Tree)child);
                } else if (child instanceof TaxonList) {
                    taxonLists.add((TaxonList)child);
                } else if( xo.getChildName(i).equals(CONSTRAINED_TAXA) ) {
                    rootHeight = -1; // ignore it? should we errror?

                    XMLObject constrainedTaxa = (XMLObject) child;

                    TaxonList taxa = (TaxonList)constrainedTaxa.getChild(TaxonList.class);
                    List allc = new ArrayList();

                    for(int nc = 0; nc < constrainedTaxa.getChildCount(); ++nc) {

                        Object object = constrainedTaxa.getChild(nc);
                        if( object instanceof XMLObject ) {
                            XMLObject constraint = (XMLObject)object;

                            if( constraint.getName().equals(TMRCA_CONSTRAINT)) {
                                TaxonList taxaSubSet = (TaxonList)constraint.getChild(TaxonList.class);
                                ParametricDistributionModel dist =
                                        (ParametricDistributionModel) constraint.getChild(ParametricDistributionModel.class);

                                allc.add(new TaxaConstraint(taxaSubSet, dist));
                            }
                        }
                    }

                    if( allc.size() == 0 ) {
                        taxonLists.add(taxa);
                    } else {

                        // collect subtrees here
                        List st = new ArrayList();

                        final String setsNotCOmpatibleMessage = "taxa sets not compatible";

                        while( allc.size() > 0 ) {
                            // pick a group of taxon-subsets where each is contained in the next
                            // ordered by set inclusion from first to last

                            List next = new ArrayList();

                            // arbitrarily pick the first of remaining ones
                            next.add(allc.remove(0));
                            TaxonList baseConstraint = ((TaxaConstraint)next.get(0)).taxons;

                            for(int k = 0; k < allc.size(); ++k) {
                                final TaxonList taxonsk = ((TaxaConstraint)allc.get(k)).taxons;
                                final int nIn = sizeOfIntersection(taxonsk, baseConstraint);

                                if( nIn > 0 ) {
                                    if( nIn == baseConstraint.getTaxonCount() || nIn == taxonsk.getTaxonCount() ) {
                                        for(int j = 0; j < next.size() ; ++j) {
                                            TaxonList jtaxons = ((TaxaConstraint)next.get(j)).taxons;
                                            int c = sizeOfIntersection(jtaxons, taxonsk);
                                            if( c == taxonsk.getTaxonCount() ) {
                                                next.add(j, allc.remove(k));
                                                break;
                                            } else if( c != jtaxons.getTaxonCount() ) {
                                               throw new XMLParseException(setsNotCOmpatibleMessage);
                                            } else if( j+1 == next.size() ) {
                                                next.add(allc.remove(k));
                                                break;
                                            }
                                        }
                                        baseConstraint = ((TaxaConstraint)next.get(0)).taxons;

                                    }  else {
                                        throw new XMLParseException(setsNotCOmpatibleMessage);
                                    }
                                    --k;
                                }
                            }

                            for(int k = 1; k < next.size(); ++k) {
                                // worry about equality!
                                final TaxaConstraint ckm1 = ((TaxaConstraint)next.get(k - 1));
                                final TaxaConstraint ck = ((TaxaConstraint)next.get(k));
                                int intersectionSize = sizeOfIntersection(ckm1.taxons, ck.taxons);
                                if( intersectionSize != ckm1.taxons.getTaxonCount() ) {
                                    throw new XMLParseException(setsNotCOmpatibleMessage);
                                }
                                if( ckm1.upper > ck.upper ) {
                                   ckm1.upper = ck.upper;
                                }
                            }

                            // build tree for first subset
                            final TaxaConstraint taxaConstraint = ((TaxaConstraint)next.get(0));
                            SimpleTree tree = simulator.simulateTree(taxaConstraint.taxons, demoModel);
                            if( taxaConstraint.realLimits() ) {
                               attemptToScaleTree(tree, (taxaConstraint.lower + taxaConstraint.upper)/2);
                            }

                            // add more trees incrementally
                            for(int k = 1; k < next.size(); ++k) {
                                final TaxaConstraint constraintj = ((TaxaConstraint)next.get(k));
                                // build tree for taxons in difference
                                final Taxa list = new Taxa();
                                for(int j = 0; j < constraintj.taxons.getTaxonCount(); ++j) {
                                    Taxon taxonj = taxa.getTaxon(j);
                                    if( tree.getTaxonIndex(taxonj) < 0 ) {
                                        list.addTaxon(taxonj);
                                    }
                                }

                                MutableTree treeForRemaining = simulator.simulateTree(list, demoModel);
                                if( constraintj.realLimits() ) {
                                    double low = Math.max(constraintj.lower, tree.getNodeHeight(tree.getRoot()));
                                    attemptToScaleTree(treeForRemaining, 0.75 * low + 0.25 * constraintj.upper);

                                    // combine the trees
                                    final SimpleNode newRoot = new SimpleNode();
                                    final SimpleNode node = new SimpleNode(tree, tree.getRoot());
                                    newRoot.addChild(node);
                                    newRoot.addChild(new SimpleNode(treeForRemaining, treeForRemaining.getRoot()));
                                    newRoot.setHeight(0.5 * low + 0.5 * constraintj.upper);
                                    tree = new SimpleTree(newRoot);
                                } else {
                                    tree = simulator.simulateTree(new Tree[]{tree, treeForRemaining} , demoModel, -1);
                                }
                            }
                            st.add(tree);
                        }

                        // add a taxon list for remaining taxa
                        final Taxa list = new Taxa();
                        for(int j = 0; j < taxa.getTaxonCount(); ++j) {
                            Taxon taxonj = taxa.getTaxon(j);
                            for(int k = 0; k < st.size(); ++k) {
                                if( ((Tree)st.get(k)).getTaxonIndex(taxonj) >= 0 ) {
                                    taxonj = null;
                                    break;
                                }
                            }
                            if( taxonj != null ){
                                list.addTaxon(taxonj);
                            }
                        }
                        if( list.getTaxonCount() > 0 ) {
                           taxonLists.add(list);
                        }
                        if( st.size() > 1 ) {
                          final Tree tree1 = simulator.simulateTree((Tree[])st.toArray(new Tree[]{}), demoModel, -1);
                          subtrees.add(tree1);
                        } else {
                           subtrees.add(st.get(0));
                        }
                    }
                }
            }

            if (taxonLists.size() == 0 && subtrees.size() < 2) throw new XMLParseException("Expected at least one taxonList or two subtrees in " + getParserName() + " element.");

            Tree tree = null;

            try {
                Tree[] trees = new Tree[taxonLists.size()+subtrees.size()];
                // simulate each taxonList separately
                for (int i = 0; i < taxonLists.size(); i++) {
                    trees[i] = simulator.simulateTree((TaxonList)taxonLists.get(i), demoModel);
                }
                // add the preset trees
                for (int i = 0; i < subtrees.size(); i++) {
                    trees[i+taxonLists.size()] = (Tree)subtrees.get(i);
                }

                tree = simulator.simulateTree(trees, demoModel, rootHeight);
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException(iae.getMessage());
            }
            return tree;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a simulated tree under the given demographic model.";
        }

        public Class getReturnType() { return Object.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(RESCALE_HEIGHT, true, "Attempt to rescale the tree to the given root height"),
                new XORRule(
                        new ElementRule(ConstantPopulationModel.class),
                        new ElementRule(ExponentialGrowthModel.class))
        };
    };

    dr.evolution.coalescent.CoalescentSimulator simulator = new dr.evolution.coalescent.CoalescentSimulator();
}
