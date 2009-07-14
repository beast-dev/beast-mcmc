/*
 * CoalescentSimulator.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.tree.*;
import dr.evolution.util.*;
import dr.evomodelxml.TreeModelParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.math.UnivariateFunction;
import dr.xml.*;

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
    public static final String ROOT_HEIGHT = TreeModelParser.ROOT_HEIGHT;
    public static final String CONSTRAINED_TAXA = "constrainedTaxa";
    public static final String TMRCA_CONSTRAINT = "tmrca";
    public static final String IS_MONOPHYLETIC = "monophyletic";

    /**
     * Simulates a coalescent tree from a set of subtrees.
     */
    public CoalescentSimulator() {
    }

    /**
     * Simulates a coalescent tree from a set of subtrees.
     *
     * @param subtrees         an array of tree to be used as subtrees
     * @param model            the demographic model to use
     * @param rootHeight       an optional root height with which to scale the whole tree
     * @param preserveSubtrees true of subtrees should be preserved
     * @return a simulated coalescent tree
     */
    public SimpleTree simulateTree(Tree[] subtrees, DemographicModel model, double rootHeight, boolean preserveSubtrees) {

        SimpleNode[] roots = new SimpleNode[subtrees.length];
        SimpleTree tree;

        dr.evolution.util.Date mostRecent = null;
        for (Tree subtree : subtrees) {
            Date date = Tree.Utils.findMostRecentDate(subtree);
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (Tree subtree : subtrees) {
                Date date = Tree.Utils.findMostRecentDate(subtree);
                if (date != null) {
                    double diff = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    for (int j = 0; j < subtree.getNodeCount(); j++) {
                        NodeRef node = subtree.getNode(j);
                        ((SimpleTree) subtree).setNodeHeight(node, subtree.getNodeHeight(node) + diff);
                    }
                }
            }
        }

        for (int i = 0; i < roots.length; i++) {
            roots[i] = new SimpleNode(subtrees[i], subtrees[i].getRoot());
        }

        // if just one taxonList then finished
        if (roots.length == 1) {
            tree = new SimpleTree(roots[0]);
        } else {
            tree = new SimpleTree(simulator.simulateCoalescent(roots, model.getDemographicFunction()));
        }

        if (rootHeight > 0.0) {
            if (preserveSubtrees) {
                limitNodes(tree, rootHeight);
            } else {
                attemptToScaleTree(tree, rootHeight);
            }
        }

        return tree;
    }


    /**
     * Simulates a coalescent tree, given a taxon list.
     *
     * @param taxa  the set of taxa to simulate a coalescent tree between
     * @param model the demographic model to use
     * @return a simulated coalescent tree
     */
    public SimpleTree simulateTree(TaxonList taxa, DemographicModel model) {

        return simulator.simulateTree(taxa, model.getDemographicFunction());
    }

    /**
     * Clip nodes height above limit.
     *
     * @param tree  to clip
     * @param limit height limit
     */
    private static void limitNodes(MutableTree tree, double limit) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef n = tree.getInternalNode(i);
            if (tree.getNodeHeight(n) > limit) {
                tree.setNodeHeight(n, limit);
            }
        }
        MutableTree.Utils.correctHeightsForTips(tree);
    }

    private static void attemptToScaleTree(MutableTree tree, double rootHeight) {
        // avoid empty tree
        if (tree.getRoot() == null) return;

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
        final boolean isMonophyletic;
        double upper;

        TaxaConstraint(TaxonList taxons, ParametricDistributionModel p, boolean isMono) {
            this.taxons = taxons;
            this.isMonophyletic = isMono;

            if (p != null) {
                final UnivariateFunction univariateFunction = p.getProbabilityDensityFunction();
                lower = univariateFunction.getLowerBound();
                upper = univariateFunction.getUpperBound();
            } else {
                lower = 0;
                upper = Double.POSITIVE_INFINITY;
            }
        }

        TaxaConstraint(TaxonList taxons, double low, double high, boolean isMono) {
            this.taxons = taxons;
            this.isMonophyletic = isMono;
            upper = high;
            lower = low;
        }

        public boolean realLimits() {
            return lower != 0 || upper != Double.POSITIVE_INFINITY;
        }
    }

    static private int sizeOfIntersection(TaxonList tl1, TaxonList tl2) {
        int nIn = 0;
        for (int j = 0; j < tl1.getTaxonCount(); ++j) {
            if (tl2.getTaxonIndex(tl1.getTaxon(j)) >= 0) {
                ++nIn;
            }
        }
        return nIn;
    }


    static private boolean contained(TaxonList taxons, TaxonList taxons1) {
        return sizeOfIntersection(taxons, taxons1) == taxons.getTaxonCount();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COALESCENT_TREE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoalescentSimulator simulator = new CoalescentSimulator();

            DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);
            List<TaxonList> taxonLists = new ArrayList<TaxonList>();
            List<Tree> subtrees = new ArrayList<Tree>();

            double rootHeight = xo.getAttribute(ROOT_HEIGHT, -1.0);

            if (xo.hasAttribute(RESCALE_HEIGHT)) {
                rootHeight = xo.getDoubleAttribute(RESCALE_HEIGHT);
            }

            // should have one child that is node
            for (int i = 0; i < xo.getChildCount(); i++) {
                final Object child = xo.getChild(i);

                // AER - swapped the order of these round because Trees are TaxonLists...
                if (child instanceof Tree) {
                    subtrees.add((Tree) child);
                } else if (child instanceof TaxonList) {
                    taxonLists.add((TaxonList) child);
                } else if (xo.getChildName(i).equals(CONSTRAINED_TAXA)) {
                    rootHeight = -1; // ignore it? should we errror?

                    XMLObject constrainedTaxa = (XMLObject) child;

                    // all taxa
                    final TaxonList taxa = (TaxonList) constrainedTaxa.getChild(TaxonList.class);

                    List<TaxaConstraint> constraints = new ArrayList<TaxaConstraint>();
                    final String setsNotCompatibleMessage = "taxa sets not compatible";

                    // pick up all constraints. order in partial order, where taxa_1 @in taxa_2 implies
                    // taxa_1 is before taxa_2.


                    for (int nc = 0; nc < constrainedTaxa.getChildCount(); ++nc) {

                        final Object object = constrainedTaxa.getChild(nc);
                        if (object instanceof XMLObject) {
                            final XMLObject constraint = (XMLObject) object;

                            if (constraint.getName().equals(TMRCA_CONSTRAINT)) {
                                TaxonList taxaSubSet = (TaxonList) constraint.getChild(TaxonList.class);
                                ParametricDistributionModel dist =
                                        (ParametricDistributionModel) constraint.getChild(ParametricDistributionModel.class);
                                boolean isMono = constraint.getAttribute(IS_MONOPHYLETIC, true);

                                final TaxaConstraint taxaConstraint = new TaxaConstraint(taxaSubSet, dist, isMono);
                                int insertPoint;
                                for (insertPoint = 0; insertPoint < constraints.size(); ++insertPoint) {
                                    // if new <= constraints[insertPoint] insert before insertPoint

                                    final TaxaConstraint iConstraint = constraints.get(insertPoint);
                                    if (iConstraint.isMonophyletic) {
                                        if (!taxaConstraint.isMonophyletic) {
                                            continue;
                                        }

                                        final TaxonList taxonsip = iConstraint.taxons;
                                        final int nIn = sizeOfIntersection(taxonsip, taxaSubSet);
                                        if (nIn == taxaSubSet.getTaxonCount()) {
                                            break;
                                        }
                                        if (nIn > 0 && nIn != taxonsip.getTaxonCount()) {
                                            throw new XMLParseException(setsNotCompatibleMessage);
                                        }
                                    } else {
                                        // reached non mono area
                                        if (!taxaConstraint.isMonophyletic) {
                                            if (iConstraint.upper >= taxaConstraint.upper) {
                                                break;
                                            }
                                        } else {
                                            break;
                                        }
                                    }
                                }
                                constraints.add(insertPoint, taxaConstraint);
                            }
                        }
                    }
                    final int nConstraints = constraints.size();

                    if (nConstraints == 0) {
                        taxonLists.add(taxa);
                    } else {
                        for (int nc = 0; nc < nConstraints; ++nc) {
                            TaxaConstraint cnc = constraints.get(nc);
                            if (!cnc.isMonophyletic) {
                                for (int nc1 = nc - 1; nc1 >= 0; --nc1) {
                                    TaxaConstraint cnc1 = constraints.get(nc1);
                                    int x = sizeOfIntersection(cnc.taxons, cnc1.taxons);
                                    if (x > 0) {
                                        Taxa combinedTaxa = new Taxa(cnc.taxons);
                                        combinedTaxa.addTaxa(cnc1.taxons);
                                        cnc = new TaxaConstraint(combinedTaxa, cnc.lower, cnc.upper, cnc.isMonophyletic);
                                        constraints.set(nc, cnc);
                                    }
                                }
                            }
                        }
                        // determine upper bound for each set.
                        double[] upper = new double[nConstraints];
                        for (int nc = nConstraints - 1; nc >= 0; --nc) {
                            final TaxaConstraint cnc = constraints.get(nc);
                            if (cnc.realLimits()) {
                                upper[nc] = cnc.upper;
                            } else {
                                upper[nc] = Double.POSITIVE_INFINITY;
                            }
                        }

                        for (int nc = nConstraints - 1; nc >= 0; --nc) {
                            final TaxaConstraint cnc = constraints.get(nc);
                            if (upper[nc] < Double.POSITIVE_INFINITY) {
                                for (int nc1 = nc - 1; nc1 >= 0; --nc1) {
                                    final TaxaConstraint cnc1 = constraints.get(nc1);
                                    if (contained(cnc1.taxons, cnc.taxons)) {
                                        upper[nc1] = Math.min(upper[nc1], upper[nc]);
                                        if (cnc1.realLimits() && cnc1.lower > upper[nc1]) {
                                            throw new XMLParseException(setsNotCompatibleMessage);
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        // collect subtrees here
                        List<Tree> st = new ArrayList<Tree>();
                        for (int nc = 0; nc < constraints.size(); ++nc) {
                            final TaxaConstraint nxt = constraints.get(nc);
                            // collect all previously built subtrees which are a subset of taxa set to be added
                            List<Tree> subs = new ArrayList<Tree>();
                            Taxa newTaxons = new Taxa(nxt.taxons);
                            for (int k = 0; k < st.size(); ++k) {
                                final Tree stk = st.get(k);
                                int x = sizeOfIntersection(stk, nxt.taxons);
                                if (x == st.get(k).getTaxonCount()) {
                                    final Tree tree = st.remove(k);
                                    --k;
                                    subs.add(tree);
                                    newTaxons.removeTaxa(tree);
                                }
                            }

                            SimpleTree tree = simulator.simulateTree(newTaxons, demoModel);
                            final double lower = nxt.realLimits() ? nxt.lower : 0;
                            if (upper[nc] < Double.MAX_VALUE) {
                                attemptToScaleTree(tree, (lower + upper[nc]) / 2);
                            }
                            if (subs.size() > 0) {
                                if (tree.getTaxonCount() > 0) subs.add(tree);
                                double h = -1;
                                if (upper[nc] < Double.MAX_VALUE) {
                                    for (Tree t : subs) {
                                        h = Math.max(h, t.getNodeHeight(t.getRoot()));
                                    }
                                    h = (h + upper[nc]) / 2;
                                }
                                tree = simulator.simulateTree(subs.toArray(new Tree[subs.size()]), demoModel, h, true);
                            }
                            st.add(tree);

                        }

                        // add a taxon list for remaining taxa
                        final Taxa list = new Taxa();
                        for (int j = 0; j < taxa.getTaxonCount(); ++j) {
                            Taxon taxonj = taxa.getTaxon(j);
                            for (Tree aSt : st) {
                                if (aSt.getTaxonIndex(taxonj) >= 0) {
                                    taxonj = null;
                                    break;
                                }
                            }
                            if (taxonj != null) {
                                list.addTaxon(taxonj);
                            }
                        }
                        if (list.getTaxonCount() > 0) {
                            taxonLists.add(list);
                        }
                        if (st.size() > 1) {
                            final Tree t = simulator.simulateTree(st.toArray(new Tree[st.size()]), demoModel, -1, false);
                            subtrees.add(t);
                        } else {
                            subtrees.add(st.get(0));
                        }
                    }
                }
            }

            if (taxonLists.size() == 0) {
                if (subtrees.size() == 1) {
                    return subtrees.get(0);
                }
                throw new XMLParseException("Expected at least one taxonList or two subtrees in "
                        + getParserName() + " element.");
            }

            try {
                Tree[] trees = new Tree[taxonLists.size() + subtrees.size()];
                // simulate each taxonList separately
                for (int i = 0; i < taxonLists.size(); i++) {
                    trees[i] = simulator.simulateTree(taxonLists.get(i), demoModel);
                }
                // add the preset trees
                for (int i = 0; i < subtrees.size(); i++) {
                    trees[i + taxonLists.size()] = subtrees.get(i);
                }
                
                return simulator.simulateTree(trees, demoModel, rootHeight, trees.length != 1);
                
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException(iae.getMessage());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a simulated tree under the given demographic model.";
        }

        public Class getReturnType() {
            return Object.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(RESCALE_HEIGHT, true, "Attempt to rescale the tree to the given root height"),
                AttributeRule.newDoubleRule(ROOT_HEIGHT, true, ""),
                new ElementRule(Tree.class, 0, Integer.MAX_VALUE),
                new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
                new ElementRule(CONSTRAINED_TAXA, new XMLSyntaxRule[]{
                        new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
                        // need more here
                }, true),
                new ElementRule(DemographicModel.class),
        };
    };

    dr.evolution.coalescent.CoalescentSimulator simulator = new dr.evolution.coalescent.CoalescentSimulator();
}
