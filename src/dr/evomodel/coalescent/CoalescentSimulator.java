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
import dr.evolution.util.TaxonList;
import dr.evolution.util.TimeScale;
import dr.xml.*;

import java.util.ArrayList;


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
    public Tree simulateTree(Tree[] subtrees, DemographicModel model, double rootHeight) {

        SimpleNode[] roots = new SimpleNode[subtrees.length];
        MutableTree tree = null;

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
    public Tree simulateTree(TaxonList taxa, DemographicModel model) {

        return simulator.simulateTree(taxa, model.getDemographicFunction());
    }

    private void attemptToScaleTree(MutableTree tree, double rootHeight) {
        double scale = rootHeight / tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef n = tree.getInternalNode(i);
            tree.setNodeHeight(n, tree.getNodeHeight(n) * scale);
        }
        MutableTree.Utils.correctHeightsForTips(tree);
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
            for (int i =0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);

                // AER - swapped the order of these round because Trees are TaxonLists...
                if (child instanceof Tree) {
                    subtrees.add(child);
                } else if (child instanceof TaxonList) {
                    taxonLists.add(child);
                }
            }

            if (taxonLists.size() == 0 && subtrees.size() < 2) throw new XMLParseException("Expected at least one taxonList or two subtrees in " + getParserName() + " element.");

            Tree tree = null;

            try {
                Tree[] trees = new Tree[taxonLists.size()+subtrees.size()];
                // simulate each taxonList separately
                for (int i =0; i < taxonLists.size(); i++) {
                    trees[i] = simulator.simulateTree((TaxonList)taxonLists.get(i), demoModel);
                }
                // add the preset trees
                for (int i =0; i < subtrees.size(); i++) {
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