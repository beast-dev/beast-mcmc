/*
 * CoalescentSimulatorParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * CoalescentSimulator parser. Replaces the one now called 'OldCoalescentSimulator' which
 * used the parser name 'coalescentTree'. Simulates a tree using arbitrarily nested subtrees
 * which can be scaled to have particular node heights to be compatible with calibrations.
 */
public class CoalescentSimulatorParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_SIMULATOR = "coalescentSimulator";
    public static final String HEIGHT = "height";
    public static final String CONSTRAINTS_TREE = "constraintsTree";

    public String getParserName() {
        return COALESCENT_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoalescentSimulator simulator = new CoalescentSimulator();

        DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);
        List<TaxonList> taxonLists = new ArrayList<TaxonList>();
        List<Tree> subtrees = new ArrayList<Tree>();

        double height = xo.getAttribute(HEIGHT, Double.NaN);

        Tree constraintsTree =null;
        if(xo.hasChildNamed((CONSTRAINTS_TREE))){
            XMLObject cxo = xo.getChild(CONSTRAINTS_TREE);
            constraintsTree = (Tree) cxo.getChild(Tree.class);
        }

        // should have one child that is node
        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);
            if(child != constraintsTree){
                // AER - swapped the order of these round because Trees are TaxonLists...
                if (child instanceof Tree) {
                    subtrees.add((Tree) child);
                } else if (child instanceof TaxonList) {
                    taxonLists.add((TaxonList) child);
                }
            }
          }

        if (taxonLists.size() == 0) {
            if (subtrees.size() == 1) {
                return subtrees.get(0);
            } if (constraintsTree==null){
                throw new XMLParseException("Expected at least one taxonList or two subtrees or a constraints tree in "
                        + getParserName() + " element.");
            }
        }

        Taxa remainingTaxa = new Taxa();
        for (int i = 0; i < taxonLists.size(); i++) {
            remainingTaxa.addTaxa(taxonLists.get(i));
        }

        for (int i = 0; i < subtrees.size(); i++) {
            remainingTaxa.removeTaxa(subtrees.get(i));
        }
        if(constraintsTree!=null){
            remainingTaxa.removeTaxa(constraintsTree);
        }

        try {
            if(constraintsTree!=null){
                subtrees.add(simulator.simulateTree(constraintsTree,constraintsTree.getRoot(),demoModel));
            }

            Tree[] trees = new Tree[subtrees.size() + remainingTaxa.getTaxonCount()];
            // add the preset trees
            for (int i = 0; i < subtrees.size(); i++) {
                trees[i] = subtrees.get(i);
            }

            // add all the remaining taxa in as single tip trees...
            for (int i = 0; i < remainingTaxa.getTaxonCount(); i++) {
                Taxa tip = new Taxa();
                tip.addTaxon(remainingTaxa.getTaxon(i));
                trees[i + subtrees.size()] = simulator.simulateTree(tip, demoModel);
            }

            return simulator.simulateTree(trees, demoModel, height, trees.length != 1);

        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a simulated tree under the given demographic model. The element can " +
                "be nested to simulate with monophyletic clades. It also accepts an optionally resolved constraints" +
                "tree that that defines the topology of the taxa it includes. The tree will be rescaled to the"+
                "given height.";
    }

    public Class getReturnType() {
        return CoalescentSimulator.class; //Object.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(HEIGHT, true, ""),
            new ElementRule(Tree.class, 0, Integer.MAX_VALUE),
            new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
            new ElementRule(DemographicModel.class, 0, Integer.MAX_VALUE),
            new ElementRule(CONSTRAINTS_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class, 1, 1),
                    },true)
    };
}
