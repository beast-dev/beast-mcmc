/*
 * RescaledTreeParser.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evoxml;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.xml.*;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Parser takes a tree and rescales the node heights to match a set of clade heights defined by taxon sets.
 * @author Andrew Rambaut
 *
 * @version $Id$
 */
public class RescaledTreeParser extends AbstractXMLObjectParser {

    //
    // Public stuff
    //
    public final static String RESCALED_TREE = "rescaleTree";
    public final static String CLADE = "clade";
    public final static String HEIGHT = "height";

    public String getParserName() { return RESCALED_TREE; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree)xo.getChild(Tree.class);
        SimpleTree rescaledTree = new SimpleTree(tree);

        // First flag all internal nodes as unset....
        for (int i = 0; i < rescaledTree.getInternalNodeCount(); i++) {
            rescaledTree.setNodeHeight(rescaledTree.getInternalNode(i), Double.NEGATIVE_INFINITY);
        }

        for (int i = 0; i < xo.getChildCount(); i++) {
            XMLObject cxo = (XMLObject)xo.getChild(i);
            if (cxo.getName().equals(CLADE)) {
                TaxonList taxa = (TaxonList)cxo.getChild(TaxonList.class);
                double height = cxo.getDoubleAttribute(HEIGHT);

                Map<String, SimpleNode> leafNodes = new HashMap<String, SimpleNode>();

                Set<String> leafSet = new HashSet<String>();
                for (Taxon taxon : taxa) {
                    leafSet.add(taxon.getId());
                }

                NodeRef mrca = Tree.Utils.getCommonAncestorNode(tree, leafSet);
                if (mrca == null ||  Tree.Utils.getLeafCount(tree, mrca) != leafSet.size()) {
                    throw new XMLParseException("Clade defined by taxon Set, " + taxa.getId() + ", is not found in the guide tree");
                }

                rescaledTree.setNodeHeight(mrca, height);
            }
        }

        if (xo.hasAttribute(HEIGHT)) {
            rescaledTree.setNodeHeight(rescaledTree.getRoot(), xo.getDoubleAttribute(HEIGHT));
        }

        interpolateHeights(rescaledTree, rescaledTree.getRoot());

        return rescaledTree;
    }

    private double interpolateHeights(MutableTree tree, NodeRef node) {
        if (!tree.isExternal(node)) {

            double maxHeight = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef child = tree.getInternalNode(i);

                double h = interpolateHeights(tree, child);

                if (h > maxHeight) {
                    maxHeight = h;
                }
            }

            double height = tree.getNodeHeight(node);
            if (Double.isInfinite(height)) {
                int count = 1;
                NodeRef parent = tree.getParent(node);
                while (parent != null && Double.isInfinite(tree.getNodeHeight(parent))) {
                    parent = tree.getParent(parent);
                    count ++;
                }

                if (parent == null) {
                    height = maxHeight + count;
                } else {
                    height = (tree.getNodeHeight(parent) + maxHeight) / 2.0;
                }
            }

            tree.setNodeHeight(node, height);

            return height;
        } else {
            return tree.getNodeHeight(node);
        }
    }

    public String getParserDescription() {
        return "This element rescales a given tree with a set of clade heights.";
    }

    public Class getReturnType() { return Tree.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(HEIGHT),
            new ElementRule(Tree.class),
            new ElementRule(CLADE, new XMLSyntaxRule[] {
                    AttributeRule.newDoubleRule(HEIGHT),
                    new ElementRule(TaxonList.class)
            }, 1, Integer.MAX_VALUE)
    };
}