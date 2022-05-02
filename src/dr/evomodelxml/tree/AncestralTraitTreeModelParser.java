/*
 * TransformedTreeModelParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.evomodel.tree.*;
import dr.inference.model.FastMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.tree.TreeModelParser.*;

/**
 * @author Marc Suchard
 */
public class AncestralTraitTreeModelParser extends AbstractXMLObjectParser {

    private static final String ANCESTRAL_TRAIT_TREE_MODEL = "ancestralTraitTreeModel";
    private static final String ANCESTOR = "ancestor";
    private static final String ANCESTRAL_PATH = "ancestralPath";
    private static final String RELATIVE_HEIGHT = "relativeToTipHeight";
    private static final String ALWAYS_CONSTRAIN_TO_ROOT = "alwaysConstrainedToRoot";

    public String getParserName() {
        return ANCESTRAL_TRAIT_TREE_MODEL;
    }

    /**
     * @return a tree object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MutableTreeModel tree = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
        List<AncestralTaxonInTree> ancestors = parseAllAncestors(tree, xo);

        int index = tree.getExternalNodeCount();
        for (AncestralTaxonInTree ancestor : ancestors) {
            ancestor.setIndex(index);
            ancestor.setNode(new NodeRef() {
                @Override
                public int getNumber() {
                    return 0;
                }

                @Override
                public void setNumber(int n) {

                }
            });
        }

        if (xo.hasChildNamed(NODE_TRAITS)) {
            for (XMLObject cxo : xo.getAllChildren(NODE_TRAITS)) {
                parseNodeTraits(cxo, tree, ancestors);
            }
        }

        return new AncestralTraitTreeModel(xo.getId(), tree, ancestors);
    }

    private static void parseNodeTraits(XMLObject cxo, Tree tree, List<AncestralTaxonInTree> ancestors)
            throws XMLParseException {

        String name = cxo.getAttribute(NAME, "trait");
        int dim = cxo.getAttribute(MULTIVARIATE_TRAIT, 1);

//        double[] initialValues = null;
//        if (cxo.hasAttribute(INITIAL_VALUE)) {
//            initialValues = cxo.getDoubleArrayAttribute(INITIAL_VALUE);
//        }


        final int rowDim = dim;
        final int colDim = tree.getExternalNodeCount() + ancestors.size();

        FastMatrixParameter parameter = new FastMatrixParameter(name, rowDim, colDim, 0.0);
        parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                rowDim * colDim));

        int parameterIndex = 0;
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
//            if (leafNodes) {
//                nodes[i].addTraitParameter(name, parameter.getParameter(parameterIndex), initialValues, firesTreeEvents);
                ++parameterIndex;
//            }
        }

        for (AncestralTaxonInTree ancestor : ancestors) {
            ++parameterIndex;
        }


        ParameterParser.replaceParameter(cxo, parameter);
    }

    private static List<AncestralTaxonInTree> parseAllAncestors(MutableTreeModel tree, XMLObject xo) throws XMLParseException {
        int index = tree.getExternalNodeCount();
        List<AncestralTaxonInTree> ancestors = new ArrayList<AncestralTaxonInTree>();
        for (XMLObject cxo : xo.getAllChildren(ANCESTOR)) {
            ancestors.add(parseAncestor(tree, cxo, index));
            ++index;
        }
        return ancestors;
    }

    private static AncestralTaxonInTree parseAncestor(MutableTreeModel tree, XMLObject xo, final int index) throws XMLParseException {

        Taxon ancestor = (Taxon) xo.getChild(Taxon.class);
        Parameter pseudoBranchLength = (Parameter) xo.getChild(Parameter.class);

        AncestralTaxonInTree ancestorInTree;

        NodeRef node = new NodeRef() {
            @Override
            public int getNumber() {
                return index;
            }

            @Override
            public void setNumber(int n) {
                throw new RuntimeException("Do not set");
            }
        };

        if (xo.hasChildNamed(MonophylyStatisticParser.MRCA)) {
            TaxonList descendants = MonophylyStatisticParser.parseTaxonListOrTaxa(
                    xo.getChild(MonophylyStatisticParser.MRCA));

            boolean alwaysContrainedToRoot = xo.getAttribute(ALWAYS_CONSTRAIN_TO_ROOT, false);

            if (alwaysContrainedToRoot) {

                try {

                if (tree.getRoot() != TreeUtils.getCommonAncestorNode(tree,
                        TreeUtils.getLeavesForTaxa(tree, descendants))) {
                    throw new XMLParseException("MRCA is not the root.");

                }
                } catch (TreeUtils.MissingTaxonException e) {
                    throw new XMLParseException(e.getMessage());
                }
            }

            try {
                ancestorInTree = new AncestralTaxonInTree(ancestor, tree, descendants, pseudoBranchLength,
                        null, node, index, 0.0, alwaysContrainedToRoot);
            } catch (TreeUtils.MissingTaxonException e) {
                throw new XMLParseException("Unable to find taxa for " + ancestor.getId());
            }
        } else {

            XMLObject cxo = xo.getChild(ANCESTRAL_PATH);

            Taxon taxon = (Taxon) cxo.getChild(Taxon.class);
            Parameter time = (Parameter) cxo.getChild(Parameter.class);

            boolean relativeHeight = cxo.getAttribute(RELATIVE_HEIGHT, false);
            double offset = 0;

            if (relativeHeight) {

                Object date = taxon.getAttribute("date");
                if (date != null && date instanceof Date) {
                    offset = taxon.getHeight();
                } else {
                    throw new XMLParseException("Taxon '" + taxon.getId() + "' has no specified date");
                }

            } else {

                if (time.getParameterValue(0) <= taxon.getHeight()) {
                    throw new XMLParseException("Ancestral path time must be > sampling time for taxon '" +
                            taxon.getId() + "'");
                }
            }

            Taxa descendent = new Taxa();
            descendent.addTaxon(taxon);

            try {
                ancestorInTree = new AncestralTaxonInTree(ancestor, tree, descendent, pseudoBranchLength,
                        time, node, index, offset, false); // TODO Refactor into separate class from MRCA version
            } catch (TreeUtils.MissingTaxonException e) {
                throw new XMLParseException("Unable to find taxa for " + ancestor.getId());
            }
        }

        return ancestorInTree;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a transformed model of the tree.";
    }

    public Class getReturnType() {
        return TransformedTreeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules =
            new XMLSyntaxRule[]{
                    new ElementRule(MutableTreeModel.class),
                    new ElementRule(ANCESTOR, new XMLSyntaxRule[] {
                            new ElementRule(Taxon.class),
                            new ElementRule(Parameter.class),
                            new XORRule(
                                    new ElementRule(MonophylyStatisticParser.MRCA, new XMLSyntaxRule[]{
                                            new XORRule(
                                                    new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                                                    new ElementRule(Taxa.class)
                                            ),
                                            AttributeRule.newBooleanRule(ALWAYS_CONSTRAIN_TO_ROOT, true),
                                    }),
                                    new ElementRule(ANCESTRAL_PATH, new XMLSyntaxRule[]{
                                            new ElementRule(Taxon.class),
                                            new ElementRule(Parameter.class),
                                            AttributeRule.newBooleanRule(RELATIVE_HEIGHT, true),
                                    })),
                    }, 0, Integer.MAX_VALUE),
                    nodeTraitsRule,
            };
}
