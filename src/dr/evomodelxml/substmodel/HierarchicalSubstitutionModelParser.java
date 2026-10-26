/*
 * HierarchicalSubstitutionModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.AbstractHierarchicalRateProvider;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.RelativeRatesProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Hierarchical Rate Substitution Model
 *
 * @author Jeff Thorne
 * @author Tae-Kun Seo
 * @author Sangchul Choi
 * @author Xiang Ji
 */
public class HierarchicalSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String HIERARCHICAL_SUBSTITUTION_MODEL = "hierarchicalSubstitutionModel";
    public static final String CHARACTER_TREE = "characterTree";
    public static final String RATES_PARAMETER = "ratesParameter";
    public static final String NORMALISED = "normalised";
    public static final String KAPPA = "kappa";
    public static final String KAPPA1 = "kappa1";
    public static final String KAPPA2 = "kappa2";
    public static final String CHARACTER_TREE_NEWICK = "characterTreeNewick";
    public static final String RATES_FROM = "ratesFrom";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final SubstitutionModel borrowed = xo.hasChildNamed(RATES_FROM)
                ? (SubstitutionModel) xo.getElementFirstChild(RATES_FROM) : null;
        FrequencyModel freqModel = xo.hasChildNamed(FrequencyModelParser.FREQUENCIES)
                ? (FrequencyModel) xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES) : null;
        if (borrowed != null) {
            if (freqModel == null) {
                freqModel = borrowed.getFrequencyModel();
            } else if (freqModel != borrowed.getFrequencyModel()) {
                throw new XMLParseException("hierarchicalSubstitutionModel: the frequency model of ratesFrom's model must "
                        + "be the hierarchicalSubstitutionModel's frequency model (give it by idref, or leave frequencies out)");
            }
        } else if (freqModel == null) {
            throw new XMLParseException("hierarchicalSubstitutionModel: frequencies is required unless ratesFrom gives the model");
        }

        List<String> selectors = new ArrayList<String>();
        for (String name : new String[]{RATES_PARAMETER, KAPPA, KAPPA1, RATES_FROM}) {
            if (xo.hasChildNamed(name)) {
                selectors.add(name.equals(KAPPA1) ? "kappa1 and kappa2" : name);
            }
        }
        if (selectors.size() > 1) {
            throw new XMLParseException("hierarchicalSubstitutionModel: give at most one of ratesParameter, kappa, "
                    + "kappa1 and kappa2, ratesFrom; found " + selectors);
        }
        if (xo.hasChildNamed(KAPPA1) != xo.hasChildNamed(KAPPA2)) {
            throw new XMLParseException("hierarchicalSubstitutionModel: kappa1 and kappa2 must be given together");
        }
        final String selector = selectors.isEmpty() ? null : selectors.get(0);
        final boolean nucleotideForm = KAPPA.equals(selector) || "kappa1 and kappa2".equals(selector)
                || RATES_FROM.equals(selector);

        try {
            TreeModel tree = characterTree(xo, freqModel, selector, nucleotideForm);
            HierarchicalSubstitutionModel.HierarchicalRateProvider provider;
            String form;
            if (borrowed != null) {
                provider = RelativeRatesProvider.of(borrowed, tree);
                form = "ratesFrom " + (borrowed instanceof HKY ? "HKY (kappa)" : "TN93 (kappa1, kappa2)");
            } else if (KAPPA.equals(selector)) {
                provider = RelativeRatesProvider.hky((Parameter) xo.getElementFirstChild(KAPPA), freqModel, tree);
                form = "HKY (kappa)";
            } else if (nucleotideForm) {
                provider = RelativeRatesProvider.tn93((Parameter) xo.getElementFirstChild(KAPPA1),
                        (Parameter) xo.getElementFirstChild(KAPPA2), freqModel, tree);
                form = "TN93 (kappa1, kappa2)";
            } else if (RATES_PARAMETER.equals(selector)) {
                XMLObject cxo = xo.getChild(RATES_PARAMETER);
                boolean normalised = cxo.getAttribute(NORMALISED, false);
                provider = new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(
                        (Parameter) cxo.getChild(Parameter.class), tree, freqModel, normalised);
                form = normalised ? "ratesParameter, normalised" : "ratesParameter";
            } else {
                provider = new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(
                        new Parameter.Default(tree.getInternalNodeCount(), 1.0), tree);
                form = "no selector: Default with all rates 1.0";
            }
            HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel(HIERARCHICAL_SUBSTITUTION_MODEL,
                    provider, freqModel);
            Logger.getLogger("dr.evomodel").info("Hierarchical substitution model" + (xo.hasId() ? " '" + xo.getId() + "'" : "")
                    + ": " + form + "; character tree clades " + clades(tree, freqModel));
            return model;
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(!xo.hasId() ? e.getMessage()
                    : "hierarchicalSubstitutionModel '" + xo.getId() + "': " + e.getMessage());
        }
    }

    private static TreeModel characterTree(XMLObject xo, FrequencyModel freqModel, String selector,
                                           boolean nucleotideForm) throws XMLParseException {
        final boolean hasTree = xo.hasChildNamed(CHARACTER_TREE);
        final boolean hasNewick = xo.hasAttribute(CHARACTER_TREE_NEWICK);
        if (hasTree && hasNewick) {
            throw new XMLParseException("hierarchicalSubstitutionModel: give characterTree or characterTreeNewick, not both");
        }
        if (hasTree) {
            return (TreeModel) xo.getElementFirstChild(CHARACTER_TREE);
        }
        if (hasNewick) {
            return AbstractHierarchicalRateProvider.characterTree(xo.getStringAttribute(CHARACTER_TREE_NEWICK));
        }
        if (nucleotideForm) {
            if (!(freqModel.getDataType() instanceof Nucleotides)) {
                throw new XMLParseException("hierarchicalSubstitutionModel: " + selector + " needs a nucleotide frequency model; found "
                        + freqModel.getDataType().getDescription());
            }
            return AbstractHierarchicalRateProvider.characterTree(AbstractHierarchicalRateProvider.NUCLEOTIDE_NEWICK);
        }
        throw new XMLParseException("hierarchicalSubstitutionModel: characterTree is required with " + (selector == null ? "no rate selector" : selector)
                + " on data type " + freqModel.getDataType().getDescription());
    }

    private static String clades(TreeModel tree, FrequencyModel freqModel) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            List<String> codes = new ArrayList<String>();
            collectCodes(tree, tree.getInternalNode(i), freqModel, codes);
            java.util.Collections.sort(codes);
            if (sb.length() > 0) {
                sb.append(' ');
            }
            sb.append('{');
            for (int k = 0; k < codes.size(); k++) {
                sb.append(k > 0 ? "," : "").append(codes.get(k));
            }
            sb.append('}');
        }
        return sb.toString();
    }

    private static void collectCodes(TreeModel tree, dr.evolution.tree.NodeRef node, FrequencyModel freqModel,
                                     List<String> codes) {
        if (tree.isExternal(node)) {
            codes.add(tree.getTaxonId(node.getNumber()));
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                collectCodes(tree, tree.getChild(node, i), freqModel, codes);
            }
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(CHARACTER_TREE_NEWICK, true),
            new ElementRule(FrequencyModelParser.FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}, true),
            new ElementRule(RATES_FROM,
                    new XMLSyntaxRule[]{new ElementRule(SubstitutionModel.class)}, true),
            new ElementRule(CHARACTER_TREE,
                    new XMLSyntaxRule[]{new ElementRule(TreeModel.class)}, true),
            new ElementRule(RATES_PARAMETER,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(NORMALISED, true),
                            new ElementRule(Parameter.class)}, true),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(KAPPA1,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(KAPPA2,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };

    @Override
    public String getParserDescription() {
        return "";
    }

    @Override
    public Class getReturnType() {
        return HierarchicalSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return HIERARCHICAL_SUBSTITUTION_MODEL;
    }
}
