/*
 * YuleModel.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evoxml.XMLUnits;
import dr.inference.distribution.ParametricDistributionModel;
import dr.xml.*;
//import dr.evomodel.tree.ARGModel;


/**
 * This class contains methods that describe a uniform over unrooted trees with an iid branch length prior
 *
 * @author Marc Suchard
 */

public class MrBayesDefaultModel extends SpeciationModel {

    public static final String MRBAYES_DEFAULT_MODEL = "mrbayesDefaultModel";
    public static final String OUTGROUP = "artificialOutgroup";
    public static final String BRANCH_MODEL = "branchLengthModel";


    public MrBayesDefaultModel(Taxon outgroup, ParametricDistributionModel branchModel, Type units) {

        super(MRBAYES_DEFAULT_MODEL, units);

        this.branchModel = branchModel;
        this.outgroup = outgroup;

        addModel(branchModel);
    }

    //
    // functions that define a speciation model
    //
    public double logTreeProbability(int taxonCount) {
        throw new RuntimeException("Why was 'logTreeProbability' called?");
    }

    //
    // functions that define a speciation model
    //
    public double logNodeProbability(Tree tree, NodeRef node) {
        if (tree.getRoot() == node) {
            // Make sure that outgroup really is an outgroup
            if ((tree.getNodeTaxon(tree.getChild(node, 0)) == outgroup) ||
                    (tree.getNodeTaxon(tree.getChild(node, 1)) == outgroup))
                return 0;
            else
                return Double.NEGATIVE_INFINITY;
        }
        NodeRef parentNode = tree.getParent(node);
        double branchLength = tree.getNodeHeight(parentNode) - tree.getNodeHeight(node);
        if (tree.getNodeTaxon(node) == outgroup) {
            NodeRef sisterNode = tree.getChild(parentNode, 0);
            if (sisterNode == node)
                sisterNode = tree.getChild(parentNode, 1);
            branchLength += tree.getNodeHeight(parentNode) - tree.getNodeHeight(sisterNode);
        }

        return branchModel.logPdf(branchLength);
    }

    public boolean includeExternalNodesInLikelihoodCalculation() { return false; }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    /**
     * Parses an element from an DOM document into a SpeciationModel. Recognises
     * YuleModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MRBAYES_DEFAULT_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLParser.Utils.getUnitsAttr(xo);

            XMLObject cxo = (XMLObject) xo.getChild(OUTGROUP);
            Taxon outgroup = (Taxon) cxo.getChild(Taxon.class);

            cxo = (XMLObject) xo.getChild(BRANCH_MODEL);
            ParametricDistributionModel branchModel = (ParametricDistributionModel) cxo.getChild(ParametricDistributionModel.class);

            return new MrBayesDefaultModel(outgroup, branchModel, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model that puts uniform prior on all possible unrooted topologies.";
        }

        public Class getReturnType() {
            return SpeciationModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(OUTGROUP,
                        new XMLSyntaxRule[]{new ElementRule(Taxon.class)}),
                new ElementRule(BRANCH_MODEL,
                        new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)}),
                XMLUnits.SYNTAX_RULES[0]
        };
    };


    //Protected stuff
    private Taxon outgroup;
    private ParametricDistributionModel branchModel;
}
