/*
 * DuplicationTreeModelParser.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.evomodel.tree.DuplicationTreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;


/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class DuplicationTreeModelParser extends AbstractXMLObjectParser {

    private static final String DUPLICATION_TREE_MODEL = "duplicationTreeModel";
    private static final String DUPLICATION = "duplication";
    private static final String DUPLICATION_TIME_RATIO = "duplicationTimeRatio";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MutableTreeModel tree = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
        List<AncestralTaxonInTree> duplications = AncestralTraitTreeModelParser.parseAllAncestors(tree, xo, DUPLICATION);
        List<Parameter> duplicationTimeRatios = parseAllDuplicationTimes(xo);

        if (duplications.size() != duplicationTimeRatios.size()) {
            throw new RuntimeException("Each duplication event needs a ratio for its timing");
        }

        int index = tree.getExternalNodeCount();
        for (AncestralTaxonInTree duplication : duplications) {
            duplication.setIndex(index);
            duplication.setNode(new NodeRef() {
                @Override
                public int getNumber() {
                    return 0;
                }

                @Override
                public void setNumber(int n) {

                }
            });
        }

        return new DuplicationTreeModel(xo.getId(), tree, duplications, duplicationTimeRatios);
    }

    private List<Parameter> parseAllDuplicationTimes(XMLObject xo) {
        List<Parameter> duplicationTimes = new ArrayList<Parameter>();
        for (XMLObject cxo : xo.getAllChildren(DUPLICATION_TIME_RATIO)) {
            duplicationTimes.add((Parameter) cxo.getChild(Parameter.class));
        }
        return duplicationTimes;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules =
        new XMLSyntaxRule[]{
            new ElementRule(MutableTreeModel.class),
            new ElementRule(DUPLICATION, new XMLSyntaxRule[] {
                new ElementRule(Taxon.class),
                new ElementRule(Parameter.class, "Branch length towards duplication node."),
                new ElementRule(MonophylyStatisticParser.MRCA, Taxa.class),
            }, 1, 1),
            new ElementRule(DUPLICATION_TIME_RATIO, Parameter.class, "A ratio for duplication timing on the original branch", 1, 1),
        };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DuplicationTreeModel.class;
    }

    @Override
    public String getParserName() {
        return DUPLICATION_TREE_MODEL;
    }
}
