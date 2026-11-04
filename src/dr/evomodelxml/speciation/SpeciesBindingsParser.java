/*
 * SpeciesBindingsParser.java
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

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class SpeciesBindingsParser extends AbstractXMLObjectParser {
    public static final String SPECIES = "species";
    public static final String GENE_TREES = "geneTrees";
    public static final String GTREE = "gtree";
    public static final String PLOIDY = "ploidy";


    public String getParserName() {
        return SPECIES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<SpeciesBindings.SPinfo> sp = new ArrayList<SpeciesBindings.SPinfo>();
        for (int k = 0; k < xo.getChildCount(); ++k) {
            final Object child = xo.getChild(k);
            if (child instanceof SpeciesBindings.SPinfo) {
                sp.add((SpeciesBindings.SPinfo) child);
            }
        }

        final XMLObject xogt = xo.getChild(GENE_TREES);
        final int nTrees = xogt.getChildCount();
        final TreeModel[] trees = new TreeModel[nTrees];
        double[] popFactors = new double[nTrees];

        for (int nt = 0; nt < trees.length; ++nt) {
            Object child = xogt.getChild(nt);
            if (!(child instanceof TreeModel)) {
                assert child instanceof XMLObject;
                popFactors[nt] = ((XMLObject) child).getDoubleAttribute(PLOIDY);
                child = ((XMLObject) child).getChild(TreeModel.class);

            } else {
                popFactors[nt] = -1;
            }
            trees[nt] = (TreeModel) child;
        }

        try {
            return new SpeciesBindings(sp.toArray(new SpeciesBindings.SPinfo[sp.size()]), trees, popFactors);
        } catch (Error e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    /* Can't be tree because XML parser supports usage of global tags only as main tags */
    ElementRule treeWithPloidy = new ElementRule(GTREE,
            new XMLSyntaxRule[]{AttributeRule.newDoubleRule(PLOIDY),
                    new ElementRule(TreeModel.class)}, 0, Integer.MAX_VALUE);
    //XMLSyntaxRule[] someTree = {new OrRule(new ElementRule(TreeModel.class), treeWithPloidy)};

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciesBindings.SPinfo.class, 2, Integer.MAX_VALUE),
                // new ElementRule(GENE_TREES, someTree,  1, Integer.MAX_VALUE )
                new ElementRule(GENE_TREES,
                        new XMLSyntaxRule[]{
                                // start at 0 for only ploidy tree cases
                                new ElementRule(TreeModel.class, 0, Integer.MAX_VALUE),
                                treeWithPloidy
                        }),
        };
    }

    public String getParserDescription() {
        return "Binds taxa in gene trees with species information.";
    }

    public Class getReturnType() {
        return SpeciesBindings.class;
    }

}
