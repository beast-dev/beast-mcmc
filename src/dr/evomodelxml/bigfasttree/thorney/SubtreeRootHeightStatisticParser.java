/*
 * SubtreeRootHeightStatisticParser.java
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

package dr.evomodelxml.bigfasttree.thorney;


import dr.evolution.util.TaxonList;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.SubtreeRootHeightStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

public class SubtreeRootHeightStatisticParser extends AbstractXMLObjectParser {
    public final static String SUBTREE_ROOT_HEIGHT_STATISTIC = "subtreeRootHeightStatistic";
    public static final String ABSOLUTE = "absolute";
    public static final String MRCA = "mrca";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ConstrainedTreeModel tree = (ConstrainedTreeModel) xo.getChild(ConstrainedTreeModel.class);

        TaxonList taxa = null;

        if (xo.hasChildNamed(MRCA)) {
            taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        }
        boolean isAbsolute = xo.getAttribute(ABSOLUTE, false);
        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        return new SubtreeRootHeightStatistic(name,tree,taxa,isAbsolute);
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ConstrainedTreeModel.class),
        };    }

    @Override
    public String getParserDescription() {
        return "A parser for logging the heights of subtrees in a constrained tree model";
    }

    @Override
    public Class getReturnType() {
        return SubtreeRootHeightStatistic.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SUBTREE_ROOT_HEIGHT_STATISTIC;
    }
}