/*
 * TransformedTreeTraitParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.TransformedTreeTraitProvider;
import dr.evolution.tree.TreeTraitProvider;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TransformedTreeTraitParser extends AbstractXMLObjectParser {

    public static final String NAME = "transformedTrait";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeTraitProvider treeTraits = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);
        Transform transform = (Transform) xo.getChild(Transform.class);

        return new TransformedTreeTraitProvider(treeTraits, transform);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TransformedTreeTraitProvider.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

    private XMLSyntaxRule[] rules = {
            new ElementRule(TreeTraitProvider.class),
            new ElementRule(Transform.class),
    };
}
