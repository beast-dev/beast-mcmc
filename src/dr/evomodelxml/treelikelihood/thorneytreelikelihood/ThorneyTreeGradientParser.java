/*
 * ThorneyTreeGradientParser.java
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

package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.treelikelihood.thorneytreelikelihood.ThorneyTreeGradient;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ThorneyTreeLikelihood;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author JT McCrone
 * @version $Id$
 */
public class ThorneyTreeGradientParser extends AbstractXMLObjectParser {

    private final String THORNEY_TREE_GRADIENT = "thorneyTreeGradient";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ThorneyTreeLikelihood likelihood = (ThorneyTreeLikelihood) xo.getChild(ThorneyTreeLikelihood.class);

        return new ThorneyTreeGradient(likelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ThorneyTreeLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Construct ThorneyTreeGradient";
    }

    @Override
    public Class getReturnType() {
        return ThorneyTreeGradient.class;
    }

    @Override
    public String getParserName() {
        return THORNEY_TREE_GRADIENT;
    }
}
