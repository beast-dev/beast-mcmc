/*
 * PopsIOTreeNodeSlideParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.PopsIOTreeNodeSlide;
import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * User: Graham Jones
 * Date: 11/05/12
 */
public class PopsIOTreeNodeSlideParser  extends AbstractXMLObjectParser {
    public static final String PIOTREE_NODESLIDE = "pioTreeNodeSlide";



    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        PopsIOSpeciesBindings piosb = (PopsIOSpeciesBindings) xo.getChild(PopsIOSpeciesBindings.class);
        PopsIOSpeciesTreeModel piostm = (PopsIOSpeciesTreeModel) xo.getChild(PopsIOSpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new PopsIOTreeNodeSlide(piostm, piosb, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(PopsIOSpeciesBindings.class),
                new ElementRule(PopsIOSpeciesTreeModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator for species tree: changes tree without breaking embedding of gene trees.";
    }

    @Override
    public Class getReturnType() {
        return PopsIOTreeNodeSlide.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getParserName() {
        return PIOTREE_NODESLIDE;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
