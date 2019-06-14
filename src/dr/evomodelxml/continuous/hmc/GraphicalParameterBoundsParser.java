/*
 * NodeHeightTransformParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightBounds;
import dr.inference.model.GraphicalParameterBound;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class GraphicalParameterBoundsParser extends AbstractXMLObjectParser {

    public static final String NAME = "graphicalParameterBounds";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        GraphicalParameterBound graphicalParameterBound;
        if (tree != null) {
            graphicalParameterBound = new NodeHeightBounds(parameter, tree);
        } else {
            graphicalParameterBound = new GraphicalParameterBound.FixedBound(parameter);
        }
        return graphicalParameterBound;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class, true),
                new ElementRule(Parameter.class)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GraphicalParameterBound.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
