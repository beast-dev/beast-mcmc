/*
 * ScaledTreeLengthRateModelParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ScaledTreeLengthRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class ScaledTreeLengthRateModelParser extends AbstractXMLObjectParser {
    public static final String MODEL_NAME = "scaledTreeLengthModel";
    public static final String SCALING_FACTOR = "scalingFactor";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter totalLength = (Parameter) xo.getElementFirstChild(SCALING_FACTOR);
        if (totalLength == null) {
            totalLength = new Parameter.Default(1, 1.0);
        }
        Logger.getLogger("dr.evomodel.branchratemodel").info("\n ---------------------------------\nCreating ScaledTreeLengthRateModel model.");
        Logger.getLogger("dr.evomodel.branchratemodel").info("\tTotal tree length will be scaled to " + totalLength.getParameterValue(0) + ".");
        Logger.getLogger("dr.evomodel.branchratemodel").info("\tIf you publish results using this rate model, please cite Alekseyenko, Lee and Suchard (2008) Syst. Biol 57: 772-784." +
                "\n---------------------------------\n");

        return new ScaledTreeLengthRateModel(tree, totalLength);
    }

    public String getParserDescription() {
        return "This element returns a branch rate model that scales the total length " +
                "of the tree to specified valued (default=1.0).";
    }

    public Class getReturnType() {
        return ScaledTreeLengthRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SCALING_FACTOR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(TreeModel.class)
    };

}
