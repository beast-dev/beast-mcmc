/*
 * AlloppHybPopSizesScaleParser.java
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

package dr.evomodel.alloppnet.parsers;

import dr.evomodel.alloppnet.operators.AlloppHybPopSizesScale;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.inference.operators.MCMCOperator;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 03/08/12
 */
public class AlloppHybPopSizesScaleParser extends AbstractXMLObjectParser {

    public static final String HYB_POP_SIZES_SCALE = "hybPopSizesScaleOperator";


    public String getParserName() {
        return HYB_POP_SIZES_SCALE;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
        AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scalingfactor = xo.getDoubleAttribute(ScaleOperatorParser.SCALE_FACTOR);
        return new AlloppHybPopSizesScale(apspnet, apsp, scalingfactor, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(ScaleOperatorParser.SCALE_FACTOR),
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator which scales the population size of a newly formed hybrid.";

    }

    @Override
    public Class getReturnType() {
        return AlloppHybPopSizesScale.class;
    }

}
