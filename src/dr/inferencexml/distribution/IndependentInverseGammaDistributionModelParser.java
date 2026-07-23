/*
 * IndependentInverseGammaDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class IndependentInverseGammaDistributionModelParser extends AbstractXMLObjectParser{
    public final static String INDEPENDENT_INVERSE_GAMMA_DISTRIBUTION_MODEL = "independentInverseGammaDistributionModel";
    public final static String ID = "id";
    public final static String SHAPE = "shape";
    public final static String SCALE = "scale";
    public final static String DATA = "data";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String id = xo.getStringAttribute(ID);
        Parameter shape = (Parameter) xo.getChild(SHAPE).getChild(Parameter.class);
        Parameter scale = (Parameter) xo.getChild(SCALE).getChild(Parameter.class);
        Parameter data = (Parameter) xo.getChild(DATA).getChild(Parameter.class);


        return new IndependentInverseGammaDistributionModel(id, shape, scale, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SHAPE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
               new ElementRule(SCALE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
               }),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            })
    };


    @Override
    public String getParserDescription() {
        return "Independent inverse gamma distributions with different parameters";
    }

    @Override
    public Class getReturnType() {
        return IndependentInverseGammaDistributionModel.class;
    }

    @Override
    public String getParserName() {
        return INDEPENDENT_INVERSE_GAMMA_DISTRIBUTION_MODEL;
    }
}
