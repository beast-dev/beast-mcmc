/*
 * LatentFactorModelParser.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class LatentFactorModelParser extends AbstractXMLObjectParser {
    public final static String LATENT_FACTOR_MODEL = "latentFactorModel";
    public final static String NUMBER_OF_FACTORS = "factorNumber";
    public final static String FACTORS = "factors";
    public final static String DATA = "data";
    public final static String LATENT = "latent";

    public String getParserName() {
        return LATENT_FACTOR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//        Parameter data = null;
//        Parameter factor  = null;
//        Parameter latent  = null;

        Parameter factor = (Parameter) xo.getChild(FACTORS).getChild(Parameter.class);

//        int factors=xo.getAttribute(NUMBER_OF_FACTORS, 4);
//        XMLObject dataXML= xo.getChild(DATA);
//        if (dataXML.getChild(0) instanceof Parameter) {
//            data = (Parameter) dataXML.getChild(Parameter.class);
//        }
//        else {
//            data = new Parameter.Default(dataXML.getDoubleChild(0));
//        }
//exit if there's no parameter? What happens?

//        int factors=xo.getAttribute(NUMBER_OF_FACTORS, 4);
//        XMLObject dataXML= xo.getChild(DATA);
//        if (dataXML.getChild(0) instanceof Parameter) {
//            data = (Parameter) dataXML.getChild(Parameter.class);
//        }
        return new LatentFactorModel(null, factor, null);
    }

    private static final XMLSyntaxRule[] rules = {
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            })
    };

//    <latentFactorModel>
//      <factors>
//         <parameter idref="factors"/>
//      </factors>
//    </latentFactorModel>


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "TODO"; // TODO
    }

    @Override
    public Class getReturnType() {
        return LatentFactorModel.class;
    }


}
