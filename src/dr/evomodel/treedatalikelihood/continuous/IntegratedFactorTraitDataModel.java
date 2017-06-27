/*
 * IntegratedFactorTraitDataModel.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class IntegratedFactorTraitDataModel extends ContinuousTraitDataModel {

    public IntegratedFactorTraitDataModel(String name,
                                          CompoundParameter parameter,
                                          List<Integer> missingIndices,
                                          int dimTrait,
                                          PrecisionType precisionType) {

        super(name, parameter, missingIndices, dimTrait, precisionType);

        if (precisionType != PrecisionType.FULL) {
            throw new IllegalArgumentException("Integrated factor model requires full precision likelihood peeling");
        }
    }

    // TODO Move into separate class file
    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter loadings = (Parameter) xo.getElementFirstChild(LOADINGS);
            Parameter precision = (Parameter) xo.getElementFirstChild(PRECISION);

            return null;
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
            return IntegratedFactorTraitDataModel.class;
        }

        @Override
        public String getParserName() {
            return INTEGRATED_FACTOR_Model;
        }
    };

    public static final String INTEGRATED_FACTOR_Model = "integratedFactorModel";
    public static final String LOADINGS = "loadings";
    public static final String PRECISION = "precision";

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(LOADINGS, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(PRECISION, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),


    };
}
