/*
 * InstantaneousMixtureSubstitutionModelParser.java
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

package dr.evomodelxml.substmodel;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.AminoAcidMixture;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.SubstitutionRateMatrixMixture;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LogLinearModel;
import dr.inference.model.Parameter;
import dr.evomodel.substmodel.FrequencyModel;
import dr.oldevomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.xml.*;
import dr.evomodel.substmodel.InstantaneousMixtureSubstitutionModel;
import java.util.ArrayList;
import java.util.List;

public class InstantaneousMixtureSubstitutionModelParser extends AbstractXMLObjectParser {

    private static final String MIXTURE_MODEL = "instantaneousMixtureSubstitutionModel";
    private static final String MIXTURE_WEIGHTS = "mixtureWeights";
    private static final String RELATIVE = "weightsRelativeToOne";

    public String getParserName() {
        return MIXTURE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean transform = xo.getAttribute(RELATIVE, false);

        List<SubstitutionModel> modelList = xo.getAllChildren(SubstitutionModel.class);
        if ( modelList.size() > 2 ) {
            throw new RuntimeException("instantaneousMixtureSubstitutionModel only implemented for 2 components.");
        }

        DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

        XMLObject cxo = xo.getChild(ComplexSubstitutionModelParser.ROOT_FREQUENCIES);
        FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        if (dataType != rootFreq.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
        }

        Parameter weights = (Parameter) xo.getElementFirstChild(MIXTURE_WEIGHTS);
        if ( weights.getDimension() != 1 ) {
            throw new RuntimeException("Invalid number of weights for instantaneousMixtureSubstitutionModel which uses p, 1-p mixture.");
        }

        if ( weights.getDimension() != 1 || weights.getBounds().getLowerLimit(0) > Double.MIN_VALUE || Math.abs(1.0 - weights.getBounds().getUpperLimit(0)) > Double.MIN_VALUE ) {
            if (!transform) {
                throw new RuntimeException("Invalid bounds of weight parameter for instantaneousMixtureSubstitutionModel which uses p, 1-p mixture.");
            }
        }

        return new InstantaneousMixtureSubstitutionModel(xo.getId(),dataType,rootFreq,modelList,weights,transform);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Mixture of substitution models by mixing Q-matrices.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(RELATIVE, true),
            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(ComplexSubstitutionModelParser.ROOT_FREQUENCIES, FrequencyModel.class),
            new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE),
            new ElementRule(MIXTURE_WEIGHTS, Parameter.class)
    };
}