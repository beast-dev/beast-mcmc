/*
 * MarkovModulatedGY94CodonModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.ColtEigenSystem;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.codon.MarkovModulatedGY94CodonModel;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.HiddenCodons;
import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenSystem;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedGY94CodonModelParser extends GY94CodonModelParser {

    public static final String MARKOV_MODULATED_YANG_MODEL = "markovModulatedYangCodonModel";
    public static final String HIDDEN_COUNT = "hiddenCount";
    public static final String SWITCHING_RATES = "switchingRates";
    public static final String DIAGONALIZATION = "diagonalization";

    public String getParserName() {
        return MARKOV_MODULATED_YANG_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);
        HiddenCodons codons;
        if (dataType instanceof HiddenCodons)
            codons = (HiddenCodons) dataType;
        else
            throw new XMLParseException("Must construct " + MARKOV_MODULATED_YANG_MODEL + " with hidden codons");

        Parameter omegaParam = (Parameter) xo.getElementFirstChild(OMEGA);
        Parameter kappaParam = (Parameter) xo.getElementFirstChild(KAPPA);
        Parameter switchingParam = (Parameter) xo.getElementFirstChild(SWITCHING_RATES);
        FrequencyModel freqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

        EigenSystem eigenSystem;
        if (xo.getAttribute(DIAGONALIZATION, "default").compareToIgnoreCase("colt") == 0)
            eigenSystem = new ColtEigenSystem(dataType.getStateCount());
        else
            eigenSystem = new DefaultEigenSystem(dataType.getStateCount());

        return new MarkovModulatedGY94CodonModel(codons, switchingParam, omegaParam, kappaParam, freqModel, eigenSystem);
    }

    public String getParserDescription() {
        return "This element represents the a Markov-modulated Yang model of codon evolution.";
    }

    public Class getReturnType() {
        return MarkovModulatedGY94CodonModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            AttributeRule.newStringRule(DataType.DATA_TYPE),
            AttributeRule.newStringRule(GeneticCode.GENETIC_CODE),
            new ElementRule(OMEGA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SWITCHING_RATES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(FrequencyModel.class),
            AttributeRule.newStringRule(DIAGONALIZATION),
    };

}
