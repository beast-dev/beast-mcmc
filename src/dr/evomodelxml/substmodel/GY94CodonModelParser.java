/*
 * GY94CodonModelParser.java
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

import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class GY94CodonModelParser extends AbstractXMLObjectParser {

    public static final String YANG_CODON_MODEL = "yangCodonModel";
    public static final String OMEGA = "omega";
    public static final String KAPPA = "kappa";


    public String getParserName() { return YANG_CODON_MODEL; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Codons codons = Codons.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            codons = Codons.findByName(codeStr);
        }

        Parameter omegaParameter = (Parameter)xo.getElementFirstChild(OMEGA);
        
        int dim = omegaParameter.getDimension();
        double value = omegaParameter.getParameterValue(dim - 1); 
        if(value < 0) {
        	throw new RuntimeException("Negative Omega parameter value " + value);
        }//END: negative check
        
        Parameter kappaParameter = (Parameter)xo.getElementFirstChild(KAPPA);
        
        dim = kappaParameter.getDimension();
        value = kappaParameter.getParameterValue(dim - 1);
        if(value < 0) {
        	throw new RuntimeException("Negative kappa parameter value value " + value);
        }//END: negative check
        
        FrequencyModel freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);
        GY94CodonModel codonModel = new GY94CodonModel(codons, omegaParameter, kappaParameter, freqModel);

//            codonModel.printRateMap();

        return codonModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the Yang model of codon evolution.";
    }

    public Class getReturnType() { return GY94CodonModel.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new StringAttributeRule(GeneticCode.GENETIC_CODE,
                    "The genetic code to use",
                    GeneticCode.GENETIC_CODE_NAMES, true),
            new ElementRule(OMEGA,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
            new ElementRule(FrequencyModel.class)
    };

}