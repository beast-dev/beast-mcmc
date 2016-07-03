/*
 * EmpiricalCodonModelParser.java
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


import java.util.logging.Logger;

import dr.evomodel.substmodel.codon.EmpiricalCodonModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneticCode;
import dr.evomodel.substmodel.EmpiricalRateMatrixReader;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * XML parser for ECM
 * 
 * @author Stefan Zoller
 */
public class EmpiricalCodonModelParser extends AbstractXMLObjectParser {
	public static final String EMPIRICAL_CODON_MODEL = "empiricalCodonModel";
    public static final String EMPIRICAL_RATE_MATRIX = "empiricalRateMatrix";
    public static final String ECM_DATA_DIR = "ecmDataDir";
    public static final String ECM_DATA_MATRIX = "ecmRateFile";
    public static final String ECM_FREQ_MATRIX = "ecmFreqFile";
    public static final String OMEGA = "omega";
    public static final String KAPPATSTV = "kappaTsTv";
    public static final String MULTI_NT_CHANGE = "multiNtChange";

    public String getParserName() { return EMPIRICAL_CODON_MODEL; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Codons codons = Codons.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
                codons = Codons.UNIVERSAL;
            } else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
                codons = Codons.VERTEBRATE_MT;
            } else if (codeStr.equals(GeneticCode.YEAST.getName())) {
                codons = Codons.YEAST;
            } else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
                codons = Codons.MOLD_PROTOZOAN_MT;
            } else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
                codons = Codons.INVERTEBRATE_MT;
            } else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
                codons = Codons.CILIATE;
            } else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
                codons = Codons.ECHINODERM_MT;
            } else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
                codons = Codons.EUPLOTID_NUC;
            } else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
                codons = Codons.BACTERIAL;
            } else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
                codons = Codons.ALT_YEAST;
            } else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
                codons = Codons.ASCIDIAN_MT;
            } else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
                codons = Codons.FLATWORM_MT;
            } else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
                codons = Codons.BLEPHARISMA_NUC;
            } else if (codeStr.equals(GeneticCode.NO_STOPS.getName())) {
                codons = Codons.NO_STOPS;
            }
        }
        
        Parameter omegaParam = (Parameter)xo.getElementFirstChild(OMEGA);
        Parameter kappaParam = null;
        Parameter mntParam = null;
        if(xo.hasChildNamed(KAPPATSTV)) {
        	kappaParam = (Parameter)xo.getElementFirstChild(KAPPATSTV);
        	if(kappaParam.getDimension() != 2 && kappaParam.getDimension() != 9) {
            	throw new XMLParseException("If you use the kappa parameter, you need to enter exactly\n" +
            			"two values for ts and tv or nine values\n" +
            			"according to the Kosiol ECM+F+omega+9k model");
            }
    	} else {
    		mntParam = (Parameter)xo.getElementFirstChild(MULTI_NT_CHANGE);
    	}
    	
        String dirString = xo.getStringAttribute(ECM_DATA_DIR);
        String freqString = xo.getStringAttribute(ECM_FREQ_MATRIX);
        String matString = xo.getStringAttribute(ECM_DATA_MATRIX);
        
        EmpiricalRateMatrixReader rateMat = new EmpiricalRateMatrixReader(EMPIRICAL_RATE_MATRIX, codons,
        													dirString, freqString, matString);

        FrequencyModel freqModel = null;
        if (xo.getChild(FrequencyModel.class) != null) {
        	freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);
        } else {
        	freqModel = createNewFreqModel(codons, rateMat);
        }

        return new EmpiricalCodonModel(codons, omegaParam, kappaParam, mntParam, rateMat, freqModel);
    }
    
    // creates new FrequencyModel from XML frequencies
    private FrequencyModel createNewFreqModel(DataType codons, EmpiricalRateMatrixReader type) throws XMLParseException {
    	double[] freqs = type.getFrequencies();
    	double sum = 0;
        for (int j = 0; j < freqs.length; j++) {
            sum += freqs[j];
        }

        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new XMLParseException("Frequencies do not sum to 1 (they sum to " + sum + ")");
        }
        
    	FrequencyModel fm = new FrequencyModel(codons, freqs);
    	
    	Logger.getLogger("dr.evomodel").info("Using frequencies from data file");
    	
    	return fm;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the empirical model of codon evolution.";
    }

    public Class getReturnType() { return EmpiricalCodonModel.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new StringAttributeRule(GeneticCode.GENETIC_CODE,
            "The genetic code to use",
            new String[] {
                GeneticCode.UNIVERSAL.getName(),
                GeneticCode.VERTEBRATE_MT.getName(),
                GeneticCode.YEAST.getName(),
                GeneticCode.MOLD_PROTOZOAN_MT.getName(),
                GeneticCode.INVERTEBRATE_MT.getName(),
                GeneticCode.CILIATE.getName(),
                GeneticCode.ECHINODERM_MT.getName(),
                GeneticCode.EUPLOTID_NUC.getName(),
                GeneticCode.BACTERIAL.getName(),
                GeneticCode.ALT_YEAST.getName(),
                GeneticCode.ASCIDIAN_MT.getName(),
                GeneticCode.FLATWORM_MT.getName(),
                GeneticCode.BLEPHARISMA_NUC.getName(),
                GeneticCode.NO_STOPS.getName()}, true),
        new StringAttributeRule(ECM_DATA_DIR,
        	"The directory with the ECM data file",
            "ecmdata", true),
        new StringAttributeRule(ECM_DATA_MATRIX,
        	"The file with the ECM data matrix",
            "matrix.csv", true),
        new StringAttributeRule(ECM_FREQ_MATRIX,
            "The csv file with the ECM frequency matrix",
            "freqs.csv", true),
        new ElementRule(OMEGA,
        	new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
        new XORRule(
        	new ElementRule(KAPPATSTV,
        		new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
        	new ElementRule(MULTI_NT_CHANGE,
        	    new XMLSyntaxRule[] { new ElementRule(Parameter.class) })),
        new ElementRule(FrequencyModel.class, true)
    };
}

