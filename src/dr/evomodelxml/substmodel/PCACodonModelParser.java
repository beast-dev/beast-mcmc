/*
 * PCACodonModelParser.java
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

import dr.evomodel.substmodel.PCARateMatrixMammalia;
import dr.evomodel.substmodel.codon.PCACodonModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneticCode;
import dr.evomodel.substmodel.AbstractPCARateMatrix;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * XML parser for PCACodonModel
 * 
 * @author Stefan Zoller
 */
public class PCACodonModelParser extends AbstractXMLObjectParser {
    public static final String PCA_CODON_MODEL = "pcaCodonModel";
    public static final String PCATYPE = "pcaType";
    public static final String PCANUMBER = "pcaNumber";
    public static final String PCA_DIMENSION = "pcaDimension";
    public static final String PCA_DATA_DIR = "pcaDataDir";

    public String getParserName() { return PCA_CODON_MODEL; }

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
        
        // get number of PCs
        Parameter pcaDimensionParameter = (Parameter)xo.getElementFirstChild(PCA_DIMENSION);
        
        // get directory with pca rate matrix files; fallback to default "pcadata"
        String dirString = "pcadata";
        if(xo.hasAttribute(PCA_DATA_DIR)) {
        	dirString = xo.getStringAttribute(PCA_DATA_DIR);
        }
        
        // get type of rate matrix; fallback to mammalia pca
        AbstractPCARateMatrix pcaType = new PCARateMatrixMammalia(pcaDimensionParameter.getDimension(), dirString);
        // check for other type of pca
        if(xo.hasAttribute(PCATYPE)) {
            String pcaTypeString = xo.getStringAttribute(PCATYPE);
            if(pcaTypeString.equals(PCARateMatrixMammalia.getName())) {
                pcaType = new PCARateMatrixMammalia(pcaDimensionParameter.getDimension(), dirString);
            }
        }

        // decide if getting frequencies from csv or estimating from MSA
        FrequencyModel freqModel = null;
        if (xo.getChild(FrequencyModel.class) != null) {
        	freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);
        } else {
        	freqModel = createNewFreqModel(codons, pcaType);
        }

        return new PCACodonModel(codons, pcaType, pcaDimensionParameter, freqModel);
    }
    
    // read frequencies from XML and return FrequencyModel object
    private FrequencyModel createNewFreqModel(DataType codons, AbstractPCARateMatrix type) throws XMLParseException {
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
        return "This element represents the PCA model of codon evolution.";
    }

    public Class getReturnType() { return PCACodonModel.class; }

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
        new StringAttributeRule(PCATYPE,
            "The PCA type to use",
            new String[] {
                PCARateMatrixMammalia.getName(),
            }, true),
        new StringAttributeRule(PCA_DATA_DIR,
        	"The directory with the PCA data files",
            "pcadata", true),
        new ElementRule(PCA_DIMENSION,
        	new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
        new ElementRule(FrequencyModel.class, true)
    };
}

