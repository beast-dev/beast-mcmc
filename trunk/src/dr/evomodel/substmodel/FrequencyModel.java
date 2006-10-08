/*
 * FrequencyModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.*;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A model of equlibrium frequencies
 *
 * @version $Id: FrequencyModel.java,v 1.26 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class FrequencyModel extends AbstractModel {

	public static final String FREQUENCIES = "frequencies";
	public static final String FREQUENCY_MODEL = "frequencyModel";
	
	public FrequencyModel(DataType dataType, Parameter frequencyParameter) {  
	
		super(FREQUENCY_MODEL);
		this.frequencyParameter = frequencyParameter;
		addParameter(frequencyParameter);
		frequencyParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, dataType.getStateCount()));
		this.dataType = dataType;
	}
	
	public final void setFrequency(int i, double value) {
		frequencyParameter.setParameterValue(i, value);
	}
	
	public final double getFrequency(int i) { return frequencyParameter.getParameterValue(i); }
	
	public int getFrequencyCount() { return frequencyParameter.getDimension(); }
	
	public final double[] getFrequencies() { 
		double[] frequencies = new double[getFrequencyCount()];
		for (int i =0; i < frequencies.length; i++) {
			frequencies[i] = getFrequency(i);		 
		}
		return frequencies;
	}
	
	public DataType getDataType() { return dataType; }
	
	// *****************************************************************
	// Interface Model
	// *****************************************************************
		
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no intermediates need recalculating....
	}
	
	protected void handleParameterChangedEvent(Parameter parameter, int index) {
		// no intermediates need recalculating....
	}
	
	protected void storeState() {} // no state apart from parameters to store 
	protected void restoreState() {} // no state apart from parameters to restore 
	protected void acceptState() {} // no state apart from parameters to accept 
	protected void adoptState(Model source) {} // no state apart from parameters to adopt 
	
	public Element createElement(Document doc) {
		throw new RuntimeException("Not implemented!");
	}

	/**
	 * Reads a frequency model from an XMLObject. 
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return FREQUENCY_MODEL; }
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
				
			FrequencyModel freqModel = null;
			DataType dataType = null;
				
			if (xo.hasAttribute(DataType.DATA_TYPE)) {
				String dataTypeStr = xo.getStringAttribute(DataType.DATA_TYPE);
				if (dataTypeStr.equals(Nucleotides.DESCRIPTION)) {
					dataType = Nucleotides.INSTANCE;
				} else if (dataTypeStr.equals(HiddenNucleotides.DESCRIPTION)) {
					dataType = HiddenNucleotides.INSTANCE;
				} else if (dataTypeStr.equals(AminoAcids.DESCRIPTION)) {
					dataType = AminoAcids.INSTANCE;
				} else if (dataTypeStr.equals(Codons.DESCRIPTION)) {
					dataType = Codons.UNIVERSAL;
				} else if (dataTypeStr.equals(TwoStates.DESCRIPTION)) {
					dataType = TwoStates.INSTANCE;
				}
			}
							
			if (dataType == null) dataType = (DataType)xo.getChild(DataType.class);
			
			if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
				String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
				if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
					dataType = Codons.UNIVERSAL;
				} else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
					dataType = Codons.VERTEBRATE_MT;
				} else if (codeStr.equals(GeneticCode.YEAST.getName())) {
					dataType = Codons.YEAST;
				} else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
					dataType = Codons.MOLD_PROTOZOAN_MT;
				} else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
					dataType = Codons.INVERTEBRATE_MT;
				} else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
					dataType = Codons.CILIATE;
				} else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
					dataType = Codons.ECHINODERM_MT;
				} else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
					dataType = Codons.EUPLOTID_NUC;
				} else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
					dataType = Codons.BACTERIAL;
				} else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
					dataType = Codons.ALT_YEAST;
				} else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
					dataType = Codons.ASCIDIAN_MT;
				} else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
					dataType = Codons.FLATWORM_MT;
				} else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
					dataType = Codons.BLEPHARISMA_NUC;
				}
			}
			
			Parameter freqsParam = (Parameter)xo.getSocketChild(FREQUENCIES);
			double[] frequencies = null;	
						
			for (int i =0; i < xo.getChildCount(); i++) {
				Object obj = xo.getChild(i);
				if (obj instanceof PatternList) {
					frequencies = ((PatternList)obj).getStateFrequencies();
				} 
			}
			
			if (frequencies != null) {
				if (freqsParam.getDimension() != frequencies.length) {
					throw new XMLParseException("dimension of frequency parameter and number of sequence states don't match!");
				}
				for (int j = 0; j < frequencies.length; j++) {
					freqsParam.setParameterValue(j, frequencies[j]);
				}
			}
			
			freqModel = new FrequencyModel(dataType, freqsParam);			
			return freqModel;
		}
		
		public String getParserDescription() {
			return "A model of equilibrium base frequencies.";
		}

		public Class getReturnType() { return FrequencyModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new XORRule(
				new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data", new String[] { Nucleotides.DESCRIPTION, HiddenNucleotides.DESCRIPTION, AminoAcids.DESCRIPTION, Codons.DESCRIPTION, TwoStates.DESCRIPTION}, false),
				new ElementRule(DataType.class)
				),
			new ElementRule(FREQUENCIES, 
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) })
		
		};
};
	
	private DataType dataType = null;
	private Parameter frequencyParameter = null;
	
}
