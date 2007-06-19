/*
 * RelaxedClockHKY.java
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

import dr.evolution.datatype.HiddenNucleotides;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A model with hidden states that represent different rates.
 *
 * @version $Id: RelaxedClockHKY.java,v 1.4 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class RelaxedClockHKY extends AbstractSubstitutionModel
{
	public static final String RELAXED_CLOCK_HKY = "relaxedHKY";
	public static final String BETA = "beta";
	public static final String KAPPA_1 = "kappa1";
	public static final String KAPPA_2 = "kappa2";
	public static final String RATE_RATIO = "rateRatio";
	public static final String FREQUENCIES = "frequencies";
	
	/** kappa */
	private Parameter kappaParameter;
	private Parameter kappa2Parameter;
	
	/** beta - the rate of change of type */
	private Parameter betaParameter;
	
	/** fast rate parameter */
	private Parameter rateRatioParameter;
	
	
	byte[] rateMap;
	
	/**
	 * Constructor
	 */
	public RelaxedClockHKY(Parameter kappaParameter, Parameter kappa2Parameter, Parameter betaParameter, Parameter rateRatioParameter, FrequencyModel freqModel) {
		super(RELAXED_CLOCK_HKY, HiddenNucleotides.INSTANCE, freqModel);
		
		this.kappaParameter = kappaParameter;
		this.kappa2Parameter = kappa2Parameter;
		this.betaParameter = betaParameter;
		this.rateRatioParameter = rateRatioParameter;
		addParameter(kappaParameter);
		addParameter(kappa2Parameter);
		addParameter(betaParameter);
		addParameter(rateRatioParameter);
		constructRateMap();
		setupRelativeRates();
	}
	
	public void frequenciesChanged() {
		// DO NOTHING
	}
	
	public void ratesChanged() {
		//constructRateMap();
		setupRelativeRates();
	}
	
	/**
	 * set kappa
	 */
	public void setKappa(double kappa) { kappaParameter.setParameterValue(0, kappa); }
	
	/**
	 * @return kappa
	 */
	public double getKappa() { return kappaParameter.getParameterValue(0); }
	
	/**
	 * setup substitution matrix
	 */
	public void setupRelativeRates()
	{
		double kappa = getKappa();
		double kappa2 = kappa2Parameter.getParameterValue(0);
		double beta = betaParameter.getParameterValue(0);
		double rateRatio = rateRatioParameter.getParameterValue(0);
		
		for (int i = 0; i < rateCount; i++) {
			switch(rateMap[i]) {
				case 0: relativeRates[i] = beta; break;				// rate of the site changes
				case 1: relativeRates[i] = 0.0; break;				// rate and state change
				case 2: relativeRates[i] = kappa; break;			// slow transition
				case 3: relativeRates[i] = kappa2 * rateRatio; break; 			// fast transition
				case 4: relativeRates[i] = 1.0; break;				// slow transversion
				case 5: relativeRates[i] = rateRatio; break;			// fast transversion
			}
		}
	}
	
	/**
	 * Construct a map of the rate classes in the rate matrix:
	 *		0: rate change only
	 *		1: rate change *and* nucleotide change
	 *		2: slow transition
	 *		3: fast transition
	 *		4: slow transversion
	 *		5: fast transversion
	 */
	private void constructRateMap() {
	
		byte rateClass;
		int fromNuc, toNuc;
		int fromRate, toRate;
		int count = 0;
		
		rateMap = new byte[rateCount];
				
		for (int i = 0; i < stateCount; i++) {
		
			for (int j = i + 1; j < stateCount; j++) {
			
				rateClass = -1;
				
				fromNuc = i % 4;
				toNuc = j % 4;
				fromRate = i / 4;
				toRate =  j / 4;
				
				if (fromNuc == toNuc) {
					// rate transition
					if (fromRate == toRate) {
						throw new RuntimeException("Shouldn't be possible");
					}
					rateClass = 0;
				} else if (fromRate != toRate) {
					rateClass = 1;
				} else if ((fromNuc == 0 && toNuc == 2) || (fromNuc == 1 && toNuc == 3) || (fromNuc == 2 && toNuc == 0) || (fromNuc == 3 && toNuc == 1)) {
					// transition
					if (fromRate == 0) {
						rateClass = 2;
					} else {
						rateClass = 3;
					}
				} else {
					// transversion
					if (fromRate == 0) {
						rateClass = 4;
					} else {
						rateClass = 5;
					}
				}
				
				rateMap[count] = rateClass;
				count++;
			}	
		}
	}
	
	// Normalize rate matrix to one expected substitution per unit time
	public double getHiddenRate() {
		
		double[][] matrix = new double[4*2][4*2];
		// Set the instantaneous rate matrix
		int k = 0;
		for (int i = 0; i < 8; i++) {
			for (int j = i+1; j < 8; j++) {
				matrix[i][j] = relativeRates[k] * freqModel.getFrequency(j);
				matrix[j][i] = relativeRates[k] * freqModel.getFrequency(i);
				//System.out.println(amat[i][j]);
				k++;
			}
		}
		makeValid(matrix, stateCount);
		normalize(matrix, freqModel.getFrequencies(), stateCount);
		
		double[] vector = freqModel.getFrequencies();
		
		double hiddenRate = 0.0;
		for (int i = 0; i < 4; i++) {
			hiddenRate += matrix[i][i+4]*vector[i+4];
			hiddenRate += matrix[i+4][i]*vector[i];
		}
		
		return hiddenRate;
	}

	/**
	 * Parses an element from an DOM document into a DemographicModel. Recognises
	 * ConstantPopulation and ExponentialGrowth.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return RELAXED_CLOCK_HKY; }
		
		public String getParserDescription() {
			return "A relaxed HKY model.";
		}
		
		public Class getReturnType() { return SubstitutionModel.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(KAPPA_1, new XMLSyntaxRule[] {new ElementRule(Parameter.class, "A parameter representing the transition transversion bias")}),
			new ElementRule(KAPPA_2, new XMLSyntaxRule[] {
				new ElementRule(Parameter.class, "A parameter representing the transition transversion bias of the faster class")}),
			new ElementRule(BETA, new XMLSyntaxRule[] {new ElementRule(Parameter.class, "A parameter representing the rate of change between the two classes")}),
			new ElementRule(RATE_RATIO, new XMLSyntaxRule[] {new ElementRule(Parameter.class, "A parameter representing relative speed up of the second class")})
		};
			
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
				
			Parameter kappaParam = null;
			Parameter kappa2Param = null;
			Parameter betaParam = null;
			Parameter rateRatioParam = null;
			FrequencyModel freqModel = null;
					
			if (xo.getChildCount() < 4) {
				throw new XMLParseException("Expected five children in relaxedHKYModel, got " + xo.getChildCount()); 
			}
			
			
			kappaParam = (Parameter)xo.getSocketChild(KAPPA_1);
			kappa2Param = (Parameter)xo.getSocketChild(KAPPA_2);
			betaParam = (Parameter)xo.getSocketChild(BETA);
			rateRatioParam = (Parameter)xo.getSocketChild(RATE_RATIO);
			freqModel = (FrequencyModel)xo.getSocketChild(FREQUENCIES);
			
			for (int i = 0; i< freqModel.getFrequencyCount(); i++) {
				System.out.println(HiddenNucleotides.INSTANCE.getChar(i) +"\t"+ freqModel.getFrequency(i));
			}

			if (!(freqModel.getDataType() instanceof HiddenNucleotides)) throw new IllegalArgumentException("Datatype must be hidden nucleotides!!");
				
			return new RelaxedClockHKY(kappaParam, kappa2Param, betaParam, rateRatioParam, freqModel);
		}
	};
}