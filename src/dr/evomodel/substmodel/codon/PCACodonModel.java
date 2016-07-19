/*
 * PCACodonModel.java
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

package dr.evomodel.substmodel.codon;

import dr.evolution.datatype.Codons;
import dr.evomodelxml.substmodel.PCACodonModelParser;
import dr.evomodel.substmodel.*;
import dr.inference.model.Parameter;

/**
 * PCA model of codon evolution
 * 
 * @author Stefan Zoller
 */
public class PCACodonModel extends BaseSubstitutionModel {
	// principal components, means, scale factors
	protected AbstractPCARateMatrix rateMatrix;

	protected byte[] rateMap;
	
	protected Parameter pcFactors;

	/**
     * constructors
     *
     * @param codonDataType				Data type as Codons.UNIVERSAL
     * @param pcaType           		Rate matrix with PCs, means, scale factors
     * @param pcaDimensionParameter		Scalars for PCs
     * @param freqModel					Frequency model
     */
	public PCACodonModel(Codons codonDataType, AbstractPCARateMatrix pcaType, Parameter pcaDimensionParameter,
            FrequencyModel freqModel) {
		this(codonDataType, pcaType, pcaDimensionParameter, freqModel,
				new DefaultEigenSystem(codonDataType.getStateCount()));
	}
	
	
	public PCACodonModel(Codons codonDataType,
						    AbstractPCARateMatrix pcaType,
						    Parameter pcaDimensionParameter,
						    FrequencyModel freqModel,
						    EigenSystem eigenSystem)
	{
		super(PCACodonModelParser.PCA_CODON_MODEL, codonDataType, freqModel, eigenSystem);
		
		this.rateMatrix = pcaType;
		
		this.pcFactors = pcaDimensionParameter;
		addVariable(pcFactors);
		pcFactors.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
				pcFactors.getDimension()));
		
		// initialize scalars for principal components
		double[] startFacs = pcaType.getStartFacs();
		double facSum = 0.0;
		for(int i=0; i<pcFactors.getDimension(); i++) {
			facSum += startFacs[i];
		}
		for(int i=0; i<pcFactors.getDimension(); i++) {
			pcFactors.setParameterValueQuietly(i, startFacs[i]/facSum);
		}
	}

    
	// setup substitution matrix
    public void setupRelativeRates(double[] rates) {
        double[] m = rateMatrix.getMeans();
        double[] sc = rateMatrix.getScales();
        for (int i = 0; i < rateCount; i++) {
            rates[i] = m[i];
        }
        for(int j = 0; j < pcFactors.getDimension(); j++) {
        	double[] pc = rateMatrix.getPCAt(j);
        	double factor = getPcFactor(j);
        	for (int i = 0; i < rateCount; i++) {
                rates[i] += factor*pc[i]*sc[i];
            }
        }
        for (int i = 0; i < rateCount; i++) {
        	if(rates[i] < Double.MIN_VALUE) {
        		rates[i] = Double.MIN_VALUE;
        	}
        }
        return;
    }
    
    protected void ratesChanged() {
	}
    
    protected void frequenciesChanged() {
    }
    
    /**
     * Getter and setter for Parameter pcFactor
     */

    public double getPcFactor(int dim) {
        return pcFactors.getParameterValue(dim);
    }
    
    public double[] getPcFactor() {
    	return pcFactors.getParameterValues();
    }
    
    public void setPcFactor(int dim, double fac) {
		pcFactors.setParameterValue(dim, fac);
		updateMatrix = true;
	}
    
    public void setPcFactor(double[] fac) {
    	for(int i=0; i<pcFactors.getDimension(); i++) {
    		pcFactors.setParameterValue(i, fac[i]);
    	}
    	updateMatrix = true;
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

	public String toXHTML() {
		StringBuffer buffer = new StringBuffer();

		buffer.append("<em>PCA Codon Model</em>");

		return buffer.toString();
	}
}
