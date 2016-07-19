/*
 * EmpiricalCodonModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.oldevomodelxml.substmodel.EmpiricalCodonModelParser;
import dr.inference.model.Parameter;
import java.util.logging.Logger;

/**
 * Empirical model of codon evolution
 * 
 * @author Stefan Zoller
 */
public class EmpiricalCodonModel extends AbstractCodonModel {
	
	protected byte[] rateMap;
	
	private Parameter omegaParameter;
	private Parameter kappaParameter;	// 2d (kappats and kappatv) or 9d
	private Parameter multintParameter;
	private EmpiricalCodonRateMatrix rateMat;
	
	private int modelType;
	private final int ECM_OMEGA_2K = 2;
	private final int ECM_OMEGA_9K = 3;
	private final int ECM_OMEGA_NU = 4;
	private final int ECM_OMEGA = 1;
	

	/**
     * constructor
     *
     * @param codonDataType				Data type as Codons.UNIVERSAL
     * @param omegaParam				Parameter: Omega
     * @param kappaParam				Parameter: Kappa (multidimensional)
     * @param mntParam					Parameter: Multi-nt
     * @param rMat						Initial rate matrix and frequencies
     * @param freqModel					Frequency model
     */
	public EmpiricalCodonModel(Codons codonDataType,
						    Parameter omegaParam,
						    Parameter kappaParam,
						    Parameter mntParam,
						    EmpiricalCodonRateMatrix rMat,
						    FrequencyModel freqModel)
	{
		super(EmpiricalCodonModelParser.EMPIRICAL_CODON_MODEL, codonDataType, freqModel);
		
		// setup parameters
		this.omegaParameter = omegaParam;
		addVariable(omegaParameter);
		omegaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
				omegaParameter.getDimension()));
		
		if(kappaParam != null) {
			this.kappaParameter = kappaParam;
			addVariable(kappaParameter);
			kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
				kappaParameter.getDimension()));
		}
		if(mntParam != null) {
			this.multintParameter = mntParam;
			addVariable(multintParameter);
			multintParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
				multintParameter.getDimension()));
		}
		
		this.rateMat = rMat;

		constructRateMap();
		
		checkForModelType();
	}
	
	// decide which model to use: ECM_OMEGA_2K, ECM_OMEGA_9K, ECM_OMEGA_NU or ECM_OMEGA
	private void checkForModelType() {
		this.modelType = 0;
		if(kappaParameter != null) {
			if(kappaParameter.getDimension() == 2) {
				this.modelType = ECM_OMEGA_2K;
				Logger.getLogger("dr.evomodel").info("Using model ECM+omega+2k");
			} else {
				this.modelType = ECM_OMEGA_9K;
				Logger.getLogger("dr.evomodel").info("Using model ECM+omega+9k");
			}
		}
		if(multintParameter != null){
			this.modelType = ECM_OMEGA_NU;
			Logger.getLogger("dr.evomodel").info("Using model ECM+omega+nu");
		}
		if(kappaParameter == null && multintParameter == null) {
			this.modelType = ECM_OMEGA;
			Logger.getLogger("dr.evomodel").info("Using model ECM+omega");
		}
	}

    
	// setup substitution matrix depending on model type
	public void setupRelativeRates() {
		switch(modelType) {
		case ECM_OMEGA:
			setupRelativeRatesECMOmega();
			break;
		case ECM_OMEGA_2K:
			setupRelativeRatesECMOmega2k();
			break;
		case ECM_OMEGA_9K:
			setupRelativeRatesECMOmega9k();
			break;
		case ECM_OMEGA_NU:
			setupRelativeRatesECMOmegaNu();
			break;
		}
	}
	
	// actual setup routines for different models
	private void setupRelativeRatesECMOmega() {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
			case 1:													// 1ts, 0tv, syn
			case 3:													// 0ts, 1tv, syn
			case 5:													// 2ts, 0tv, syn
			case 7:													// 1ts, 1tv, syn
			case 9:													// 0ts, 2tv, syn
			case 11:												// 3ts, 0tv, syn
			case 13:												// 2ts, 1tv, syn
			case 15:												// 1ts, 2tv, syn
			case 17: relativeRates[i] = initRateMatrix[i]; break;	// 0ts, 3tv, syn
			
			case 2:															// 1ts, 0tv, nonsyn
			case 4:															// 0ts, 1tv, nonsyn
			case 6:															// 2ts, 0tv, nonsyn
			case 8:															// 1ts, 1tv, nonsyn
			case 10:														// 0ts, 2tv, nonsyn
			case 12:														// 3ts, 0tv, nonsyn
			case 14:														// 2ts, 1tv, nonsyn
			case 16:														// 1ts, 2tv, nonsyn
			case 18: relativeRates[i] = initRateMatrix[i] * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
	private void setupRelativeRatesECMOmega2k() {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		double kts = getKappaTs();
		double ktv = getKappaTv();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
				case 1: relativeRates[i] = initRateMatrix[i] * kts; break;						// 1ts, 0tv, syn
				case 2: relativeRates[i] = initRateMatrix[i] * kts * omega; break;				// 1ts, 0tv, nonsyn
				case 3: relativeRates[i] = initRateMatrix[i] * ktv; break;						// 0ts, 1tv, syn
				case 4: relativeRates[i] = initRateMatrix[i] * ktv * omega; break;				// 0ts, 1tv, nonsyn
				case 5: relativeRates[i] = initRateMatrix[i] * kts * kts; break;				// 2ts, 0tv, syn
				case 6: relativeRates[i] = initRateMatrix[i] * kts * kts * omega; break;		// 2ts, 0tv, nonsyn
				case 7: relativeRates[i] = initRateMatrix[i] * kts * ktv; break;				// 1ts, 1tv, syn
				case 8: relativeRates[i] = initRateMatrix[i] * kts * ktv * omega; break;		// 1ts, 1tv, nonsyn
				case 9: relativeRates[i] = initRateMatrix[i] * ktv * ktv; break;				// 0ts, 2tv, syn
				case 10: relativeRates[i] = initRateMatrix[i] * ktv * ktv * omega; break;		// 0ts, 2tv, nonsyn
				case 11: relativeRates[i] = initRateMatrix[i] * kts * kts * kts; break;			// 3ts, 0tv, syn
				case 12: relativeRates[i] = initRateMatrix[i] * kts * kts * kts * omega; break;	// 3ts, 0tv, nonsyn
				case 13: relativeRates[i] = initRateMatrix[i] * kts * kts * ktv; break;			// 2ts, 1tv, syn
				case 14: relativeRates[i] = initRateMatrix[i] * kts * kts * ktv * omega; break;	// 2ts, 1tv, nonsyn
				case 15: relativeRates[i] = initRateMatrix[i] * kts * ktv * ktv; break;			// 1ts, 2tv, syn
				case 16: relativeRates[i] = initRateMatrix[i] * kts * ktv * ktv * omega; break;	// 1ts, 2tv, nonsyn
				case 17: relativeRates[i] = initRateMatrix[i] * ktv * ktv * ktv; break;			// 0ts, 3tv, syn
				case 18: relativeRates[i] = initRateMatrix[i] * ktv * ktv * ktv * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
	private void setupRelativeRatesECMOmega9k() {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		double[] kappa = getKappa();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
				case 1: relativeRates[i] = initRateMatrix[i] * kappa[0]; break;						// 1ts, 0tv, syn
				case 2: relativeRates[i] = initRateMatrix[i] * kappa[0] * omega; break;				// 1ts, 0tv, nonsyn
				case 3: relativeRates[i] = initRateMatrix[i] * kappa[1]; break;						// 0ts, 1tv, syn
				case 4: relativeRates[i] = initRateMatrix[i] * kappa[1] * omega; break;				// 0ts, 1tv, nonsyn
				case 5: relativeRates[i] = initRateMatrix[i] * kappa[2]; break;				// 2ts, 0tv, syn
				case 6: relativeRates[i] = initRateMatrix[i] * kappa[2] * omega; break;		// 2ts, 0tv, nonsyn
				case 7: relativeRates[i] = initRateMatrix[i] * kappa[3]; break;				// 1ts, 1tv, syn
				case 8: relativeRates[i] = initRateMatrix[i] * kappa[3] * omega; break;		// 1ts, 1tv, nonsyn
				case 9: relativeRates[i] = initRateMatrix[i] * kappa[4]; break;				// 0ts, 2tv, syn
				case 10: relativeRates[i] = initRateMatrix[i] * kappa[4] * omega; break;		// 0ts, 2tv, nonsyn
				case 11: relativeRates[i] = initRateMatrix[i] * kappa[5]; break;			// 3ts, 0tv, syn
				case 12: relativeRates[i] = initRateMatrix[i] * kappa[5] * omega; break;	// 3ts, 0tv, nonsyn
				case 13: relativeRates[i] = initRateMatrix[i] * kappa[6]; break;			// 2ts, 1tv, syn
				case 14: relativeRates[i] = initRateMatrix[i] * kappa[6] * omega; break;	// 2ts, 1tv, nonsyn
				case 15: relativeRates[i] = initRateMatrix[i] * kappa[7]; break;			// 1ts, 2tv, syn
				case 16: relativeRates[i] = initRateMatrix[i] * kappa[7] * omega; break;	// 1ts, 2tv, nonsyn
				case 17: relativeRates[i] = initRateMatrix[i] * kappa[8]; break;			// 0ts, 3tv, syn
				case 18: relativeRates[i] = initRateMatrix[i] * kappa[8] * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
	private void setupRelativeRatesECMOmegaNu() {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		double mnt = getMultiNt();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
			case 1:																// 1ts, 0tv, syn
			case 3: relativeRates[i] = initRateMatrix[i]; break;				// 0ts, 1tv, syn
			case 2: 															// 1ts, 0tv, nonsyn
			case 4: relativeRates[i] = initRateMatrix[i] * omega; break;		// 0ts, 1tv, nonsyn
			
			case 5:																// 2ts, 0tv, syn
			case 7: 															// 1ts, 1tv, syn
			case 9: 															// 0ts, 2tv, syn
			case 11: 															// 3ts, 0tv, syn
			case 13: 															// 2ts, 1tv, syn
			case 15: 															// 1ts, 2tv, syn
			case 17: relativeRates[i] = initRateMatrix[i] * mnt; break;			// 0ts, 3tv, syn
			
			case 6: 															// 2ts, 0tv, nonsyn
			case 8: 															// 1ts, 1tv, nonsyn
			case 10: 															// 0ts, 2tv, nonsyn
			case 12: 															// 3ts, 0tv, nonsyn
			case 14: 															// 2ts, 1tv, nonsyn
			case 16: 															// 1ts, 2tv, nonsyn
			case 18: relativeRates[i] = initRateMatrix[i] * mnt * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
    
    protected void ratesChanged() {
	}
    
    protected void frequenciesChanged() {
	}
    
    // getter and setter for parameters

    public void setOmega(double omega) {
		omegaParameter.setParameterValue(0, omega);
		updateMatrix = true;
	}

    public double getOmega() { return omegaParameter.getParameterValue(0); }
    
    
    public void setKappa(double kts, double ktv) {
    	if(kappaParameter != null) {
    		kappaParameter.setParameterValue(0, kts);
    		kappaParameter.setParameterValue(1, ktv);
    		updateMatrix = true;
    	}
	}
    
    public double getKappaTs() { 
    	if(kappaParameter != null) {
    		return kappaParameter.getParameterValue(0);
    	} else {
    		return 0.0;
    	}
    }
    
    public double getKappaTv() { 
    	if(kappaParameter != null) {
    		return kappaParameter.getParameterValue(1);
    	} else {
    		return 0.0;
    	}
    }
    
    public double[] getKappa() {
    	if(kappaParameter != null) {
    		return kappaParameter.getParameterValues();
    	} else {
    		return new double[9];
    	}
    }
    
    
    public void setMultiNt(double mnt) {
    	if(multintParameter != null) {
    		multintParameter.setParameterValue(0, mnt);
    		updateMatrix = true;
    	}
	}

    public double getMultiNt() { 
    	if(multintParameter != null) {
    		return multintParameter.getParameterValue(0);
    	} else {
    		return 0.0;
    	}
	}
    
    
    /**
	 * Construct a map of the rate classes in the rate matrix using the current
	 * genetic code. Classes are:
	 * 	1-2: 1ts, 0tv (syn/nonsyn)
	 *  3-4: 0ts, 1tv
	 *  5-6: 2ts, 0tv
	 *  7-8: 1ts, 1tv
	 *  9-10: 0ts, 2tv
	 *  11-12: 3ts, 0tv
	 *  13-14: 2ts, 1tv
	 *  15-16: 1ts, 2tv
	 *  17-18: 0ts, 3tv
	 */
	protected void constructRateMap()
	{
		int u, v, i1, j1, k1, i2, j2, k2, ts, tv, non;
		byte rateClass;
		int[] codon;
		int cs1, cs2, aa1, aa2;

		int i = 0;

		rateMap = new byte[rateCount];

		for (u = 0; u < stateCount; u++) {

			codon = codonDataType.getTripletStates(u);
			i1 = codon[0];
			j1 = codon[1];
			k1 = codon[2];

			cs1 = codonDataType.getState(i1, j1, k1);
			aa1 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs1));

			for (v = u + 1; v < stateCount; v++) {
				
				ts = 0;
				tv = 0;
				non = 0;
				rateClass = -1;

				codon = codonDataType.getTripletStates(v);
				i2 = codon[0];
				j2 = codon[1];
				k2 = codon[2];

				cs2 = codonDataType.getState(i2, j2, k2);
				aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));

				if (i1 != i2) {
					if ( (i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
						 (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1) ) { // C <-> T
						ts++; // Transition at position 1
					} else {
						tv++; // Transversion at position 1
					}
				}
				if (j1 != j2) {
					if ( (j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
						(j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1) ) { // C <-> T
						ts++; // Transition
					} else {
						tv++; // Transversion
					}
				}
				if (k1 != k2) {
					if ( (k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
						(k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1) ) { // C <-> T
						ts++; // Transition
					} else {
						tv++; // Transversion
					}
				}

	 			if (aa1 != aa2) {
					non = 1; // Is a non-synonymous change
				}

	 			// decide for rateClass
	 			switch(ts) {
	 				case 0:
	 					switch(tv) {
	 						case 1: rateClass = 3; break;	// 0ts, 1tv
							case 2: rateClass = 9; break;	// 0ts, 2tv
							case 3: rateClass = 17; break;	// 0ts, 3tv
							default: break;
	 					}
	 					break;
	 				case 1:
	 					switch(tv) {
	 						case 0: rateClass = 1; break;	// 1ts, 0tv
							case 1: rateClass = 7; break;	// 1ts, 1tv
							case 2: rateClass = 15; break;	// 1ts, 2tv
							default: break;
	 					}
	 					break;
	 				case 2:
	 					switch(tv) {
 							case 0: rateClass = 5; break;	// 2ts, 0tv
 							case 1: rateClass = 13; break;	// 2ts, 1tv
 							default: break;
	 					}
	 					break;
	 				case 3:
	 					rateClass = 11; break;	// 3ts, 0tv
	 				default: break;
	 			}
	 			
	 			if(non == 1) {
	 				rateClass += 1;
	 			}
				rateMap[i] = rateClass;
				i++;
			}

		}
	}

    public void printRateMap()
    {
        int u, v, i1, j1, k1, i2, j2, k2, ts, tv, non;
        byte rateClass;
        int[] codon;
        int cs1, cs2, aa1, aa2;

        System.out.print("\t");
        for (v = 0; v < stateCount; v++) {
            codon = codonDataType.getTripletStates(v);
            i2 = codon[0];
            j2 = codon[1];
            k2 = codon[2];

            System.out.print("\t" + Nucleotides.INSTANCE.getChar(i2));
            System.out.print(Nucleotides.INSTANCE.getChar(j2));
            System.out.print(Nucleotides.INSTANCE.getChar(k2));
        }
        System.out.println();

        System.out.print("\t");
        for (v = 0; v < stateCount; v++) {
            codon = codonDataType.getTripletStates(v);
            i2 = codon[0];
            j2 = codon[1];
            k2 = codon[2];

            cs2 = codonDataType.getState(i2, j2, k2);
            aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));
            System.out.print("\t" + AminoAcids.INSTANCE.getChar(aa2));
        }
        System.out.println();

        for (u = 0; u < stateCount; u++) {

            codon = codonDataType.getTripletStates(u);
            i1 = codon[0];
            j1 = codon[1];
            k1 = codon[2];

            System.out.print(Nucleotides.INSTANCE.getChar(i1));
            System.out.print(Nucleotides.INSTANCE.getChar(j1));
            System.out.print(Nucleotides.INSTANCE.getChar(k1));

            cs1 = codonDataType.getState(i1, j1, k1);
            aa1 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs1));

            System.out.print("\t" + AminoAcids.INSTANCE.getChar(aa1));

            for (v = 0; v < stateCount; v++) {
	
            		ts = 0;
            		tv = 0;
            		non = 0;
	            	rateClass = -1;
	            	
	                codon = codonDataType.getTripletStates(v);
	                i2 = codon[0];
	                j2 = codon[1];
	                k2 = codon[2];
	
	                cs2 = codonDataType.getState(i2, j2, k2);
	                aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));
	
	                if (i1 != i2) {
						if ( (i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
							 (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1) ) { // C <-> T
							ts += 1; // Transition at position 1
						} else {
							tv += 1; // Transversion at position 1
						}
					}
					if (j1 != j2) {
						if ( (j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
							(j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1) ) { // C <-> T
							ts += 1; // Transition
						} else {
							tv += 1; // Transversion
						}
					}
					if (k1 != k2) {
						if ( (k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
							(k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1) ) { // C <-> T
							ts += 1; // Transition
						} else {
							tv += 1; // Transversion
						}
					}
	
		 			if (aa1 != aa2) {
						non = 1; // Is a non-synonymous change
					}
	
		 			// decide for rateClass
		 			switch(ts) {
		 				case 0:
		 					switch(tv) {
		 						case 1: rateClass = 3; break;	// 0ts, 1tv
								case 2: rateClass = 9; break;	// 0ts, 2tv
								case 3: rateClass = 17; break;	// 0ts, 3tv
								default: break;
		 					}
		 					break;
		 				case 1:
		 					switch(tv) {
		 						case 0: rateClass = 1; break;	// 1ts, 0tv
								case 1: rateClass = 7; break;	// 1ts, 1tv
								case 2: rateClass = 15; break;	// 1ts, 2tv
								default: break;
		 					}
		 					break;
		 				case 2:
		 					switch(tv) {
	 							case 0: rateClass = 5; break;	// 2ts, 0tv
	 							case 1: rateClass = 13; break;	// 2ts, 1tv
	 							default: break;
		 					}
		 					break;
		 				case 3:
		 					rateClass = 11; break;	// 3ts, 0tv
		 				default: break;
		 			}
		 			
		 			if(non == 1) {
		 				rateClass += 1;
		 			}
	
	                System.out.print("\t" + rateClass);

            }
            System.out.println();

        }
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

	public String toXHTML() {
		StringBuffer buffer = new StringBuffer();

		buffer.append("<em>Empirical Codon Model</em> omega = ");
		buffer.append(getOmega());
		buffer.append(", kappa_ts = ");
		buffer.append(getKappaTs());
		buffer.append(", kappa_tv = ");
		buffer.append(getKappaTv());
		buffer.append(", multi_nt = ");
		buffer.append(getMultiNt());
		buffer.append(", initial matrix = " + rateMat.getDirName() + "/" + rateMat.getMatName());
		buffer.append(", initial freqs = " + rateMat.getDirName() + "/" + rateMat.getFreqName());
		return buffer.toString();
	}

	
	static String format1 = "%2.1e";
	static String format2 = "%2.4e";
	
	public String printQ() {
		double[][] myQ = getQ();
		if (myQ != null) {
			StringBuffer sb = new StringBuffer();
			for (int i = 0; i < myQ.length; i++) {
				for(int j = 0; j < myQ[0].length; j++) {
					sb.append(String.format(format1, myQ[i][j]));
					sb.append("\t");
				}
				sb.append("\n");
			}
			return sb.toString();
		} else {
			return "No Q ready.";
		}
    }
	
	public String printRelRates() {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < relativeRates.length; i++) {
			sb.append(String.format(format2, relativeRates[i]));
			sb.append("\t");
		}
		sb.append("\n");
		return sb.toString();
	}
}
