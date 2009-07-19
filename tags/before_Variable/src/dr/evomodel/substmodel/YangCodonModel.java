/*
 * YangCodonModel.java
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

import dr.evolution.datatype.*;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;
import jebl.evolution.sequences.NucleotideState;

/**
 * Yang model of codon evolution
 *
 * @version $Id: YangCodonModel.java,v 1.21 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class YangCodonModel extends AbstractCodonModel
{
	public static final String YANG_CODON_MODEL = "yangCodonModel";
	public static final String OMEGA = "omega";
	public static final String KAPPA = "kappa";

	/** kappa */
	protected Parameter kappaParameter;

	/** omega */
	protected Parameter omegaParameter;

	protected byte[] rateMap;

	/**
	 * Constructor
	 */
	public YangCodonModel(Codons codonDataType,
							Parameter omegaParameter,
							Parameter kappaParameter,
							FrequencyModel freqModel)
	{
		super(YANG_CODON_MODEL, codonDataType, freqModel);

		this.omegaParameter = omegaParameter;
		addParameter(omegaParameter);
		omegaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
				omegaParameter.getDimension()));

		this.kappaParameter = kappaParameter;
		addParameter(kappaParameter);
		kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
				kappaParameter.getDimension()));

		constructRateMap();

        addStatistic(synonymousRateStatistic);
	}

	/**
	 * set kappa
	 */
	public void setKappa(double kappa) {
		kappaParameter.setParameterValue(0, kappa);
		updateMatrix = true;
	}

	/**
	 * @return kappa
	 */
	public double getKappa() { return kappaParameter.getParameterValue(0); }

	/**
	 * set dN/dS
	 */
	public void setOmega(double omega) {
		omegaParameter.setParameterValue(0, omega);
		updateMatrix = true;
	}


	/**
	 * @return dN/dS
	 */
	public double getOmega() { return omegaParameter.getParameterValue(0); }

    public double getSynonymousRate() {
        double k = getKappa();
        double o = getOmega();
        return ((31.0 * k) + 36.0) / ((31.0 * k) + 36.0 + (138.0 * o) + (58.0 * o * k));
    }

    public double getNonSynonymousRate() {
        return 0;
    }

	/**
	 * setup substitution matrix
	 */
	public void setupRelativeRates()
	{
		double kappa = getKappa();
		double omega = getOmega();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
				case 0: relativeRates[i] = 0.0; break;			// codon changes in more than one codon position
				case 1: relativeRates[i] = kappa; break;		// synonymous transition
				case 2: relativeRates[i] = 1.0; break;			// synonymous transversion
				case 3: relativeRates[i] = kappa * omega; break;// non-synonymous transition
				case 4: relativeRates[i] = omega; break;		// non-synonymous transversion
			}
		}


	}

	/**
	 * Construct a map of the rate classes in the rate matrix using the current
	 * genetic code. Classes:
	 *		0: codon changes in more than one codon position (or stop codons)
	 *		1: synonymous transition
	 *		2: synonymous transversion
	 *		3: non-synonymous transition
	 *		4: non-synonymous transversion
	 */
	protected void constructRateMap()
	{
		int u, v, i1, j1, k1, i2, j2, k2;
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

				codon = codonDataType.getTripletStates(v);
				i2 = codon[0];
				j2 = codon[1];
				k2 = codon[2];

				cs2 = codonDataType.getState(i2, j2, k2);
				aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));

				rateClass = -1;
				if (i1 != i2) {
					if ( (i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
						 (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1) ) { // C <-> T
						rateClass = 1; // Transition at position 1
					} else {
						rateClass = 2; // Transversion at position 1
					}
				}
				if (j1 != j2) {
					if (rateClass == -1) {
						if ( (j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
							 (j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1) ) { // C <-> T
							rateClass = 1; // Transition
						} else {
							rateClass = 2; // Transversion
						}
					} else
						rateClass = 0; // Codon changes at more than one position
				}
				if (k1 != k2) {
					if (rateClass == -1) {
						if ( (k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
							 (k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1) ) { // C <-> T
							rateClass = 1; // Transition
						} else {
							rateClass = 2; // Transversion
						}
					} else
						rateClass = 0; // Codon changes at more than one position
				}

	 			if (rateClass != 0) {
					if (aa1 != aa2) {
						rateClass += 2; // Is a non-synonymous change
					}
				}

				rateMap[i] = rateClass;
				i++;
			}

		}
	}

    public void printRateMap()
    {
        int u, v, i1, j1, k1, i2, j2, k2;
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

                codon = codonDataType.getTripletStates(v);
                i2 = codon[0];
                j2 = codon[1];
                k2 = codon[2];

                cs2 = codonDataType.getState(i2, j2, k2);
                aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));

                rateClass = -1;
                if (i1 != i2) {
                    if ( (i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
                         (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1) ) { // C <-> T
                        rateClass = 1; // Transition at position 1
                    } else {
                        rateClass = 2; // Transversion at position 1
                    }
                }
                if (j1 != j2) {
                    if (rateClass == -1) {
                        if ( (j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
                             (j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1) ) { // C <-> T
                            rateClass = 1; // Transition
                        } else {
                            rateClass = 2; // Transversion
                        }
                    } else
                        rateClass = 0; // Codon changes at more than one position
                }
                if (k1 != k2) {
                    if (rateClass == -1) {
                        if ( (k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
                             (k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1) ) { // C <-> T
                            rateClass = 1; // Transition
                        } else {
                            rateClass = 2; // Transversion
                        }
                    } else
                        rateClass = 0; // Codon changes at more than one position
                }

                 if (rateClass != 0) {
                    if (aa1 != aa2) {
                        rateClass += 2; // Is a non-synonymous change
                    }
                }

                System.out.print("\t" + rateClass);

            }
            System.out.println();

        }
    }

	// *****************************************************************
	// Interface ModelComponent
	// *****************************************************************

	/**
	 * Reads a the model from an XMLObject.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return YANG_CODON_MODEL; }

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
			Parameter kappaParam = (Parameter)xo.getElementFirstChild(KAPPA);
			FrequencyModel freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);
			YangCodonModel codonModel = new YangCodonModel(codons, omegaParam, kappaParam, freqModel);

//            codonModel.printRateMap();

            return codonModel;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the Yang model of codon evolution.";
		}

		public Class getReturnType() { return YangCodonModel.class; }

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
			new ElementRule(OMEGA,
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(KAPPA,
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(FrequencyModel.class)
		};
	};


    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

	public String toXHTML() {
		StringBuffer buffer = new StringBuffer();

		buffer.append("<em>Yang Codon Model</em> kappa = ");
		buffer.append(getKappa());
		buffer.append(", omega = ");
		buffer.append(getOmega());

		return buffer.toString();
	}

    private Statistic synonymousRateStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "synonymousRate";
        }

        public int getDimension() { return 1; }

        public double getStatisticValue(int dim) {
            return getSynonymousRate();
        }

    };

   /* private Statistic nonsynonymousRateStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "nonSynonymousRate";
        }

        public int getDimension() { return 1; }

        public double getStatisticValue(int dim) {
            return getNonSynonymousRate();
        }

    };*/
}
