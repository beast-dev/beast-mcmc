package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.HiddenCodons;
import dr.evoxml.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedYangModel extends YangCodonModel {

	public static final String MARKOV_MODULATED_YANG_MODEL = "markovModulatedYangCodonModel";
	public static final String HIDDEN_COUNT = "hiddenCount";
	public static final String SWITCHING_RATES = "switchingRates";

	private static final byte RATE = 5;


	public MarkovModulatedYangModel(
	                         HiddenCodons codonDataType,
	                         Parameter omegaParameter,
	                         Parameter kappaParameter,
	                         Parameter switchingRates,
	                         FrequencyModel freqModel) {

		super(codonDataType, omegaParameter, kappaParameter, freqModel);

		this.hiddenClassCount = codonDataType.getHiddenClassCount();
		this.switchingRates = switchingRates;

//		super(YANG_CODON_MODEL, codonDataType, freqModel);

//		this.codonDataType = codonDataType;
//		this.geneticCode = codonDataType.getGeneticCode();

	}

//	public void setupRelativeRates() {
//		double kappa = getKappa();
//		double[] omega = omegaParameter.getParameterValues();
//		double[] rates = switchingRates.getParameterValues();
//
//		int stateCount = this.stateCount/hiddenClassCount;
//
//		int index = 0;
//		for(int i=0; i<stateCount; i++) {
//			for(int j=i+1; j<stateCount; j++) {
//				for(int h=0; h<hiddenClassCount; h++) {
//					int d = getIndex(h*stateCount+i, h*stateCount+j, this.stateCount);
//					switch (rateMap[index]) {
//						case 0: relativeRates[d] = 0.0; break;			    // codon changes in more than one codon position
//						case 1: relativeRates[d] = kappa; break;		    // synonymous transition
//						case 2: relativeRates[d] = 1.0; break;			    // synonymous transversion
//						case 3: relativeRates[d] = kappa * omega[h]; break; // non-synonymous transition
//						case 4: relativeRates[d] = omega[h]; break;		    // non-synonymous transversion
//					}
//				}
//				index++;
//			}
//		}
//		// Add the switching class rates
//		int rateIndex = 0;
//		for(int g=0; g<hiddenClassCount; g++) {
//			for(int h=g+1; h<hiddenClassCount; h++ ) {  // from g -> h
//				for(int i=0; i<stateCount; i++) {
//					int d = getIndex(g*stateCount+i, h*stateCount+i, this.stateCount);
//					relativeRates[d] = rates[rateIndex];
//				}
//				rateIndex++;
//			}
//		}
//	}

	public void setupRelativeRates() {
		for(int i=0; i<relativeRates.length; i++)
			relativeRates[i] = 1.0;
	}

//    public double[] getEigenValues() {
//        synchronized (this) {
//            if (updateMatrix) {
//                setupMatrix();
//            }
//        }
//	    System.err.println("length = "+Eval.length);
//	    System.err.println("Eigenvalues = "+ new Vector(Eval));
//        return Eval;
//    }

	// Mapping: Matrix[i][j] = Compressed vector[i*(S - 3/2) - i^2 / 2 + j - 1]
	private int getIndex(int i, int j, int S) {
		return (i*(2*S - 3) - i*i)/2 + j - 1;
	}

		protected void constructRateMap() {
			// Construct map for non-hidden states only
			hiddenClassCount = ((HiddenCodons)codonDataType).getHiddenClassCount();
			stateCount /= hiddenClassCount;
			super.constructRateMap();			
			stateCount *= hiddenClassCount;
		}

		public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return MARKOV_MODULATED_YANG_MODEL; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//			Codons codons = Codons.UNIVERSAL;
//			if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
//				String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
//				if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
//					codons = Codons.UNIVERSAL;
//				} else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
//					codons = Codons.VERTEBRATE_MT;
//				} else if (codeStr.equals(GeneticCode.YEAST.getName())) {
//					codons = Codons.YEAST;
//				} else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
//					codons = Codons.MOLD_PROTOZOAN_MT;
//				} else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
//					codons = Codons.INVERTEBRATE_MT;
//				} else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
//					codons = Codons.CILIATE;
//				} else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
//					codons = Codons.ECHINODERM_MT;
//				} else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
//					codons = Codons.EUPLOTID_NUC;
//				} else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
//					codons = Codons.BACTERIAL;
//				} else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
//					codons = Codons.ALT_YEAST;
//				} else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
//					codons = Codons.ASCIDIAN_MT;
//				} else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
//					codons = Codons.FLATWORM_MT;
//				} else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
//					codons = Codons.BLEPHARISMA_NUC;
//				} else if (codeStr.equals(GeneticCode.NO_STOPS.getName())) {
//					codons = Codons.NO_STOPS;
//				}
//			}                        ÷

			DataType dataType = DataTypeUtils.getDataType(xo);
			HiddenCodons codons = null;
			if (dataType instanceof HiddenCodons)
				codons = (HiddenCodons) dataType;
			else
				throw new XMLParseException("Must construct "+MARKOV_MODULATED_YANG_MODEL+" with hidden codons");

			Parameter omegaParam = (Parameter)xo.getElementFirstChild(OMEGA);
			Parameter kappaParam = (Parameter)xo.getElementFirstChild(KAPPA);
			Parameter switchingParam = (Parameter)xo.getElementFirstChild(SWITCHING_RATES);
			FrequencyModel freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);
			return new MarkovModulatedYangModel(codons, omegaParam, kappaParam, switchingParam, freqModel);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the a Markov-modulated Yang model of codon evolution.";
		}

		public Class getReturnType() { return MarkovModulatedYangModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
//			new ElementRule(DataType.DATA_TYPE,
//				new XMLSyntaxRule[] { new ElementRule(HiddenCodons.class) }),
			AttributeRule.newStringRule(DataType.DATA_TYPE),
			new ElementRule(OMEGA,
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(KAPPA,
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(SWITCHING_RATES,
				new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
			new ElementRule(FrequencyModel.class)
		};
	};

	private int hiddenClassCount;
	private Parameter switchingRates;

}
