package dr.evomodel.substmodel;

import dr.evolution.datatype.Codons;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.AbstractCodonModel;
import dr.evomodelxml.substmodel.PCACodonModelParser;
import dr.evomodel.substmodel.AbstractPCARateMatrix;
import dr.inference.model.Parameter;

/**
 * PCA model of codon evolution
 * 
 * @author Stefan Zoller
 */
public class PCACodonModel extends AbstractCodonModel {
	// principal components, means, scale factors
	protected AbstractPCARateMatrix rateMatrix;

	protected byte[] rateMap;
	
	private Parameter pcFactors;

    /**
     * constructor
     *
     * @param codonDataType				Data type as Codons.UNIVERSAL
     * @param pcaType           		Rate matrix with PCs, means, scale factors
     * @param pcaDimensionParameter		Scalars for PCs
     * @param freqModel					Frequency model
     */
	public PCACodonModel(Codons codonDataType,
						    AbstractPCARateMatrix pcaType,
						    Parameter pcaDimensionParameter,
						    FrequencyModel freqModel)
	{
		super(PCACodonModelParser.PCA_CODON_MODEL, codonDataType, freqModel);
		
		this.rateMatrix = pcaType;
		
		// initialize scalars for principal components
		this.pcFactors = pcaDimensionParameter;
		double[] startFacs = pcaType.getStartFacs();
		double facSum = 0.0;
		for(int i=0; i<pcFactors.getDimension(); i++) {
			facSum += startFacs[i];
		}
		for(int i=0; i<pcFactors.getDimension(); i++) {
			pcFactors.setParameterValueQuietly(i, startFacs[i]/facSum);
		}
		addVariable(pcFactors);
		pcFactors.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
				pcFactors.getDimension()));
	}

    
	// setup substitution matrix
    public void setupRelativeRates() {
        double[] m = rateMatrix.getMeans();
        double[] sc = rateMatrix.getScales();
        for (int i = 0; i < rateCount; i++) {
            relativeRates[i] = m[i];
        }
        for(int j = 0; j < pcFactors.getDimension(); j++) {
        	double[] pc = rateMatrix.getPCAt(j);
        	double f = pcFactors.getParameterValue(j);
        	for (int i = 0; i < rateCount; i++) {
                relativeRates[i] += f*pc[i]*sc[i];
            }
        }
        for (int i = 0; i < rateCount; i++) {
        	if(relativeRates[i]<0) {
        		relativeRates[i] = 0.0;
        	}
        }
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
