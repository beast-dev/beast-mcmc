package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.Codons;
import dr.app.beagle.evomodel.substmodel.BaseSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.DefaultEigenSystem;
import dr.app.beagle.evomodel.substmodel.EigenSystem;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.parsers.PCACodonModelParser;
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
     * @param AbstractPCARateMatrix		Rate matrix with PCs, means, scale factors
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
