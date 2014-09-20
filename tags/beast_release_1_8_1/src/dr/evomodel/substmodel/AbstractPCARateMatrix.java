package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;

/**
 * Abstract baseclass for PCA rate matrix
 *
 * @author Stefan Zoller
 *
 */
public abstract class AbstractPCARateMatrix {
    
    /**
     * constructor
     *
     * @param name      Name of matrix
     * @param dataType	Data type as Codons.UNIVERSAL
     * @param dir		Directory which includes the rate matrix csv files
     */
    public AbstractPCARateMatrix(String name, DataType dataType, String dir) {
		this.name = name;
		this.dataType = dataType;
		this.dataDir = dir;
	}
	
	public static final String getName() { return name; }
	public final DataType getDataType() { return dataType; }
	
	public double[] getPCAt(int i) { return pcs[i]; }
	public double[] getFrequencies() { return frequencies; }
	public double[] getMeans() { return means; }
	public double[] getScales() { return scales; }
	public double[] getStartFacs() { return startFacs; }
	
	public void setFrequencies(double[] f) {
	    this.frequencies = f;
	}
	
	public void setMeans(double[] m) {
	    this.means = m;
	}
	
	public void setScales(double[] s) {
	    this.scales = s;
	}
	
	public void setPCs(double[][] p) {
	    this.pcs = p;
	}
	
	public void setStartFacs(double[] sf) {
	    this.startFacs = sf;
	}
			
	protected double[][] pcs;
	protected double[] means;
	protected double[] frequencies;
	protected double[] scales;
	protected double[] startFacs;
	
	private static String name;
	protected static String dataDir;
	protected DataType dataType;
}
