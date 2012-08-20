package test.dr.integration;

/**
 * GTR Parameter Estimation Tester.
 *
 * @author Walter Xie
 * @version 1.0
 * @since <pre>08/06/2009</pre>
 */
public class BEASTInputFile {

	public static final double ct = 1; // fixed	
	public static final double birthRate = 50;
	
	private double ac;
	private double ag;
	private double at;
	private double cg;
	private double gt;
    private String fileNamePrefix;    

    
	public double getAc() {
		return ac;
	}
	public void setAc(double ac) {
		this.ac = ac;
	}
	public double getAg() {
		return ag;
	}
	public void setAg(double ag) {
		this.ag = ag;
	}
	public double getAt() {
		return at;
	}
	public void setAt(double at) {
		this.at = at;
	}
	public double getCg() {
		return cg;
	}
	public void setCg(double cg) {
		this.cg = cg;
	}
	public double getGt() {
		return gt;
	}
	public void setGt(double gt) {
		this.gt = gt;
	}
	public String getFileNamePrefix() {
		return fileNamePrefix;
	}
	public void setFileNamePrefix (String fileNamePrefix) {
		this.fileNamePrefix = fileNamePrefix;
	}

 
}
