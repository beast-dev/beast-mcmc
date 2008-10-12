package dr.evomodel.newtreelikelihood;

/**
 * @author Marc A. Suchard
 */
public class GPUInfo {

	GPUInfo(String name, long memorySize) {
		this.name = name;
		this.memorySize = memorySize;
	}

	String name;
	long memorySize;
	
}
