package dr.evomodel.newtreelikelihood;

import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */
public class GPUInfo {
	
//	GPUInfo() {
//		
//	}
	
//	GPUInfo(String name, long memorySize) {
//		numberDevices = 1;
//		this.name = new String[] {name};
//		this.memorySize = new long[] {memorySize};
//	}
//	
	GPUInfo(int numberDevices) {
		this.numberDevices = numberDevices;
		nameList = new String[numberDevices];
		memorySizeList = new int[numberDevices];
		clockSpeedList = new int[numberDevices];
	}
	
	public void addInfo(int deviceNumber, String name, int memorySize, int clockSpeed) {
		nameList[deviceNumber] = name;
		memorySizeList[deviceNumber] = memorySize;
		clockSpeedList[deviceNumber] = clockSpeed;
	}

	int numberDevices;
	String[] nameList;
	int[]	 memorySizeList;
	int[]	 clockSpeedList;

	
	public String getName(int deviceNumber) { return nameList[deviceNumber]; }
	public int getMemorySize(int deviceNumber) { return (int)Math.round((double)memorySizeList[deviceNumber] / (1024.0 * 1024.0)); }
	public double getClockSpeed(int deviceNumber) { return (double)clockSpeedList[deviceNumber] / 1000000.0; }
	
	public String toString() {
		
		StringBuffer bf = new StringBuffer();
		bf.append("Total number of GPU devices: "+numberDevices+"\n");
		for(int i=0; i<numberDevices; i++) {
			bf.append("\t#"+(i+1)+":\t"+getName(i)+"\n");
			bf.append("\t\tGlobal memory (MB) "+getMemorySize(i)+"\n");
			bf.append("\t\tClock speed (GHz) "+getClockSpeed(i)+"\n");
		}						
		return bf.toString();				
	}
	
}
