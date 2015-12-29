/*
 * GPUInfo.java
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
