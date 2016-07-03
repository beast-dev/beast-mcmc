/*
 * AbstractPCARateMatrix.java
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
