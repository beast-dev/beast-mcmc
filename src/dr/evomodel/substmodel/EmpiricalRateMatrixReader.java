/*
 * EmpiricalRateMatrix.java
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

package dr.evomodel.substmodel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.logging.Logger;

import dr.evolution.datatype.DataType;

/**
 * Empirical rate matrix
 *
 * @author Stefan Zoller
 *
 */
public class EmpiricalRateMatrixReader {
    
	/**
     * constructor
     *
     * @param name						Name of matrix
     * @param dataType  				Data type as Codons.UNIVERSAL
     * @param dir						Directory that contains csv files with matrices
     * @param fName						Name of frequency csv file
     * @param mName						Name of initial matrix csv file
     */
    public EmpiricalRateMatrixReader(String name, DataType dataType, String dir, String fName, String mName) {
		this.name = name;
		this.dataType = dataType;
		this.dataDir = dir;
		this.matName = mName;
		this.freqName = fName;
		
		setupRates();
		setupFreqs();
	}
	
	public final String getName() { return name; }
	public final DataType getDataType() { return dataType; }
	
	public double[] getRates() { return rates; }
	public double[] getFrequencies() { return frequencies; }
	
	public String getFreqName() { return freqName; }
	public String getMatName() { return matName; }
	public String getDirName() { return dataDir; }
	
	public void setFrequencies(double[] f) {
	    this.frequencies = f;
	}
	
	public void setRates(double[] p) {
	    this.rates = p;
	}
	
	private void setupRates() {
		this.rates = new double[1830];
		readFile(this.rates, matName, 1830);
	}
	
	private void setupFreqs() {
		this.frequencies = new double[61];
		readFile(this.frequencies, freqName, 61);
	}
	
	// read csv data
	private void readFile(double[] r, String f, int nr) {
		File file = new File(dataDir, f);
	    
		try {
            BufferedReader bufRdr  = new BufferedReader(new FileReader(file));
            String line = null;
            int col = 0;

            while((line = bufRdr.readLine()) != null) {
        	    StringTokenizer st = new StringTokenizer(line,",");
        	    while (st.hasMoreTokens()) {
        		    r[col] = Double.valueOf(st.nextToken()).doubleValue();
        		    col++;
        	    }
            }
            bufRdr.close();
            if(col != nr) {
            	Logger.getLogger("dr.evomodel").severe("Matrix does not contain " + nr + " values but " + col);
            }
        } catch (FileNotFoundException e) {
        	Logger.getLogger("dr.evomodel").severe("Caught FileNotFoundException: " + e.getMessage());
        } catch (IOException e) {
        	Logger.getLogger("dr.evomodel").severe("Caught IOException: " + e.getMessage());
        } finally {
        }
	}
			
	protected double[] rates;
	protected double[] frequencies;
	
	private String name;
	protected String dataDir;
	protected String matName;
	protected String freqName;
	protected DataType dataType;
}
