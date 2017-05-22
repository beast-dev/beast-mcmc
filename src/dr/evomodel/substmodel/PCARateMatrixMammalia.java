/*
 * PCARateMatrixMammalia.java
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

import dr.evolution.datatype.Codons;
import java.util.*;
import java.io.*;

/**
 * Mammalia rate matrix for PCACodonModel
 *
 * @author Stefan Zoller
 *
 */
public class PCARateMatrixMammalia extends AbstractPCARateMatrix{
	
	// filenames for csv files including actual rate matrices
	public static final String FREQS_FILE = "freqs.csv";
	public static final String SCALES_FILE = "scales.csv";
	public static final String MEANS_FILE = "means.csv";
	public static final String FACTORS_FILE = "startfacs.csv";
	public static final String PCS_FILE = "pcs.csv";
    
    public PCARateMatrixMammalia(int nrOfPCs, String dir) {
		super("mammalia", Codons.UNIVERSAL, dir);
		
		// reading data files
		setupMeans();
		setupFreqs();
		setupScales();
		setupStartFacs();
		setupPCs(nrOfPCs);
	}
	
	private void setupFreqs() {
	    // 61 frequencies for the codons
	    double[] f = readSingleArray(FREQS_FILE, 61);
	    setFrequencies(f);
	}
	
	private void setupScales() {
	    // 1830 scale factors for the PCs
	    double[] s = readSingleArray(SCALES_FILE, 1830);
	    setScales(s);
	}
	
	private void setupMeans() {
	    // 1830 mean values for the PCs
	    double[] m = readSingleArray(MEANS_FILE, 1830);
	    setMeans(m);
	}
	
	private void setupStartFacs() {
	    // 1830 start factors for the PCs
	    double[] sf = readSingleArray(FACTORS_FILE, 1830);
	    setStartFacs(sf);
	}
	
	private void setupPCs(int nr) {
	    double[][] p = new double[nr][1830];
	    
	    File file = new File(dataDir, PCS_FILE);
	    
	    try {
            BufferedReader bufRdr  = new BufferedReader(new FileReader(file));
            
            String line = null;
            int row = 0;
            int col = 0;

            //read each line of text file
            while(((line = bufRdr.readLine()) != null) && (row < nr))
            {
            	col = 0;
            	StringTokenizer st = new StringTokenizer(line,",");
            	while (st.hasMoreTokens())
            	{
            		//get next token and store it in the array
            		p[row][col] = Double.valueOf(st.nextToken()).doubleValue();
            		col++;
            	}
            	row++;
            }
            bufRdr.close();
        } catch (FileNotFoundException e) {
            System.err.println("Caught FileNotFoundException: " + e.getMessage());
        } catch (IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
        } finally {
        }
	    setPCs(p);
    }
    
    private double[] readSingleArray(String filename, int dim) {
        File file = new File(dataDir, filename);
        double[] arr = new double[dim];
        
        try {
            BufferedReader bufRdr  = new BufferedReader(new FileReader(file));
            String line = null;
            int col = 0;

            while((line = bufRdr.readLine()) != null) {
        	    StringTokenizer st = new StringTokenizer(line,",");
        	    while (st.hasMoreTokens()) {
        		    arr[col] = Double.valueOf(st.nextToken()).doubleValue();
        		    col++;
        	    }
            }
            bufRdr.close();
        } catch (FileNotFoundException e) {
            System.err.println("Caught FileNotFoundException: " + e.getMessage());
        } catch (IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
        } finally {
        }
        return arr;
    }
}
