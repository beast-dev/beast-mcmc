/*
 * DistributionLikelihood.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.distribution;

import dr.util.Attribute;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * A class that returns the log likelihood of a set of data (statistics)
 * being distributed according to a distribution generated empirically from some data.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id:$
 */

public class EmpiricalDistributionLikelihood extends AbstractDistributionLikelihood {

    public static final String EMPIRICAL_DISTRIBUTION_LIKELIHOOD = "empricalDistributionLikelihood";

    private int from = -1;
    private int to = Integer.MAX_VALUE;

    public EmpiricalDistributionLikelihood(String fileName, boolean inverse) {
        super(null);

        readFile(fileName);

        this.inverse = inverse;    
    }

    protected void readFile(String fileName) {

        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));

            String line1 = reader.readLine();
            StringTokenizer st = new StringTokenizer(line1," ");
            values = new double[st.countTokens()];
            for(int i=0; i<values.length; i++)
                values[i] = Double.valueOf(st.nextToken());
            String line2 = reader.readLine();
            st = new StringTokenizer(line2," ");
            density = new double[st.countTokens()];
            for(int i=0; i<density.length; i++)
                density[i] = Double.valueOf(st.nextToken());

            reader.close();

        } catch (FileNotFoundException e) {
            System.err.println("File not found: "+fileName);
            System.exit(-1);
        } catch (IOException e) {
            System.err.println("IO exception reading: "+fileName);
            System.exit(-1);
        }

    }


    public void setRange(int from, int to) {
        this.from = from;
        this.to = to;
    }
  
    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        double logL = 0.0;

        for (Attribute<double[]> data : dataList) {
            for (int j = Math.max(0, from); j < Math.min(data.getAttributeValue().length, to); j++) {

                double value = data.getAttributeValue()[j];
                logL += logPDF(value);
            }
        }
        return logL;
    }

    protected double logPDF(double value) {
        return 0.0;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    protected  double[] values;
    protected  double[] density;

    protected boolean inverse;
}

