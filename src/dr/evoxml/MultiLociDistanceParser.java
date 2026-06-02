/*
 * MultiLociDistanceParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evoxml;

import dr.xml.*;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * Calculates overall mean pairwise distances across multiple loci
 * 
 */
public class MultiLociDistanceParser extends AbstractXMLObjectParser {

    public static final String MULTI_LOCI_DISTANCE = "multiLociDistance";
    public String getParserName() {
        return MULTI_LOCI_DISTANCE; 
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int childNum = xo.getChildCount();
        ArrayList<DistanceMatrix> distMatList = new ArrayList<DistanceMatrix>();
        for(int i = 0; i < childNum; i++){
            distMatList.add((DistanceMatrix)xo.getChild(i));
        }
        double[] meanDists = new double[childNum];
        for(int i = 0; i < childNum; i++){
            meanDists[i] = distMatList.get(i).getMeanDistance();
        }

        double sum = 0.0;
        for(int i = 0; i < childNum; i++){
            sum = sum+meanDists[i];
        }
        double mean = sum/(double)childNum;

        printMeans(distMatList,mean);

        return mean;

    }

    public void printMeans(ArrayList<DistanceMatrix> distMat, double mean){
        System.out.println("Individual mean distances:");
        for(int i = 0; i < distMat.size(); i++){
            System.out.print(distMat.get(i).getId()+": ");
            System.out.println(distMat.get(i).getMeanDistance());
        }
        System.out.println("overallMean: "+mean);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(DistanceMatrix.class, 1, Integer.MAX_VALUE)
    };

    public String getParserDescription() {
        return "Constructs a distance matrix from a pattern list or alignment";
    }

    public Class getReturnType() { return Double.class; }
}
