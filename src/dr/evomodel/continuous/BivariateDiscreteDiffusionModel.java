/*
 * BivariateDiscreteDiffusionModel.java
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

package dr.evomodel.continuous;

import dr.inference.model.Parameter;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 14, 2007
 * Time: 7:01:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class BivariateDiscreteDiffusionModel extends MultivariateDiffusionModel {

    public static final String DISCRETE_DIFFUSION_PROCESS = "multivariateDiscreteDiffusionModel";
    public static final String GRID_X_DIMENSION = "xGridDimension";
    public static final String GRID_Y_DIMENSION = "yGridDimension";
    //	public static final String EIGEN_FILE = "eigensystemFileName";
    public static final String EVEC_NAME = "eigenvectorsFileName";
    public static final String EVAL_NAME = "eigenvaluesFileName";

    /**
     * Construct a discrete diffusion model.
     */
    public BivariateDiscreteDiffusionModel(Parameter graphRate,
                                           int xDim, int yDim,
                                           double[] eVal, double[][] eVec) {
        super();
//		super(diffusionPrecisionMatrixParameter);
        this.graphRate = graphRate;
        this.xDim = xDim;
        this.yDim = yDim;
        this.totalDim = xDim * yDim;
        this.eVal = eVal;
        this.eVec = eVec;
//		probabilityCache = new HashMap<Integer,Double>();

        addVariable(graphRate);

        System.err.println("TEST00 = " + getCTMCProbability(0, 0, 0.0));
        System.err.println("TEST01 = " + getCTMCProbability(0, 1, 0.0));
    }

    private int getIndex(int x, int y) {
        return x * yDim + y;
    }

    private int[] getXY(int index) {
        int[] xy = new int[2];
        xy[0] = index / yDim;
        xy[1] = index - xy[0] * yDim;
        return xy;
    }

    private double getProbability(int I, int J, double time) {
        double probability = 0;
        for (int k = 0; k < totalDim; k++) {
            probability += eVec[I][k] * Math.exp(time * eVal[k]) * eVec[J][k];
        }
        return probability;

    }

    public void handleParameterChangedEvent(Parameter parameter, int index) {
        // Clear cached probabilities
//		probabilityCache.clear();
    }

    private double getCTMCProbability(int I, int J, double time) {
        double probability = 0;
        for (int k = 0; k < totalDim; k++) {
            probability += eVec[I][k] * Math.exp(time * eVal[k]) * eVec[J][k];
        }
        return probability;

    }

    /**
     * @return the log likelihood of going from start to stop in the given time
     */
    public double getLogLikelihood(double[] start, double[] stop, double time) {
//		System.err.println("xy0 = "+start[0]+":"+start[1]);
//		System.err.println("xy1 = "+ stop[0]+":"+ stop[1]);


        int startIndex = getIndex((int) start[0], (int) start[1]);
        int stopIndex = getIndex((int) stop[0], (int) stop[1]);

//		System.err.println("i0 = "+startIndex);
//		System.err.println("i1 = "+stopIndex);
//		System.err.println("");

        return Math.log(getCTMCProbability(startIndex, stopIndex, time));

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DISCRETE_DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int xDim = xo.getIntegerAttribute(GRID_X_DIMENSION);
            int yDim = xo.getIntegerAttribute(GRID_Y_DIMENSION);

            String evecFileName = xo.getStringAttribute(EVEC_NAME);
            String evalFileName = xo.getStringAttribute(EVAL_NAME);

            File evecFile;
            File evalFile;

            try {
                File file = new File(evecFileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                evecFile = new File(parent, name);
                new FileReader(evecFile);
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + evecFileName + "' can not be opened for " + getParserName() + " element.");
            }

            try {
                File file = new File(evalFileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                evalFile = new File(parent, name);
                new FileReader(evalFile);
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + evalFileName + "' can not be opened for " + getParserName() + " element.");
            }

            double[] eval = null;
            double[][] evec = null;

            try {
                eval = TopographicalMap.readEigenvalues(evalFile.getAbsolutePath());
                evec = TopographicalMap.readEigenvectors(evecFile.getAbsolutePath());
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                assert false;
            }

//			XMLObject cxo = (XMLObject) xo.getChild(DIFFUSION_CONSTANT);
//			MatrixParameter diffusionParam = (MatrixParameter) cxo.getChild(MatrixParameter.class);
//			MatrixParameter diffusionParam = null;
            Parameter graphRate = (Parameter) xo.getChild(Parameter.class);

//			if (diffusionParam.getRowDimension()>1 && diffusionParam.getColumnDimension()>1)
//				throw new XMLParseException("The bivariate discrete diffusion model currently uses a single rate.");

            final int totalDim = xDim * yDim;
            if (totalDim != evec.length || totalDim != eval.length)
                throw new XMLParseException("Number of eigenvalues and eigenvectors must match the map grid dimensions");

            return new BivariateDiscreteDiffusionModel(graphRate, xDim, yDim, eval, evec);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Describes a multivariate discrete diffusion process.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(EVEC_NAME),
                AttributeRule.newStringRule(EVAL_NAME),
                AttributeRule.newIntegerRule(GRID_X_DIMENSION),
                AttributeRule.newIntegerRule(GRID_Y_DIMENSION),
                new ElementRule(Parameter.class),
//				new ElementRule(DIFFUSION_CONSTANT,
//						new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
        };

        public Class getReturnType() {
            return BivariateDiscreteDiffusionModel.class;
        }

    };

    private final Parameter graphRate;
    private final int xDim;
    private final int yDim;
    private final int totalDim;
    HashMap<Integer, Double> probabilityCache;

    private final double[] eVal;
    private final double[][] eVec;

}
