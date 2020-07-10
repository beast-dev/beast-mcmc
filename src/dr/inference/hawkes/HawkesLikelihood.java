/*
 * HawkesLikelihood.java
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

package dr.inference.hawkes;

import dr.evomodel.antigenic.MultidimensionalScalingLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class HawkesLikelihood extends AbstractModelLikelihood implements Reportable,
        GradientWrtParameterProvider {

    private final static String REQUIRED_FLAGS_PROPERTY = "hph.required.flags";

    @Override
    public String getReport() {
        return getId() + ": " + getLogLikelihood();
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return locationsParameter;
    }

    @Override
    public int getDimension() {
        return locationsParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        // TODO Cache !!!
        if (gradient == null) {
            gradient = new double[locationsParameter.getDimension()];
        }

        hphCore.getGradient(gradient);

        return gradient; // TODO Do not expose internals
    }

    public enum ObservationType {
        POINT,
        UPPER_BOUND,
        LOWER_BOUND,
        MISSING
    }

    private final static String HAWKES_LIKELIHOOD = "HawkesLikelihood";

    /**
     * Something
     * @param hphDimension
     * @param hphPrecision
     * @param locationsParameter
     */
    public HawkesLikelihood(
            int hphDimension,
            Parameter hphPrecision,
            MatrixParameterInterface locationsParameter,
            DataTable<double[]> dataTable,
            boolean isLeftTruncated,
            boolean reorderData) {

        super(HAWKES_LIKELIHOOD);

        this.hphDimension = hphDimension;

//        // construct a compact data table
//        String[] rowLabelsOriginal = dataTable.getRowLabels();
////        String[] columnLabels = dataTable.getRowLabels();
//
//        int rowCount = dataTable.getRowCount();
//        locationCount = rowCount;
//
//        int[] permute;
//        if (reorderData) {
//            permute = getPermutation(rowLabelsOriginal, locationsParameter);
//        } else {
//            permute = new int[locationCount];
//            for (int i = 0; i < locationCount; ++i) {
//                permute[i] = i; // identity
//            }
//        }
//
//        String[] rowLabels = new String[locationCount];
//
//        int observationCount = rowCount * rowCount;
////        double[] observations = new double[observationCount];
//        observations = new double[observationCount];
//        ObservationType[] observationTypes = new ObservationType[observationCount];
//
//        double[][] tmp = new double[rowCount][rowCount];
//
//        for (int i = 0; i < rowCount; i++) {
//            rowLabels[i] = rowLabelsOriginal[permute[i]];
//
//            double[] dataRow = dataTable.getRow(permute[i]);
//
//            for (int j = i + 1; j < rowCount; j++) {
//                tmp[i][j] = tmp[j][i] = dataRow[permute[j]];
//            }
//        }
//
//        int u = 0;
//        for (int i = 0; i < rowCount; i++) {
//            for (int j = 0; j < rowCount; j++) {
//                if (i == j) {
//                    observations[u] = 0.0;
//                    observationTypes[u] = ObservationType.POINT;
//                } else {
//                    observations[u] = tmp[i][j];
//                    if (Double.isNaN(observations[u])) {
//                        observationTypes[u] = ObservationType.MISSING;
//                    } else {
//                        observationTypes[u] = ObservationType.POINT;
//                    }
//                }
//                u++;
//            }
//        }
//
//        this.vectorDimension = initialize(hphDimension, hphPrecision, isLeftTruncated, locationsParameter,
//                rowLabels, observations, observationTypes);
    }

//    private class Data {
//        int observationCount;
//        double[] observations;
//        ObservationType[] observationTypes;
//
//        Data(int observationCount, double[] observations, ObservationType[] observationTypes) {
//            this.observationCount = observationCount;
//            this.observations = observations;
//            this.observationTypes = observationTypes;
//        }
//    }

    public MatrixParameterInterface getMatrixParameter() { return locationsParameter; }

//    private int[] getPermutation(String[] source, MatrixParameterInterface destination) {
//
//        if (source.length != destination.getColumnDimension()) {
//            throw new IllegalArgumentException("Dimension mismatch");
//        }
//
//        final int length = source.length;
//
//        Map<String,Integer> map = new HashMap<String, Integer>(destination.getColumnDimension());
//        for (int i = 0; i < length; ++i) {
//            map.put(source[i],i);
//        }
//
//        int[] permute = new int[length];
//        for (int i = 0; i < length; ++i) {
//            Integer p = map.get(destination.getParameter(i).getParameterName());
//            if (p == null) {
//                Logger.getLogger("dr.app.beagle").info("Missing label!!!");
//            } else {
//                permute[i] = p;
//            }
//        }
//
//        return permute;
//    }

    private HawkesCore getCore() {
        long computeMode = 0;
        String r = System.getProperty(REQUIRED_FLAGS_PROPERTY);
        if (r != null) {
            computeMode = Long.parseLong(r.trim());
        }

        HawkesCore core;
        //if (computeMode >= HawkesCore.USE_NATIVE_HPH) {
            System.err.println("Attempting to use a native HPH core with flag: " + computeMode + "; may the force be with you ....");
            core = new MassivelyParallelHPHImpl();
            flags = computeMode;
        //} else {
            //System.err.println("Computer mode found: " + computeMode + " vs. " + r);
            //core = new HawkesCoreImpl();

        //}
        return core;
    }

    public int getHphDimension() { return hphDimension; }

    public int getLocationCount() { return locationCount; }

    protected int initialize(
            final int hphDimension,
            final Parameter tauXprec,
            final Parameter sigmaXprec,
            final Parameter tauTprec,
            final Parameter omega,
            final Parameter theta,
            final Parameter mu0,
            final MatrixParameterInterface locationsParameter,
            final String[] locationLabels,
            final double[] timesData) {

        this.hphCore = getCore();

        System.err.println("Initializing with flags: " + flags);

        this.hphCore.initialize(hphDimension, locationCount, flags);
        this.locationLabels = locationLabels;

        this.locationsParameter = locationsParameter;
        int internalDimension = hphCore.getInternalDimension();
        setupLocationsParameter(this.locationsParameter);
        addVariable(locationsParameter);

        this.hphModelParameters[0] = sigmaXprec;
        addVariable(sigmaXprec);
        this.hphModelParameters[1] = tauXprec;
        addVariable(tauXprec);
        this.hphModelParameters[2] = tauTprec;
        addVariable(tauTprec);
        this.hphModelParameters[3] = omega;
        addVariable(omega);
        this.hphModelParameters[4] = theta;
        addVariable(theta);
        this.hphModelParameters[5] = mu0;
        addVariable(mu0);


        hphCore.setParameters(hphModelParameters.getParameterValues());

        updateAllLocations(locationsParameter);

        // make sure everything is calculated on first evaluation
        makeDirty();

        return internalDimension;
    }

    private void updateAllLocations(MatrixParameterInterface locationsParameter) {
        // TODO Can make more efficient (if necessary) using vectorDimension padding
        hphCore.updateLocation(-1, locationsParameter.getParameterValues());
    }

    private void setupLocationsParameter(MatrixParameterInterface locationsParameter) {
        final boolean exisitingParameter = locationsParameter.getColumnDimension() > 0;

        if (exisitingParameter){
            if (locationsParameter.getColumnDimension() != locationCount){
                throw new RuntimeException("locationsParameter column dimension ("+locationsParameter.getColumnDimension()+") is not equal to the locationCount ("+locationCount+")");
            }
            if (locationsParameter.getRowDimension() != hphDimension){
                throw new RuntimeException("locationsParameter row dimension ("+locationsParameter.getRowDimension()+") is not equal to the hphDimension ("+hphDimension+")");
            }
        } else {
            throw new IllegalArgumentException("Dimensions on matrix must be set");
        }

        for (int i = 0; i < locationLabels.length; i++) {
            if (locationsParameter.getParameter(i).getParameterName().compareTo(locationLabels[i]) != 0) {
                throw new RuntimeException("Mismatched trait parameter name (" + locationsParameter.getParameter(i).getParameterName() +
                        ") and data dimension name (" + locationLabels[i] + ")");
            }
        }

        for (int i = 0; i < locationsParameter.getColumnDimension(); ++i) {
            Parameter param = locationsParameter.getParameter(i);
            try {
                param.getBounds();
            } catch (NullPointerException exception) {
                param.addBounds(new Parameter.DefaultBounds(
                        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, param.getDimension()));
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or hphPrecision need to be updated

        if (variable == locationsParameter) {

            if (index == -1) {
                updateAllLocations(locationsParameter);
            } else {

                int locationIndex = index / hphDimension;
                hphCore.updateLocation(locationIndex, locationsParameter.getColumnValues(locationIndex));
            }
        } else if (variable == hphPrecisionParameter) {
            hphCore.setParameters(hphPrecisionParameter.getParameterValues());
        }

        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        hphCore.storeState();
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
        hphCore.restoreState();
    }

    @Override
    protected void acceptState() {
        hphCore.acceptState();
        // do nothing
    }

    public void makeDirty() {
        likelihoodKnown = false;
        hphCore.makeDirty();
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = hphCore.calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        final static String FILE_NAME = "fileName";
        final static  String LOCATIONS = "locations";
        final static String HPH_DIMENSION = "hphDimension";
        final static String HPH_PRECISION = "hphPrecision";
        final static String INCLUDE_TRUNCATION = "includeTruncation";
        final static String USE_OLD = "useOld";
        final static String FORCE_REORDER = "forceReorder";

        public String getParserName() {
            return HAWKES_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<double[]> distanceTable;
            try {
                distanceTable = DataTable.Double.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }

            if (distanceTable.getRowCount() != distanceTable.getColumnCount()) {
                throw new XMLParseException("Data table is not symmetrical.");
            }

            int hphDimension = xo.getIntegerAttribute(HPH_DIMENSION);

            MatrixParameterInterface locationsParameter = (MatrixParameterInterface) xo.getElementFirstChild(LOCATIONS);

            Parameter hphPrecision = (Parameter) xo.getElementFirstChild(HPH_PRECISION);

            boolean useOld = xo.getAttribute(USE_OLD, false);

            boolean includeTruncation = xo.getAttribute(INCLUDE_TRUNCATION, false);

            boolean forceReorder = xo.getAttribute(FORCE_REORDER, false);

//            if (useOld) {
//                System.err.println("USE OLD");
//                return new HawkesLikelihood(hphDimension, includeTruncation, hphPrecision, (MatrixParameter)locationsParameter, distanceTable);
//            } else {
                return new HawkesLikelihood(hphDimension, hphPrecision, locationsParameter,
                        distanceTable, includeTruncation, forceReorder);
//            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of pairwise distance given vectors of coordinates" +
                    "for points according to the multidimensional scaling scheme of XXX & Rafferty (to fill in).";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(HPH_DIMENSION, false, "The dimension of the space for HPH"),
                new ElementRule(LOCATIONS, MatrixParameterInterface.class),
                AttributeRule.newBooleanRule(USE_OLD, true),
                AttributeRule.newBooleanRule(INCLUDE_TRUNCATION, true),
                AttributeRule.newBooleanRule(FORCE_REORDER, true),
                new ElementRule(HPH_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return HawkesLikelihood.class;
        }
    };

    public double getHPHPrecision() {
        return hphPrecisionParameter.getParameterValue(0);
    }

    private final int hphDimension;
    private final int locationCount;

    private HawkesCore hphCore;

    private String[] locationLabels;

    private Parameter[] hphModelParameters;
    private MatrixParameterInterface locationsParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private long flags = 0;

    private double[] gradient;
}
