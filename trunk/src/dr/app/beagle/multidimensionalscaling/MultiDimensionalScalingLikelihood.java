/*
 * MultidimensionalScalingLikelihood.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.multidimensionalscaling;

import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class MultiDimensionalScalingLikelihood extends AbstractModelLikelihood {


    public enum ObservationType {
        POINT,
        UPPER_BOUND,
        LOWER_BOUND,
        MISSING
    }

    public final static String MULTIDIMENSIONAL_SCALING_LIKELIHOOD = "multiDimensionalScalingLikelihood";

    /**
     * A simple constructor for a fully specified symmetrical data matrix
     * @param mdsDimension
     * @param mdsPrecision
     * @param locationsParameter
     * @param dataTable
     */
    public MultiDimensionalScalingLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            MatrixParameter locationsParameter,
            DataTable<double[]> dataTable) {

        super(MULTIDIMENSIONAL_SCALING_LIKELIHOOD);

        this.mdsDimension = mdsDimension;

        // construct a compact data table
        String[] rowLabels = dataTable.getRowLabels();
        String[] columnLabels = dataTable.getRowLabels();

        int rowCount = dataTable.getRowCount();
        locationCount = rowCount;

        int observationCount = rowCount * rowCount;
        double[] observations = new double[observationCount];
        ObservationType[] observationTypes = new ObservationType[observationCount];

        double[][] tmp = new double[rowCount][rowCount];

        for (int i = 0; i < rowCount; i++) {
            double[] dataRow = dataTable.getRow(i);

            for (int j = i + 1; j < rowCount; j++) {
                tmp[i][j] = tmp[j][i] = dataRow[j];
            }
        }

        int u = 0;
        for (int i = 0; i < rowCount; i++) {
            for (int j = 0; j < rowCount; j++) {
                observations[u] = (i == j ? 0 : tmp[i][j]);
                observationTypes[u] = ObservationType.POINT;
                u++;
            }

        }

        initialize(mdsDimension, mdsPrecision, locationsParameter, rowLabels, observations, observationTypes);
    }

    protected void initialize(
            final int mdsDimension,
            final Parameter mdsPrecision,
            final MatrixParameter locationsParameter,
            final String[] locationLabels,
            final double[] observations,
            final ObservationType[] observationTypes) {

        this.mdsCore = new MultiDimensionalScalingCoreImpl2();
        this.mdsCore.initialize(mdsDimension, locationCount);
        this.locationLabels = locationLabels;

        this.locationsParameter = locationsParameter;
        setupLocationsParameter(this.locationsParameter);
        addVariable(locationsParameter);

        this.mdsPrecisionParameter = mdsPrecision;
        addVariable(mdsPrecision);

        mdsCore.setParameters(mdsPrecisionParameter.getParameterValues());
        mdsCore.setPairwiseData(observations);
        for (int i = 0; i < locationCount; i++) {
            mdsCore.updateLocation(i, locationsParameter.getColumnValues(i));
        }

        // make sure everything is calculated on first evaluation
        makeDirty();
    }

    protected void setupLocationsParameter(MatrixParameter locationsParameter) {
        if (locationsParameter.getColumnDimension() > 0){
            if (locationsParameter.getColumnDimension() != locationCount){
                throw new RuntimeException("locationsParameter column dimension ("+locationsParameter.getColumnDimension()+") is not equal to the locationCount ("+locationCount+")");
            }
            if (locationsParameter.getRowDimension() != mdsDimension){
                throw new RuntimeException("locationsParameter row dimension ("+locationsParameter.getRowDimension()+") is not equal to the mdsDimension ("+mdsDimension+")");
            }
        } else{
            locationsParameter.setColumnDimension(mdsDimension);
            locationsParameter.setRowDimension(locationCount);
        }
        for (int i = 0; i < locationLabels.length; i++) {
            locationsParameter.getParameter(i).setId(locationLabels[i]);
        }

        for (int i = 0; i < locationsParameter.getParameterCount(); ++i) {
            Parameter param = locationsParameter.getParameter(i);
            try {
                if (param.getBounds() != null) {
                    // Do nothing
                }
            } catch (NullPointerException exception) {
                param.addBounds(new Parameter.DefaultBounds(
                        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, param.getDimension()));
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated

        if (variable == locationsParameter) {
            int locationIndex = index / mdsDimension;

            mdsCore.updateLocation(locationIndex, locationsParameter.getColumnValues(locationIndex));
        } else if (variable == mdsPrecisionParameter) {
            mdsCore.setParameters(mdsPrecisionParameter.getParameterValues());
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        mdsCore.storeState();
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
        mdsCore.restoreState();
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        likelihoodKnown = false;
        mdsCore.makeDirty();
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = mdsCore.calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";

        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MDS_PRECISION = "mdsPrecision";

        public String getParserName() {
            return MULTIDIMENSIONAL_SCALING_LIKELIHOOD;
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

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            return new MultiDimensionalScalingLikelihood(mdsDimension, mdsPrecision, locationsParameter, distanceTable);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of pairwise distance given vectors of coordinates" +
                    "for points according to the multidimensional scaling scheme of XXX & Rafferty (xxxx).";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return MultiDimensionalScalingLikelihood.class;
        }
    };

    private final int mdsDimension;
    private final int locationCount;

    private MultiDimensionalScalingCore mdsCore;

    private String[] locationLabels;

    private Parameter mdsPrecisionParameter;
    private MatrixParameter locationsParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;
}
