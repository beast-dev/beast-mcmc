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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.*;

/**
 * @author Andrew Holbrook
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class HawkesLikelihood extends AbstractModelLikelihood implements Reportable,
        GradientWrtParameterProvider {

    private final static String REQUIRED_FLAGS_PROPERTY = "hph.required.flags";
    private final static String HAWKES_LIKELIHOOD = "hawkesLikelihood";


    public HawkesLikelihood(
            int hphDimension,
            HawkesParameters hawkesParameters) {

        super(HAWKES_LIKELIHOOD);

        this.hawkesParameters = hawkesParameters;

        this.hphDimension = hphDimension;
        this.locationCount = hawkesParameters.getLocationCount();

        initialize(hphDimension, hawkesParameters);
    }

    public static class HawkesParameters {

        final Parameter tauXprec;
        final Parameter sigmaXprec;
        final Parameter tauTprec;
        final Parameter omega;
        final Parameter theta;
        final Parameter mu0;
        final MatrixParameterInterface locationsParameter;
        final Parameter times;
        final CompoundParameter allParameters;

        public HawkesParameters(final Parameter tauXprec,
                                final Parameter sigmaXprec,
                                final Parameter tauTprec,
                                final Parameter omega,
                                final Parameter theta,
                                final Parameter mu0,
                                final MatrixParameterInterface locationsParameter,
                                final Parameter times) {
            this.tauXprec = tauXprec;
            this.sigmaXprec = sigmaXprec;
            this.tauTprec = tauTprec;
            this.omega = omega;
            this.theta = theta;
            this.mu0 = mu0;
            this.locationsParameter = locationsParameter;
            this.times = times;
            this.allParameters = new CompoundParameter("hphModelParameter", new Parameter[]{sigmaXprec, tauXprec, tauTprec, omega, theta, mu0});

            checkDimensions();
        }

        private void checkDimensions() {
            if (times.getDimension() != getLocationCount()) {
                throw new RuntimeException("Times dimension doesn't match location count.");
            }
        }

        public CompoundParameter getCompoundParameter() {
            return allParameters;
        }

        public MatrixParameterInterface getLocationsParameter() {
            return locationsParameter;
        }

        public double[] getTimes() {
            return times.getParameterValues();
        }

        public double[] getParameterValues() {
            return allParameters.getParameterValues();
        }

        public int getLocationCount() {
            return locationsParameter.getColumnDimension();
        }
    }

    protected int initialize(
            final int hphDimension,
            final HawkesParameters hawkesParameters) {

        this.hphCore = getCore();

        System.err.println("Initializing with flags: " + flags);

        this.hphCore.initialize(hphDimension, locationCount, flags);
        this.hawkesParameters = hawkesParameters;
        int internalDimension = hphCore.getInternalDimension();
        setupLocationsParameter(hawkesParameters.getLocationsParameter());
        addVariable(hawkesParameters.getCompoundParameter());


        hphCore.setParameters(hawkesParameters.getParameterValues());

        updateAllLocations(hawkesParameters.getLocationsParameter());

        // make sure everything is calculated on first evaluation
//        makeDirty();

        return internalDimension;
    }


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
        return hawkesParameters.getLocationsParameter();
    }

    @Override
    public int getDimension() {
        return hawkesParameters.getLocationsParameter().getDimension();
    }
    @Override
    public double[] getGradientLogDensity() {
        // TODO Cache !!!
        if (gradient == null) {
            gradient = new double[hawkesParameters.getLocationsParameter().getDimension()];
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

    public MatrixParameterInterface getMatrixParameter() { return hawkesParameters.getLocationsParameter(); }

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

//        for (int i = 0; i < locationLabels.length; i++) {
//            if (locationsParameter.getParameter(i).getParameterName().compareTo(locationLabels[i]) != 0) {
//                throw new RuntimeException("Mismatched trait parameter name (" + locationsParameter.getParameter(i).getParameterName() +
//                        ") and data dimension name (" + locationLabels[i] + ")");
//            }
//        }

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

//        if (variable == locationsParameter) {
//
//            if (index == -1) {
//                updateAllLocations(locationsParameter);
//            } else {
//
//                int locationIndex = index / hphDimension;
//                hphCore.updateLocation(locationIndex, locationsParameter.getColumnValues(locationIndex));
//            }
//        }
//        else if (variable == hphPrecisionParameter) {
//            hphCore.setParameters(hphPrecisionParameter.getParameterValues());
//        }

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
            updateAllLocations(hawkesParameters.getLocationsParameter());
            hphCore.setTimesData(hawkesParameters.getTimes());
            hphCore.setParameters(hawkesParameters.getParameterValues());
            logLikelihood = hphCore.calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        final static String LOCATIONS = "locations";
        final static String TIMES = "times";
        final static String HPH_DIMENSION = "hphDimension";
        final static String SIGMA_PRECISON = "sigmaXprec";
        final static String TAU_X_PRECISION = "tauXprec";
        final static String TAU_T_PRECISION = "tauTprec";
        final static String OMEGA = "omega";
        final static String THETA = "theta";
        final static String MU = "mu0";

        public String getParserName() {
            return HAWKES_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int hphDimension = xo.getIntegerAttribute(HPH_DIMENSION);

            MatrixParameterInterface locationsParameter = (MatrixParameterInterface) xo.getElementFirstChild(LOCATIONS);
            Parameter times = (Parameter) xo.getElementFirstChild(TIMES);

            Parameter sigmaXprec = (Parameter) xo.getElementFirstChild(SIGMA_PRECISON);
            Parameter tauXprec = (Parameter) xo.getElementFirstChild(TAU_X_PRECISION);
            Parameter tauTprec = (Parameter) xo.getElementFirstChild(TAU_T_PRECISION);
            Parameter omega = (Parameter) xo.getElementFirstChild(OMEGA);
            Parameter theta = (Parameter) xo.getElementFirstChild(THETA);
            Parameter mu0 = (Parameter) xo.getElementFirstChild(MU);

            HawkesParameters hawkesParameters = new HawkesParameters(tauXprec, sigmaXprec, tauTprec, omega, theta, mu0, locationsParameter, times);

            return new HawkesLikelihood(hphDimension, hawkesParameters);

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
                AttributeRule.newIntegerRule(HPH_DIMENSION, false, "The dimension of the space for HPH"),
                new ElementRule(LOCATIONS, MatrixParameterInterface.class),
                new ElementRule(TIMES, Parameter.class),
                new ElementRule(SIGMA_PRECISON, Parameter.class),
                new ElementRule(TAU_X_PRECISION, Parameter.class),
                new ElementRule(TAU_T_PRECISION, Parameter.class),
                new ElementRule(OMEGA, Parameter.class),
                new ElementRule(THETA, Parameter.class),
                new ElementRule(MU, Parameter.class)
        };

        public Class getReturnType() {
            return HawkesLikelihood.class;
        }
    };

//    public Parameter getHPHPrecision() {
//        return hphPrecisionParameter;
//    }

    private final int hphDimension;
    private final int locationCount;

    private HawkesCore hphCore;
    private HawkesParameters hawkesParameters;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private long flags = 0;

    private double[] gradient;
}
