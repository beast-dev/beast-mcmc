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

import dr.evolution.util.Taxa;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.util.HeapSort;
import dr.xml.*;
import static dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Andrew Holbrook
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class HawkesLikelihood extends AbstractModelLikelihood implements Reportable,
        GradientWrtParameterProvider {

    private final static String REQUIRED_FLAGS_PROPERTY = "hph.required.flags";
    private final static String HAWKES_LIKELIHOOD = "hawkesLikelihood";
    private final Double tolerance;


    public HawkesLikelihood(int hphDimension,
                            final Parameter tauXprec,
                            final Parameter sigmaXprec,
                            final Parameter tauTprec,
                            final Parameter omega,
                            final Parameter theta,
                            final Parameter mu0,
                            final MatrixParameterInterface locationsParameter,
                            final double[] times,
                            final Double tolerance,
                            boolean byIncrement) {

        super(HAWKES_LIKELIHOOD);

        this.hawkesModel = new HawkesModel(tauXprec, sigmaXprec, tauTprec, omega, theta, mu0, locationsParameter, times, byIncrement);

        this.hphDimension = hphDimension;
        this.locationCount = hawkesModel.getLocationCount();
        this.tolerance = tolerance;

        initialize(hphDimension, hawkesModel);
    }

    private class HawkesModel extends AbstractModel{

        final Parameter tauXprec;
        final Parameter sigmaXprec;
        final Parameter tauTprec;
        final Parameter omega;
        final Parameter theta;
        final Parameter mu0;
        final MatrixParameterInterface locationsParameter;
        final double[] times;
        final static String HAWKES_MODEL = "HawkesModel";
        final boolean byIncrement;

        public HawkesModel(final Parameter tauXprec,
                           final Parameter sigmaXprec,
                           final Parameter tauTprec,
                           final Parameter omega,
                           final Parameter theta,
                           final Parameter mu0,
                           final MatrixParameterInterface locationsParameter,
                           final double[] times,
                           boolean byIncrement) {
            super(HAWKES_MODEL);
            this.byIncrement = byIncrement;
            this.tauXprec = tauXprec;
            this.sigmaXprec = sigmaXprec;
            this.tauTprec = tauTprec;
            this.omega = omega;
            this.theta = theta;
            this.mu0 = mu0;
            this.locationsParameter = locationsParameter;
            this.times = times;

            checkDimensions();
            addVariable(tauXprec);
            addVariable(tauTprec);
            addVariable(sigmaXprec);
            addVariable(omega);
            addVariable(theta);
            addVariable(mu0);
            addVariable(locationsParameter);
        }

        private void checkDimensions() {
            if (times.length != getLocationCount()) {
                throw new RuntimeException("Times dimension doesn't match location count.");
            }

            if (getTotalDimension() != 6) {
                throw new RuntimeException("Parameter dimension is wrong.");
            }
        }

        public MatrixParameterInterface getLocationsParameter() {
            return locationsParameter;
        }

        public double[] getTimes() {
            return times;
        }

        private int getTotalDimension() {
            return sigmaXprec.getDimension() + tauXprec.getDimension() + tauTprec.getDimension() + omega.getDimension() + theta.getDimension() + mu0.getDimension();
        }

        public double[] getParameterValues() {
            double[] parameterValues = new double[]{
                    byIncrement ? sigmaXprec.getParameterValue(0) + tauXprec.getParameterValue(0) : sigmaXprec.getParameterValue(0),
                    tauXprec.getParameterValue(0),
                    tauTprec.getParameterValue(0),
                    byIncrement ? tauTprec.getParameterValue(0) + omega.getParameterValue(0) : omega.getParameterValue(0),
                    theta.getParameterValue(0),
                    mu0.getParameterValue(0)
            };
            return parameterValues;
        }

        public int getLocationCount() {
            return locationsParameter.getColumnDimension();
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {

        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        }

        @Override
        protected void storeState() {

        }

        @Override
        protected void restoreState() {

        }

        @Override
        protected void acceptState() {

        }
    }

    protected int initialize(
            final int hphDimension,
            final HawkesModel hawkesModel) {

        this.hphCore = getCore();

        System.err.println("Initializing with flags: " + flags);

        this.hphCore.initialize(hphDimension, locationCount, flags);
        this.hawkesModel = hawkesModel;
        int internalDimension = hphCore.getInternalDimension();
        setupLocationsParameter(hawkesModel.getLocationsParameter());


        hphCore.setParameters(hawkesModel.getParameterValues());

        updateAllLocations(hawkesModel.getLocationsParameter());

        addModel(hawkesModel);

        // make sure everything is calculated on first evaluation
//        makeDirty();

        return internalDimension;
    }


    @Override
    public String getReport() {
        return getId() + ": " + getLogLikelihood() + "\n" +
                GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, tolerance);
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return hawkesModel.getLocationsParameter();
    }

    @Override
    public int getDimension() {
        return hawkesModel.getLocationsParameter().getDimension();
    }
    @Override
    public double[] getGradientLogDensity() {
        // TODO Cache !!!
        if (gradient == null) {
            gradient = new double[hawkesModel.getLocationsParameter().getDimension()];
        }

        getLogLikelihood();

        hphCore.getGradient(gradient);

        return gradient; // TODO Do not expose internals
    }

    public enum ObservationType {
        POINT,
        UPPER_BOUND,
        LOWER_BOUND,
        MISSING

    }

    public MatrixParameterInterface getMatrixParameter() { return hawkesModel.getLocationsParameter(); }

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
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or hphPrecision need to be updated
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
            updateAllLocations(hawkesModel.getLocationsParameter());
            hphCore.setTimesData(hawkesModel.getTimes());
            hphCore.setParameters(hawkesModel.getParameterValues());
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
        final static String TIME_ATTRIBUTE_NAME = "timeTrait";
        final static String LOCATION_ATTRIBUTE_NAME = "locationTrait";
        final static String BY_INCREMENT = "byIncrement";
        final static String HPH_DIMENSION = "hphDimension";
        final static String SIGMA_PRECISON = "sigmaXprec";
        final static String TAU_X_PRECISION = "tauXprec";
        final static String TAU_T_PRECISION = "tauTprec";
        final static String OMEGA = "omega";
        final static String THETA = "theta";
        final static String MU = "mu0";
        final static String TOLERANCE = GRADIENT_CHECK_TOLERANCE;
        final static String JITTER = TreeTraitParserUtilities.JITTER;

        public String getParserName() {
            return HAWKES_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int hphDimension = xo.getIntegerAttribute(HPH_DIMENSION);

            MatrixParameterInterface locationsParameter = (MatrixParameterInterface) xo.getElementFirstChild(LOCATIONS);
            String timeTraitName = xo.getStringAttribute(TIME_ATTRIBUTE_NAME);
            String locationTraitName = xo.getStringAttribute(LOCATION_ATTRIBUTE_NAME);
            Taxa taxa = (Taxa) xo.getElementFirstChild(TIMES);
            double[] times = parseTimes(taxa, timeTraitName, locationsParameter, locationTraitName);


            if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
                TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
                List<Integer> missingIndices = TreeTraitParserUtilities.parseMissingIndices(locationsParameter, locationsParameter.getParameterValues());
                utilities.jitter(xo, hphDimension, missingIndices);
            }

            Parameter sigmaXprec = (Parameter) xo.getElementFirstChild(SIGMA_PRECISON);
            Parameter tauXprec = (Parameter) xo.getElementFirstChild(TAU_X_PRECISION);
            Parameter tauTprec = (Parameter) xo.getElementFirstChild(TAU_T_PRECISION);
            Parameter omega = (Parameter) xo.getElementFirstChild(OMEGA);
            Parameter theta = (Parameter) xo.getElementFirstChild(THETA);
            Parameter mu0 = (Parameter) xo.getElementFirstChild(MU);

            boolean byIncrement = xo.getAttribute(BY_INCREMENT, false);
            Double tolerance = xo.getAttribute(TOLERANCE, 1e-4);

            return new HawkesLikelihood(hphDimension, tauXprec, sigmaXprec, tauTprec, omega, theta, mu0, locationsParameter, times, tolerance, byIncrement);

        }

        private double[] parseTimes(Taxa taxa, String timeTraitName, MatrixParameterInterface locationsParameter, String locationName) {
            double[] times = new double[taxa.getTaxonCount()];
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                times[i] = Double.valueOf((String) taxa.getTaxon(i).getAttribute(timeTraitName));
                locationsParameter.getParameter(i).setId(taxa.getTaxonId(i));
                StringTokenizer st = new StringTokenizer((String) taxa.getTaxon(i).getAttribute(locationName));
                Parameter singleLocation = locationsParameter.getParameter(i);
                for (int j = 0; j < singleLocation.getDimension(); j++) {
                    singleLocation.setParameterValue(j, Double.valueOf(st.nextToken()));
                }
            }
            HeapSort.sort(times);
            final double tmp = times[0];
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                times[i] -= tmp;
            }

            return times;
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
                new ElementRule(TIMES, Taxa.class),
                AttributeRule.newStringRule(TIME_ATTRIBUTE_NAME),
                AttributeRule.newStringRule(LOCATION_ATTRIBUTE_NAME),
                AttributeRule.newBooleanRule(BY_INCREMENT, true),
                new ElementRule(SIGMA_PRECISON, Parameter.class),
                new ElementRule(TAU_X_PRECISION, Parameter.class),
                new ElementRule(TAU_T_PRECISION, Parameter.class),
                new ElementRule(OMEGA, Parameter.class),
                new ElementRule(THETA, Parameter.class),
                new ElementRule(MU, Parameter.class),
                AttributeRule.newDoubleRule(TOLERANCE, true),
                TreeTraitParserUtilities.jitterRules(true),
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
    private HawkesModel hawkesModel;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private long flags = 0;

    private double[] gradient;
}
