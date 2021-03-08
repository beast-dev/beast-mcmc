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
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.util.HeapSort;
import dr.xml.*;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import static dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;

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
                            final HawkesRateProvider rateProvider,
                            final MatrixParameterInterface locationsParameter,
                            final double[] times,
                            final Double tolerance,
                            boolean byIncrement) {

        super(HAWKES_LIKELIHOOD);

        this.hawkesModel = new HawkesModel(tauXprec, sigmaXprec, tauTprec, omega, theta, mu0, rateProvider, locationsParameter, times, byIncrement);

        this.hphDimension = hphDimension;
        this.locationCount = hawkesModel.getLocationCount();
        this.tolerance = tolerance;

        initialize(hphDimension, hawkesModel);
    }

    public class HawkesModel extends AbstractModel{

        final Parameter tauXprec;
        final Parameter sigmaXprec;
        final Parameter tauTprec;
        final Parameter omega;
        final Parameter theta;
        final Parameter mu0;
        final HawkesRateProvider rateProvider;
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
                           final HawkesRateProvider rateProvider,
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
            this.rateProvider = rateProvider;
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
            if (rateProvider instanceof Model) {
                addModel((Model) rateProvider);
            }
        }

        private void checkDimensions() {
            if (times.length != getLocationCount()) {
                throw new RuntimeException("Times dimension doesn't match location count.");
            }

            if (getTotalDimension() != 6) {
                throw new RuntimeException("Parameter dimension is wrong.");
            }
        }

        public HawkesRateProvider getRateProvider() {
            return rateProvider;
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
            likelihoodKnown = false;
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
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
        hawkesModel.getRateProvider().setRandomRates(hphCore);

        updateAllLocations(hawkesModel.getLocationsParameter());

        addModel(hawkesModel);

        // make sure everything is calculated on first evaluation
//        makeDirty();

        return internalDimension;
    }


    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        double loglik = getLogLikelihood();
        sb.append(getClass().getName()).append("(").append(loglik).append(")");

        return sb.toString();
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    public HawkesModel getHawkesModel() {
        return hawkesModel;
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

        hphCore.getLocationGradient(gradient);

        return gradient; // TODO Do not expose internals
    }

    public double[] getRandomRateGradient() {

        double[] gradient = new double[hawkesModel.getRateProvider().getParameter().getDimension()];

        getLogLikelihood();

        hphCore.getRandomRatesGradient(gradient);

        hawkesModel.getRateProvider().updateRateGradient(gradient);

        return gradient;
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
            hawkesModel.getRateProvider().setRandomRates(hphCore);
            logLikelihood = hphCore.calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    enum UnitConversion {

        KM_TO_DEGREE("kmToDegree") {
            @Override
            double convert(double input) {
                return Math.pow(110.5,2) * 6.0 * Math.PI / input;
            }
        };

        UnitConversion(String name) {
            this.name = name;
        }

        private final String name;

        abstract double convert(double input);

        public static UnitConversion factory(String match) {
            for (UnitConversion type : UnitConversion.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }

    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        final static String LOCATIONS = "locations";
        final static String TIMES = "times";
        final static String TIME_ATTRIBUTE_NAME = "timeTrait";
        final static String LOCATION_ATTRIBUTE_NAME = "locationTrait";
        final static String LOCATION_MEAN = "locationMean";
        final static String LOCATION_MEAN_NAME = "locationMeanTrait";
        final static String LOCATION_VARIANCE = "locationVariance";
        final static String LOCATION_VARIANCE_NAME = "locationVarianceTrait";
        final static String LOCATION_VARIANCE_CONVERSION = "conversion";
        final static String BY_INCREMENT = "byIncrement";
        final static String HPH_DIMENSION = "hphDimension";
        final static String SIGMA_PRECISON = "sigmaXprec";
        final static String TAU_X_PRECISION = "tauXprec";
        final static String TAU_T_PRECISION = "tauTprec";
        final static String OMEGA = "omega";
        final static String THETA = "theta";
        final static String MU = "mu0";
        final static String RANDOM_RATES = "randomRates";
        final static String ON_TREE = "onTree";
        final static String NOT_ON_TREE = "notOnTree";
        final static String TOLERANCE = GRADIENT_CHECK_TOLERANCE;
        final static String JITTER = TreeTraitParserUtilities.JITTER;
        final static String TIME_EFFECT = "timeEffect";
        final static String GLM_COEFFICIENTS = "glmCoefficients";
        final static String HAS_INTERCEPT = "hasIntercept";

        public String getParserName() {
            return HAWKES_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int hphDimension = xo.getIntegerAttribute(HPH_DIMENSION);

            XMLObject cxo = xo.getChild(LOCATIONS);
            MatrixParameterInterface locationsParameter = (MatrixParameterInterface) cxo.getChild(MatrixParameterInterface.class);

            String timeTraitName = xo.getStringAttribute(TIME_ATTRIBUTE_NAME);
            String locationTraitName = xo.getStringAttribute(LOCATION_ATTRIBUTE_NAME);
            Taxa taxa = (Taxa) xo.getElementFirstChild(TIMES);
            double[] times = parseTimes(taxa, timeTraitName, locationsParameter, locationTraitName);
            MatrixParameterInterface locationMeans = null;
            MatrixParameterInterface locationPrecision = null;
            if (cxo.hasChildNamed(LOCATION_MEAN) && cxo.hasChildNamed(LOCATION_VARIANCE)) {
                locationMeans = (MatrixParameterInterface) cxo.getElementFirstChild(LOCATION_MEAN);
                String locationMeanTraitName = cxo.getChild(LOCATION_MEAN).getStringAttribute(LOCATION_MEAN_NAME);
                locationPrecision = (MatrixParameterInterface) cxo.getElementFirstChild(LOCATION_VARIANCE);
                String locationVarianceTraitName = cxo.getChild(LOCATION_VARIANCE).getStringAttribute(LOCATION_VARIANCE_NAME);
                parseTimes(taxa, timeTraitName, locationMeans, locationMeanTraitName);

                FastMatrixParameter tmpLocationVariance = new FastMatrixParameter("tmpLocationVariance", 1, locationPrecision.getColumnDimension(), 0.0, false);
                parseTimes(taxa, timeTraitName, tmpLocationVariance, locationVarianceTraitName);

                UnitConversion conversion = UnitConversion.factory((String) cxo.getChild(LOCATION_VARIANCE).getAttribute(LOCATION_VARIANCE_CONVERSION));

                for (int row = 0; row < locationPrecision.getRowDimension(); row++) {
                    for (int col = 0; col < locationPrecision.getColumnDimension(); col++) {
                        locationPrecision.setParameterValue(row, col, conversion.convert(tmpLocationVariance.getParameterValue(0, col)));
                    }
                }
            }

            if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
                TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
                List<Integer> missingIndices = TreeTraitParserUtilities.parseMissingIndices(locationsParameter, locationsParameter.getParameterValues());
                utilities.jitter(xo, hphDimension, missingIndices, locationsParameter.getDimension());
            }

            Parameter sigmaXprec = (Parameter) xo.getElementFirstChild(SIGMA_PRECISON);
            Parameter tauXprec = (Parameter) xo.getElementFirstChild(TAU_X_PRECISION);
            Parameter tauTprec = (Parameter) xo.getElementFirstChild(TAU_T_PRECISION);
            Parameter omega = (Parameter) xo.getElementFirstChild(OMEGA);
            Parameter theta = (Parameter) xo.getElementFirstChild(THETA);
            Parameter mu0 = (Parameter) xo.getElementFirstChild(MU);

            HawkesRateProvider rateProvider;
            if (xo.hasChildNamed(RANDOM_RATES)) {
                XMLObject dxo = xo.getChild(RANDOM_RATES);
                Parameter onTreeRates = (Parameter) dxo.getChild(ON_TREE).getChild(Parameter.class);
                Parameter notOnTreeRates = dxo.getChild(NOT_ON_TREE) == null ? null : (Parameter) dxo.getChild(NOT_ON_TREE).getChild(Parameter.class);
                Parameter rates = notOnTreeRates == null ? onTreeRates : new CompoundParameter("HawkesRates", new Parameter[]{onTreeRates, notOnTreeRates});
                Taxa taxaOnTree = (Taxa) dxo.getChild(Taxa.class);
                TreeModel tree = (TreeModel) dxo.getChild(TreeModel.class);
                boolean[] onTree = new boolean[rates.getDimension()];
                int[] indices = new int[rates.getDimension()];
                double mostRecentTipHeight = mapRandomRates(taxa, timeTraitName, taxaOnTree, tree, onTree, indices);
                if (dxo.hasChildNamed(GLM_COEFFICIENTS)) {
                    boolean hasIntercept = dxo.getAttribute(HAS_INTERCEPT, false);
                    boolean timeEffect = dxo.getAttribute(TIME_EFFECT, false);
                    Parameter coefficients = (Parameter) dxo.getChild(GLM_COEFFICIENTS).getChild(Parameter.class);
                    double[] orderedHeights = new double[times.length];
                    for (int i = 0; i < times.length; i++) {
                        orderedHeights[i] = mostRecentTipHeight + times[i];
                    }
                    rateProvider = new HawkesRateProvider.GLM(rates, coefficients, tree, indices, orderedHeights, onTree, timeEffect, hasIntercept);
                } else {
                    rateProvider = new HawkesRateProvider.Default(rates, indices, onTree);
                }
            } else {
                rateProvider = new HawkesRateProvider.None();
            }

            boolean byIncrement = xo.getAttribute(BY_INCREMENT, false);
            Double tolerance = xo.getAttribute(TOLERANCE, 1e-4);

            return new HawkesLikelihood(hphDimension, tauXprec, sigmaXprec, tauTprec, omega, theta, mu0, rateProvider, locationsParameter, times, tolerance, byIncrement);

        }

        private double mapRandomRates(Taxa timeTaxa, String timeTraitName, Taxa taxaOnTree, TreeModel tree, boolean[] onTree, int[] indices) {
            double[] times = new double[timeTaxa.getTaxonCount()];
            for (int i = 0; i < timeTaxa.getTaxonCount(); i++) {
                times[i] = Double.valueOf((String) timeTaxa.getTaxon(i).getAttribute(timeTraitName));
            }
            int[] timeIndices = new int[times.length];
            HeapSort.sort(times, timeIndices);

            Map<Taxon, Integer> onTreeParameterIndexMap = new HashMap<>();
            int parameterIndex = 0;
            Double mostRecentTip = Double.POSITIVE_INFINITY;
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                Taxon currentTip = tree.getTaxon(i);
                if (taxaOnTree.contains(currentTip)) {
                    onTreeParameterIndexMap.put(currentTip, parameterIndex);
                    if (mostRecentTip > tree.getNodeHeight(tree.getNode(i))) {
                        mostRecentTip = tree.getNodeHeight(tree.getNode(i));
                    }
                    parameterIndex++;
                }
            }

            int offTreeParameterIndex = taxaOnTree.getTaxonCount();
            for (int i = 0; i < timeTaxa.getTaxonCount(); i++) {
                Taxon currentTaxon = timeTaxa.getTaxon(timeIndices[i]);
                if (taxaOnTree.contains(currentTaxon)) {
                    onTree[i] = true;
                    indices[i] = onTreeParameterIndexMap.get(currentTaxon);
                } else {
                    onTree[i] = false;
                    indices[i] = offTreeParameterIndex;
                    offTreeParameterIndex++;
                }
            }
            return mostRecentTip;
        }

        private double[] parseTimes(Taxa taxa, String timeTraitName, MatrixParameterInterface locationsParameter, String locationName) {
            double[] times = new double[taxa.getTaxonCount()];
            double[] orderedTimes = new double[taxa.getTaxonCount()];
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                times[i] = Double.valueOf((String) taxa.getTaxon(i).getAttribute(timeTraitName));
            }
            int[] timeIndices = new int[times.length];
            HeapSort.sort(times, timeIndices);
            final double tmp = times[timeIndices[0]];
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                orderedTimes[i] = times[timeIndices[i]] - tmp;
                locationsParameter.getParameter(i).setId(taxa.getTaxonId(timeIndices[i]));
                StringTokenizer st = new StringTokenizer((String) taxa.getTaxon(timeIndices[i]).getAttribute(locationName));
                Parameter singleLocation = locationsParameter.getParameter(i);
                for (int j = 0; j < singleLocation.getDimension(); j++) {
                    singleLocation.setParameterValue(j, Double.valueOf(st.nextToken()));
                }
            }

            return orderedTimes;
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
                new ElementRule(LOCATIONS,
                        new XMLSyntaxRule[]{
                                new ElementRule(MatrixParameterInterface.class),
                                new ElementRule("Optional location prior related trait construction",
                                        new XMLSyntaxRule[] {
                                                new AndRule( new XMLSyntaxRule[]{
                                                        new ElementRule(LOCATION_MEAN, MatrixParameterInterface.class),
                                                        AttributeRule.newStringRule(LOCATION_MEAN_NAME)}),
                                                new AndRule( new XMLSyntaxRule[]{
                                                        new ElementRule(LOCATION_VARIANCE, MatrixParameterInterface.class),
                                                        AttributeRule.newStringRule(LOCATION_VARIANCE_NAME)}),
                                                AttributeRule.newStringRule(LOCATION_VARIANCE_CONVERSION)
                                        },
                                        true
                                )
                        }),
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
                new ElementRule(RANDOM_RATES,
                        new XMLSyntaxRule[]{
                                new ElementRule(ON_TREE, Parameter.class, "The random rates that are on the tree."),
                                new ElementRule(NOT_ON_TREE, Parameter.class, "The random rates that are not on the tree."),
                                new ElementRule(Taxa.class, "Taxa that contains taxons with known locations on the tree"),
                                new ElementRule(TreeModel.class, "TreeModel that contains the trait parameter."),
                                new ElementRule(GLM_COEFFICIENTS, Parameter.class, "GLM coefficient parameter.", true),
                                AttributeRule.newStringRule(TIME_EFFECT, true),
                                AttributeRule.newStringRule(HAS_INTERCEPT, true)
                        }, true),
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
