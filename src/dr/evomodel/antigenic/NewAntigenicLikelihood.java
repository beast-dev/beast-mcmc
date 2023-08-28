/*
 * AntigenicLikelihood.java
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

package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.LogTricks;
import dr.math.distributions.NormalDistribution;
import dr.util.*;
import dr.xml.*;


import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
/*
    Both virus locations and serum locations are shifted by the parameter locationDrift.
    A location is increased by locationDrift x offset.
    Offset is set to 0 for the earliest virus and increasing with difference in date from earliest virus.
    This makes the raw virusLocations and serumLocations parameters not directly interpretable.
*/
public class AntigenicLikelihood extends AbstractModelLikelihood implements Citable {
    private static final boolean CHECK_INFINITE = false;
    private static final boolean USE_THRESHOLDS = true;
    private static final boolean USE_INTERVALS = true;

    public final static String ANTIGENIC_LIKELIHOOD = "antigenicLikelihood";

    // column indices in table
    private static final int VIRUS_ISOLATE = 0;
    private static final int VIRUS_STRAIN = 1;
    private static final int VIRUS_DATE = 2;
    private static final int SERUM_ISOLATE = 3;
    private static final int SERUM_STRAIN = 4;
    private static final int SERUM_DATE = 5;
    private static final int TITRE = 6;

    public enum MeasurementType {
        INTERVAL,
        POINT,
        THRESHOLD,
        MISSING
    }

    public AntigenicLikelihood(
            int mdsDimension,
            Parameter mdsPrecisionParameter,
            Parameter locationDriftParameter,
            Parameter virusDriftParameter,
            Parameter serumDriftParameter,
            MatrixParameter virusLocationsParameter,
            MatrixParameter serumLocationsParameter,
            CompoundParameter tipTraitsParameter,
            Parameter virusOffsetsParameter,
            Parameter serumOffsetsParameter,
            Parameter serumPotenciesParameter,
            Parameter serumBreadthsParameter,
            Parameter virusAviditiesParameter,
            DataTable<String[]> dataTable,
            boolean mergeSerumIsolates,
            double intervalWidth,
            double driftInitialLocations) {

        super(ANTIGENIC_LIKELIHOOD);

        this.intervalWidth = intervalWidth;
        boolean useIntervals = USE_INTERVALS && intervalWidth > 0.0;

        int thresholdCount = 0;


        double earliestDate = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dataTable.getRowCount(); i++) {

            String[] values = dataTable.getRow(i);

            String virusName = values[VIRUS_STRAIN];
            double virusDate = Double.parseDouble(values[VIRUS_DATE]);
            int virus = virusNames.indexOf(virusName);
            if (virus == -1) {
                virusNames.add(virusName);
                virusDates.add(virusDate);
                virus = virusNames.size() - 1;
            }

            String serumName = "";
            if (mergeSerumIsolates) {
                serumName = values[SERUM_STRAIN];
            } else {
                serumName = values[SERUM_ISOLATE];
            }
            double serumDate = Double.parseDouble(values[SERUM_DATE]);
            int serum = serumNames.indexOf(serumName);
            if (serum == -1) {
                serumNames.add(serumName);
                serumDates.add(serumDate);
                serum = serumNames.size() - 1;
            }

            boolean isThreshold = false;
            boolean isLowerThreshold = false;
            double rawTitre = Double.NaN;
            if (values[TITRE].length() > 0) {
                try {
                    rawTitre = Double.parseDouble(values[TITRE]);
                } catch (NumberFormatException nfe) {
                    // check if threshold below
                    if (values[TITRE].contains("<")) {
                        rawTitre = Double.parseDouble(values[TITRE].replace("<",""));
                        isThreshold = true;
                        isLowerThreshold = true;
                        thresholdCount++;
                    }
                    // check if threshold above
                    if (values[TITRE].contains(">")) {                	
                        rawTitre = Double.parseDouble(values[TITRE].replace(">",""));
                        isThreshold = true;
                        isLowerThreshold = false;
                        thresholdCount++;                    	
                        //throw new IllegalArgumentException("Error in measurement: unsupported greater than threshold at row " + (i+1));
                    }
                }
            }

            if (serumDate < earliestDate) {
                earliestDate = serumDate;
            }

            if (virusDate < earliestDate) {
                earliestDate = virusDate;
            }

            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(virus, serum, virusDate, serumDate, type, rawTitre, isLowerThreshold);

            if (USE_THRESHOLDS || !isThreshold) {
                measurements.add(measurement);
            }

        }

        double[] maxColumnTitres = new double[serumNames.size()];
        for (Measurement measurement : measurements) {
            double titre = measurement.log2Titre;
            if (Double.isNaN(titre)) {
                titre = measurement.log2Titre;
            }
            if (titre > maxColumnTitres[measurement.serum]) {
                maxColumnTitres[measurement.serum] = titre;
            }
        }

        this.mdsDimension = mdsDimension;

        this.mdsPrecisionParameter = mdsPrecisionParameter;
        addVariable(mdsPrecisionParameter);

        this.locationDriftParameter = locationDriftParameter;
        if (this.locationDriftParameter != null) {
            addVariable(locationDriftParameter);
        }

        this.virusDriftParameter = virusDriftParameter;
        if (this.virusDriftParameter != null) {
            addVariable(virusDriftParameter);
        }

        this.serumDriftParameter = serumDriftParameter;
        if (this.serumDriftParameter != null) {
            addVariable(serumDriftParameter);
        }

        this.virusLocationsParameter = virusLocationsParameter;
        if (this.virusLocationsParameter != null) {
            setupLocationsParameter(virusLocationsParameter, virusNames);
        }

        this.serumLocationsParameter = serumLocationsParameter;
        if (this.serumLocationsParameter != null) {
            setupLocationsParameter(serumLocationsParameter, serumNames);
        }

        this.tipTraitsParameter = tipTraitsParameter;
        if (tipTraitsParameter != null) {
            setupTipTraitsParameter(this.tipTraitsParameter, virusNames);
        }

        this.virusOffsetsParameter = virusOffsetsParameter;
        if (virusOffsetsParameter != null) {
            setupOffsetsParameter(virusOffsetsParameter, virusNames, virusDates, earliestDate);
        }

        this.serumOffsetsParameter = serumOffsetsParameter;
        if (serumOffsetsParameter != null) {
            setupOffsetsParameter(serumOffsetsParameter, serumNames, serumDates, earliestDate);
        }

        this.serumPotenciesParameter = setupSerumPotencies(serumPotenciesParameter, maxColumnTitres);
        this.serumBreadthsParameter = setupSerumBreadths(serumBreadthsParameter);
        this.virusAviditiesParameter = setupVirusAvidities(virusAviditiesParameter);

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicLikelihood:\n");
        sb.append("\t\t" + virusNames.size() + " viruses\n");
        sb.append("\t\t" + serumNames.size() + " sera\n");
        sb.append("\t\t" + measurements.size() + " assay measurements\n");
        if (USE_THRESHOLDS) {
            sb.append("\t\t" + thresholdCount + " thresholded measurements\n");
        }
        if (useIntervals) {
            sb.append("\n\t\tAssuming a log 2 measurement interval width of " + intervalWidth + "\n");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());

        virusLocationChanged = new boolean[this.virusLocationsParameter.getParameterCount()];
        serumLocationChanged = new boolean[this.serumLocationsParameter.getParameterCount()];
        virusEffectChanged = new boolean[virusNames.size()];
        serumEffectChanged = new boolean[serumNames.size()];
        logLikelihoods = new double[measurements.size()];
        storedLogLikelihoods = new double[measurements.size()];

        setupInitialLocations(driftInitialLocations);

        makeDirty();
    }

    private Parameter setupVirusAvidities(Parameter virusAviditiesParameter) {
        // If no row parameter is given, then we will only use the serum effects
        if (virusAviditiesParameter != null) {
            virusAviditiesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, Double.MIN_VALUE, 1));
            virusAviditiesParameter.setDimension(virusNames.size());
            addVariable(virusAviditiesParameter);
            String[] labelArray = new String[virusNames.size()];
            virusNames.toArray(labelArray);
            virusAviditiesParameter.setDimensionNames(labelArray);
            for (int i = 0; i < virusNames.size(); i++) {
                virusAviditiesParameter.setParameterValueQuietly(i, 0.0);
            }
        }
        return virusAviditiesParameter;
    }

    private Parameter setupSerumPotencies(Parameter serumPotenciesParameter, double[] maxColumnTitres) {
        // If no serum potencies parameter is given, make one to hold maximum values for scaling titres...
        if (serumPotenciesParameter == null) {
            serumPotenciesParameter = new Parameter.Default("serumPotencies");
        } else {
            serumPotenciesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            addVariable(serumPotenciesParameter);
        }

        serumPotenciesParameter.setDimension(serumNames.size());
        String[] labelArray = new String[serumNames.size()];
        serumNames.toArray(labelArray);
        serumPotenciesParameter.setDimensionNames(labelArray);
        for (int i = 0; i < maxColumnTitres.length; i++) {
            serumPotenciesParameter.setParameterValueQuietly(i, maxColumnTitres[i]);
        }

        return serumPotenciesParameter;
    }

    private Parameter setupSerumBreadths(Parameter serumBreadthsParameter) {
        // If no serum breadths parameter is given, then we will only use the serum potencies
        if (serumBreadthsParameter != null) {
            serumBreadthsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            serumBreadthsParameter.setDimension(serumNames.size());
            addVariable(serumBreadthsParameter);
            String[] labelArray = new String[serumNames.size()];
            serumNames.toArray(labelArray);
            serumBreadthsParameter.setDimensionNames(labelArray);
            for (int i = 0; i < serumNames.size(); i++) {
                serumBreadthsParameter.setParameterValueQuietly(i, 1.0);
            }
        }
        return serumBreadthsParameter;
    }

    protected void setupLocationsParameter(MatrixParameter locationsParameter, List<String> strains) {
        locationsParameter.setColumnDimension(mdsDimension);
        locationsParameter.setRowDimension(strains.size());
        for (int i = 0; i < strains.size(); i++) {
            locationsParameter.getParameter(i).setId(strains.get(i));
        }
        addVariable(locationsParameter);
    }

    private void setupOffsetsParameter(Parameter offsetsParameter, List<String> strainNames, List<Double> strainDates, double earliest) {
        offsetsParameter.setDimension(strainNames.size());
        String[] labelArray = new String[strainNames.size()];
        strainNames.toArray(labelArray);
        offsetsParameter.setDimensionNames(labelArray);
        for (int i = 0; i < strainNames.size(); i++) {
            Double offset = strainDates.get(i) - new Double(earliest);
            if (offset == null) {
                throw new IllegalArgumentException("Date missing for strain: " + strainNames.get(i));
            }
            offsetsParameter.setParameterValue(i, offset);
        }
        addVariable(offsetsParameter);
    }


    private void setupTipTraitsParameter(CompoundParameter tipTraitsParameter, List<String> strainNames) {
        tipIndices = new int[strainNames.size()];

        for (int i = 0; i < strainNames.size(); i++) {
            tipIndices[i] = -1;
        }

        for (int i = 0; i < tipTraitsParameter.getParameterCount(); i++) {
            Parameter tip = tipTraitsParameter.getParameter(i);
            String label = tip.getParameterName();
            int index = findStrain(label, strainNames);
            if (index != -1) {
                if (tipIndices[index] != -1) {
                    throw new IllegalArgumentException("Duplicated tip name: " + label);
                }

                tipIndices[index] = i;

                // rather than setting these here, we set them when the locations are set so the changes propagate
                // through to the diffusion model.
//                Parameter location = locationsParameter.getParameter(index);
//                for (int dim = 0; dim < mdsDimension; dim++) {
//                    tip.setParameterValue(dim, location.getParameterValue(dim));
//                }
            } else {
                // The tree may contain viruses not present in the assay data
                //       throw new IllegalArgumentException("Unmatched tip name in assay data: " + label);
            }
        }
        // we are only setting this parameter not listening to it:
//        addVariable(this.tipTraitsParameter);
    }

    private final int findStrain(String label, List<String> strainNames) {
        int index = 0;
        for (String strainName : strainNames) {
            if (label.startsWith(strainName)) {
                return index;
            }

            index ++;
        }
        return -1;
    }

    private void setupInitialLocations(double drift) {
        for (int i = 0; i < virusLocationsParameter.getParameterCount(); i++) {
            double offset = 0.0;
            if (virusOffsetsParameter != null) {
                offset = drift * virusOffsetsParameter.getParameterValue(i);
            }
            double r = MathUtils.nextGaussian() + offset;
            virusLocationsParameter.getParameter(i).setParameterValue(0, r);
            if (mdsDimension > 1) {
                for (int j = 1; j < mdsDimension; j++) {
                    r = MathUtils.nextGaussian();
                    virusLocationsParameter.getParameter(i).setParameterValue(j, r);
                }
            }
        }
        for (int i = 0; i < serumLocationsParameter.getParameterCount(); i++) {
            double offset = 0.0;
            if (serumOffsetsParameter != null) {
                offset = drift * serumOffsetsParameter.getParameterValue(i);
            }
            double r = MathUtils.nextGaussian() + offset;
            serumLocationsParameter.getParameter(i).setParameterValue(0, r);
            if (mdsDimension > 1) {
                for (int j = 1; j < mdsDimension; j++) {
                    r = MathUtils.nextGaussian();
                    serumLocationsParameter.getParameter(i).setParameterValue(j, r);
                }
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == virusLocationsParameter) {
            if (index != -1) {
                int loc = index / mdsDimension;
                virusLocationChanged[loc] = true;
                if (tipTraitsParameter != null && tipIndices[loc] != -1) {
                    Parameter location = virusLocationsParameter.getParameter(loc);
                    Parameter tip = tipTraitsParameter.getParameter(tipIndices[loc]);
                    int dim = index % mdsDimension;
                    tip.setParameterValue(dim, location.getParameterValue(dim));
                }
            } else {
                if (false) { // Just for debugging purposes.
//                    for (int idx = 0; idx < virusLocationsParameter.getDimension(); ++idx) {
//                        int loc = idx / mdsDimension;
//                        virusLocationChanged[loc] = true;
//                        if (tipTraitsParameter != null && tipIndices[loc] != -1) {
//                            Parameter location = virusLocationsParameter.getParameter(loc);
//                            Parameter tip = tipTraitsParameter.getParameter(tipIndices[loc]);
//                            int dim = idx % mdsDimension;
//                            tip.setParameterValue(dim, location.getParameterValue(dim));
//                        }
//                    }
                } else {
                    Arrays.fill(virusLocationChanged, true);
                    if (tipTraitsParameter != null) {
                        for (int pindex = 0; pindex < virusLocationsParameter.getParameterCount(); ++pindex) {
                            Parameter location = virusLocationsParameter.getParameter(pindex);
                            Parameter tip = tipTraitsParameter.getParameter(tipIndices[pindex]);
                            for (int i = 0; i < tip.getDimension(); ++i) {
                                tip.setParameterValueQuietly(i, location.getParameterValue(i));//
                            }
                        }
                        tipTraitsParameter.fireParameterChangedEvent();
                    }
                }
            }
        } else if (variable == serumLocationsParameter) {
            int loc = index / mdsDimension;
            serumLocationChanged[loc] = true;
        } else if (variable == mdsPrecisionParameter) {
            setLocationChangedFlags(true);
        } else if (variable == locationDriftParameter) {
            setLocationChangedFlags(true);
        } else if (variable == virusDriftParameter) {
                setLocationChangedFlags(true);
        } else if (variable == serumDriftParameter) {
                setLocationChangedFlags(true);
        } else if (variable == serumPotenciesParameter) {
            serumEffectChanged[index] = true;
        } else if (variable == serumBreadthsParameter) {
            serumEffectChanged[index] = true;
        } else if (variable == virusAviditiesParameter) {
            virusEffectChanged[index] = true;
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        System.arraycopy(logLikelihoods, 0, storedLogLikelihoods, 0, logLikelihoods.length);
    }

    @Override
    protected void restoreState() {
        double[] tmp = logLikelihoods;
        logLikelihoods = storedLogLikelihoods;
        storedLogLikelihoods = tmp;

        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    private double computeLogLikelihood() {

        double precision = mdsPrecisionParameter.getParameterValue(0);
        double sd = 1.0 / Math.sqrt(precision);

        logLikelihood = 0.0;
        int i = 0;

        for (Measurement measurement : measurements) {

            if (virusLocationChanged[measurement.virus] || serumLocationChanged[measurement.serum] || virusEffectChanged[measurement.virus] || serumEffectChanged[measurement.serum]) {

                double expectation = calculateBaseline(measurement.virus, measurement.serum) - computeDistance(measurement.virus, measurement.serum);

                switch (measurement.type) {
                    case INTERVAL: {
                        double minTitre = measurement.log2Titre;
                        double maxTitre = measurement.log2Titre + intervalWidth;
                        logLikelihoods[i] = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
                    } break;
                    case POINT: {
                        logLikelihoods[i] = computeMeasurementLikelihood(measurement.log2Titre, expectation, sd);
                    } break;
                    case THRESHOLD: {
                    	if(measurement.isLowerThreshold){
                    		logLikelihoods[i] = computeMeasurementThresholdLikelihood(measurement.log2Titre, expectation, sd);
                    	}
                    	else{
                    		logLikelihoods[i] = computeMeasurementUpperThresholdLikelihood(measurement.log2Titre, expectation, sd);                  		
                    	}
                    } break;
                    case MISSING:
                        break;
                }
            }
            logLikelihood += logLikelihoods[i];
            i++;
        }

        likelihoodKnown = true;

        setLocationChangedFlags(false);
        setSerumEffectChangedFlags(false);
        setVirusEffectChangedFlags(false);

        return logLikelihood;
    }

    private void setLocationChangedFlags(boolean flag) {
        for (int i = 0; i < virusLocationChanged.length; i++) {
            virusLocationChanged[i] = flag;
        }
        for (int i = 0; i < serumLocationChanged.length; i++) {
            serumLocationChanged[i] = flag;
        }
    }

    private void setSerumEffectChangedFlags(boolean flag) {
        for (int i = 0; i < serumEffectChanged.length; i++) {
            serumEffectChanged[i] = flag;
        }
    }

    private void setVirusEffectChangedFlags(boolean flag) {
        for (int i = 0; i < virusEffectChanged.length; i++) {
            virusEffectChanged[i] = flag;
        }
    }

    // offset virus and serum location when computing
    protected double computeDistance(int virus, int serum) {

        Parameter vLoc = virusLocationsParameter.getParameter(virus);
        Parameter sLoc = serumLocationsParameter.getParameter(serum);
        double sum = 0.0;

        // first dimension is shifted
        double vxOffset = 0.0;
        double sxOffset = 0.0;
        if (locationDriftParameter != null && virusOffsetsParameter != null && serumOffsetsParameter != null) {
            vxOffset = locationDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
            sxOffset = locationDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
        }
        if (virusDriftParameter != null && virusOffsetsParameter != null) {
            vxOffset = virusDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
        }
        if (serumDriftParameter != null && serumOffsetsParameter != null) {
            sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
        }

        double vxLoc = vLoc.getParameterValue(0) + vxOffset;
        double sxLoc = sLoc.getParameterValue(0) + sxOffset;

        double difference = vxLoc - sxLoc;
        sum += difference * difference;

        // other dimensions are not
        for (int i = 1; i < mdsDimension; i++) {
            difference = vLoc.getParameterValue(i) - sLoc.getParameterValue(i);
            sum += difference * difference;
        }

        double dist = Math.sqrt(sum);

        if (serumBreadthsParameter != null) {
            double serumBreadth = serumBreadthsParameter.getParameterValue(serum);
            dist /= serumBreadth;
        }

        return dist;
    }


    // Calculates the expected log2 titre when mapDistance = 0
    private double calculateBaseline(int virus, int serum) {
        double baseline = serumPotenciesParameter.getParameterValue(serum);
        if (virusAviditiesParameter != null) {
            baseline += virusAviditiesParameter.getParameterValue(virus);
        }
        return baseline;
    }

    private static double computeMeasurementLikelihood(double titre, double expectation, double sd) {

        double lnL = NormalDistribution.logPdf(titre, expectation, sd);

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite point measurement");
        }
        return lnL;
    }

    private static double computeMeasurementThresholdLikelihood(double titre, double expectation, double sd) {

        // real titre is somewhere between -infinity and measured 'titre'
        // want the lower tail of the normal CDF

        double lnL = NormalDistribution.cdf(titre, expectation, sd, true);          // returns logged CDF

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite threshold measurement");
        }
        return lnL;
    }


    private static double computeMeasurementUpperThresholdLikelihood(double titre, double expectation, double sd) {

        // real titre is somewhere between -infinity and measured 'titre'
        // want the lower tail of the normal CDF
    	double L = NormalDistribution.cdf(titre, expectation, sd, false);          // returns  CDF
    	double lnL = Math.log(1-L);  //get the upper tail probability, then log it

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite threshold measurement");
        }
        return lnL;
    }    
    
    
    private static double computeMeasurementIntervalLikelihood(double minTitre, double maxTitre, double expectation, double sd) {

        // real titre is somewhere between measured minTitre and maxTitre

        double cdf1 = NormalDistribution.cdf(maxTitre, expectation, sd, true);     // returns logged CDF
        double cdf2 = NormalDistribution.cdf(minTitre, expectation, sd, true);     // returns logged CDF
        double lnL = LogTricks.logDiff(cdf1, cdf2);

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            // this occurs when the interval is in the far tail of the distribution, cdf1 == cdf2
            // instead return logPDF of the point
            lnL = NormalDistribution.logPdf(minTitre, expectation, sd);
            if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
                throw new RuntimeException("infinite interval measurement");
            }
        }
        return lnL;
    }

    public void makeDirty() {
        likelihoodKnown = false;
        setLocationChangedFlags(true);
    }

    private class Measurement {
        private Measurement(final int virus, final int serum, final double virusDate, final double serumDate, final MeasurementType type, final double titre, final boolean isLowerThreshold) {
            this.virus = virus;
            this.serum = serum;
            this.virusDate = virusDate;
            this.serumDate = serumDate;
            this.type = type;
            this.titre = titre;
            this.log2Titre = Math.log(titre) / Math.log(2);
            this.isLowerThreshold = isLowerThreshold;
        }

        final int virus;
        final int serum;
        final double virusDate;
        final double serumDate;
        final MeasurementType type;
        final double titre;
        final double log2Titre;
        final boolean isLowerThreshold;

    };

    private final List<Measurement> measurements = new ArrayList<Measurement>();
    private final List<String> virusNames = new ArrayList<String>();
    private final List<String> serumNames = new ArrayList<String>();
    private final List<Double> virusDates = new ArrayList<Double>();
    private final List<Double> serumDates = new ArrayList<Double>();

    private final int mdsDimension;
    private final double intervalWidth;
    private final Parameter mdsPrecisionParameter;
    private final Parameter locationDriftParameter;
    private final Parameter virusDriftParameter;
    private final Parameter serumDriftParameter;

    private final MatrixParameter virusLocationsParameter;
    private final MatrixParameter serumLocationsParameter;

    private final Parameter virusOffsetsParameter;
    private final Parameter serumOffsetsParameter;

    private final CompoundParameter tipTraitsParameter;
    private int[] tipIndices;

    private final Parameter virusAviditiesParameter;
    private final Parameter serumPotenciesParameter;
    private final Parameter serumBreadthsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    private final boolean[] virusLocationChanged;
    private final boolean[] serumLocationChanged;
    private final boolean[] serumEffectChanged;
    private final boolean[] virusEffectChanged;
    private double[] logLikelihoods;
    private double[] storedLogLikelihoods;

// **************************************************************
// XMLObjectParser
// **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String FILE_NAME = "fileName";
        public final static String TIP_TRAIT = "tipTrait";
        public final static String VIRUS_LOCATIONS = "virusLocations";
        public final static String SERUM_LOCATIONS = "serumLocations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MERGE_SERUM_ISOLATES = "mergeSerumIsolates";
        public static final String DRIFT_INITIAL_LOCATIONS = "driftInitialLocations";
        public static final String INTERVAL_WIDTH = "intervalWidth";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String LOCATION_DRIFT = "locationDrift";
        public static final String VIRUS_DRIFT = "virusDrift";
        public static final String SERUM_DRIFT = "serumDrift";
        public static final String VIRUS_AVIDITIES = "virusAvidities";
        public static final String SERUM_POTENCIES = "serumPotencies";
        public static final String SERUM_BREADTHS = "serumBreadths";
        public final static String VIRUS_OFFSETS = "virusOffsets";
        public final static String SERUM_OFFSETS = "serumOffsets";

        public String getParserName() {
            return ANTIGENIC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<String[]> assayTable;
            try {
                assayTable = DataTable.Text.parse(new FileReader(fileName), true, false);
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }
            System.out.println("Loaded HI table file: " + fileName);

            boolean mergeSerumIsolates = xo.getAttribute(MERGE_SERUM_ISOLATES, false);

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);
            double intervalWidth = 0.0;
            if (xo.hasAttribute(INTERVAL_WIDTH)) {
                intervalWidth = xo.getDoubleAttribute(INTERVAL_WIDTH);
            }

            double driftInitialLocations = 0.0;
            if (xo.hasAttribute(DRIFT_INITIAL_LOCATIONS)) {
                driftInitialLocations = xo.getDoubleAttribute(DRIFT_INITIAL_LOCATIONS);
            }

            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter virusLocationsParameter = null;
            if (xo.hasChildNamed(VIRUS_LOCATIONS)) {
                virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
            }

            MatrixParameter serumLocationsParameter = null;
            if (xo.hasChildNamed(SERUM_LOCATIONS)) {
                serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);
            }

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter locationDrift = null;
            if (xo.hasChildNamed(LOCATION_DRIFT)) {
                locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
            }

            Parameter virusDrift = null;
            if (xo.hasChildNamed(VIRUS_DRIFT)) {
                virusDrift = (Parameter) xo.getElementFirstChild(VIRUS_DRIFT);
            }

            Parameter serumDrift = null;
            if (xo.hasChildNamed(SERUM_DRIFT)) {
                serumDrift = (Parameter) xo.getElementFirstChild(SERUM_DRIFT);
            }

            Parameter virusOffsetsParameter = null;
            if (xo.hasChildNamed(VIRUS_OFFSETS)) {
                virusOffsetsParameter = (Parameter) xo.getElementFirstChild(VIRUS_OFFSETS);
            }

            Parameter serumOffsetsParameter = null;
            if (xo.hasChildNamed(SERUM_OFFSETS)) {
                serumOffsetsParameter = (Parameter) xo.getElementFirstChild(SERUM_OFFSETS);
            }

            Parameter serumPotenciesParameter = null;
            if (xo.hasChildNamed(SERUM_POTENCIES)) {
                serumPotenciesParameter = (Parameter) xo.getElementFirstChild(SERUM_POTENCIES);
            }

            Parameter serumBreadthsParameter = null;
            if (xo.hasChildNamed(SERUM_BREADTHS)) {
                serumBreadthsParameter = (Parameter) xo.getElementFirstChild(SERUM_BREADTHS);
            }

            Parameter virusAviditiesParameter = null;
            if (xo.hasChildNamed(VIRUS_AVIDITIES)) {
                virusAviditiesParameter = (Parameter) xo.getElementFirstChild(VIRUS_AVIDITIES);
            }

            AntigenicLikelihood AGL = new AntigenicLikelihood(
                    mdsDimension,
                    mdsPrecision,
                    locationDrift,
                    virusDrift,
                    serumDrift,
                    virusLocationsParameter,
                    serumLocationsParameter,
                    tipTraitParameter,
                    virusOffsetsParameter,
                    serumOffsetsParameter,
                    serumPotenciesParameter,
                    serumBreadthsParameter,
                    virusAviditiesParameter,
                    assayTable,
                    mergeSerumIsolates,
                    intervalWidth,
                    driftInitialLocations);

            Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Utils.getCitationString(AGL));

            return AGL;
        }

//************************************************************************
// AbstractXMLObjectParser implementation
//************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of immunological assay data such as Hemagglutinin inhibition (HI) given vectors of coordinates" +
                    "for viruses and sera/antisera in some multidimensional 'antigenic' space.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                AttributeRule.newBooleanRule(MERGE_SERUM_ISOLATES, true, "Should multiple serum isolates from the same strain have their locations merged (defaults to false)"),
                AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
                AttributeRule.newDoubleRule(DRIFT_INITIAL_LOCATIONS, true, "The degree to drift initial virus and serum locations, defaults to 0.0"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "Optional parameter of tip locations from the tree", true),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "Parameter of locations of all virus"),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "Parameter of locations of all sera"),
                new ElementRule(VIRUS_OFFSETS, Parameter.class, "Optional parameter for virus dates to be stored", true),
                new ElementRule(SERUM_OFFSETS, Parameter.class, "Optional parameter for serum dates to be stored", true),
                new ElementRule(SERUM_POTENCIES, Parameter.class, "Optional parameter for serum potencies", true),
                new ElementRule(SERUM_BREADTHS, Parameter.class, "Optional parameter for serum breadths", true),
                new ElementRule(VIRUS_AVIDITIES, Parameter.class, "Optional parameter for virus avidities", true),
                new ElementRule(MDS_PRECISION, Parameter.class, "Parameter for precision of MDS embedding"),
                new ElementRule(LOCATION_DRIFT, Parameter.class, "Optional parameter for drifting locations with time", true),
                new ElementRule(VIRUS_DRIFT, Parameter.class, "Optional parameter for drifting only virus locations, overrides locationDrift", true),
                new ElementRule(SERUM_DRIFT, Parameter.class, "Optional parameter for drifting only serum locations, overrides locationDrift", true)
        };

        public Class getReturnType() {
            return AntigenicLikelihood.class;
        }
    };

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian Antigenic Cartography framework";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(CommonCitations.BEDFORD_2015_INTEGRATING);
    }

    public static void main(String[] args) {
        double[] titres = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0};

        System.out.println("titre\tpoint\tinterval(tail)\tinterval(cdf)\tthreshold");
        for (double titre : titres) {
            double point = AntigenicLikelihood.computeMeasurementLikelihood(titre, 0.0, 1.0);
            double interval = AntigenicLikelihood.computeMeasurementIntervalLikelihood(titre + 1.0, titre, 0.0, 1.0);
            double threshold = AntigenicLikelihood.computeMeasurementThresholdLikelihood(titre, 0.0, 1.0);

            System.out.println(titre + "\t" + point + "\t" + interval + "\t" + threshold);
        }
    }
}
