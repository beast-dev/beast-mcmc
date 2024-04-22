/*
 * TreeTraitParserUtilities.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.TaxonList;
import dr.evomodel.continuous.StandardizeTraits;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.getTipTraitName;

/**
 * @author Marc A. Suchard
 */

public class TreeTraitParserUtilities {

    public static final String TRAIT_PARAMETER = "traitParameter";
    public static final String TRAIT_NAME = "traitName";
    public static final String MISSING = "missingIndicator";
    public static final String RANDOM_SAMPLE = "randomSample";
    public static final String DEFAULT_TRAIT_NAME = "trait";

    public static final String RANDOMIZE = "randomize";
    public static final String RANDOMIZE_LOWER = "lower";
    public static final String RANDOMIZE_UPPER = "upper";

    public static final String ALLOW_IDENTICAL = "allowIdentical";
    public static final String JITTER = "jitter";
    public static final String WINDOW = "window";
    public static final String DUPLICATES = "duplicatesOnly";
    public static final String STANDARDIZE = "standardize";
    public static final String TARGET_SD = "targetSd";
    public static final String SAMPLE_MISSING_TRAITS = "sampleMissingTraits";

    public static final String LATENT_FROM = "latentFrom";
    public static final String LATENT_TO = "latentTo";

    public void randomize(Parameter trait, double[] lower, double[] upper) {
        // Draws each dimension in each trait from U[lower, upper)
        for (int i = 0; i < trait.getDimension(); i++) {
            final int whichLower = i % lower.length;
            final int whichUpper = i % upper.length;
            final double newValue = MathUtils.uniform(lower[whichLower], upper[whichUpper]);
            trait.setParameterValue(i, newValue);
        }
    }

//    public void standardize(Parameter trait) {
//        for (int i = 0; i < trait.)
//    }

    public static ElementRule randomizeRules(boolean optional) {
        return new ElementRule(TreeTraitParserUtilities.RANDOMIZE, new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(TreeTraitParserUtilities.RANDOMIZE_LOWER, true),
                AttributeRule.newDoubleRule(TreeTraitParserUtilities.RANDOMIZE_UPPER, true),
                new ElementRule(Parameter.class)
        }, optional);
    }

    public static ElementRule jitterRules(boolean optional) {
        return new ElementRule(JITTER, new XMLSyntaxRule[]{
                AttributeRule.newDoubleArrayRule(WINDOW),
                AttributeRule.newBooleanRule(DUPLICATES, true),
                new ElementRule(Parameter.class),
                new ElementRule(TaxonList.class, true),
        }, optional);
    }

    public void jitter(XMLObject xo, int length, List<Integer> missingIndices, int paramDim) throws XMLParseException {
        jitter(xo, length, ContinuousTraitPartialsProvider.indicesToIndicator(missingIndices, paramDim));
    }

    public void jitter(XMLObject xo, int length, boolean[] missingIndicators) throws XMLParseException {
        XMLObject cxo = xo.getChild(TreeTraitParserUtilities.JITTER);
        Parameter traits = (Parameter) cxo.getChild(Parameter.class);
        double[] window = cxo.getDoubleArrayAttribute(TreeTraitParserUtilities.WINDOW); // Must be included, no default value
        boolean duplicates = cxo.getAttribute(TreeTraitParserUtilities.DUPLICATES, true); // default = true

        TaxonList taxonList = (TaxonList) cxo.getChild(TaxonList.class);
        jitter(traits, length, missingIndicators, window, taxonList, duplicates, true);
    }

    public void randomize(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(RANDOMIZE);
        Parameter traits = (Parameter) cxo.getChild(Parameter.class);
        double[] randomizeLower;
        double[] randomizeUpper;
        if (cxo.hasAttribute(RANDOMIZE_LOWER)) {
            randomizeLower = cxo.getDoubleArrayAttribute(RANDOMIZE_LOWER);
        } else {
            randomizeLower = new double[]{-90.0};
        }
        if (cxo.hasAttribute(RANDOMIZE_UPPER)) {
            randomizeUpper = cxo.getDoubleArrayAttribute(RANDOMIZE_UPPER);
        } else {
            randomizeUpper = new double[]{+90.0};
        }
        randomize(traits, randomizeLower, randomizeUpper);
    }

    private class DoubleArray implements Comparable {

        double[] value;
        int index;

        DoubleArray(double[] value, int index) {
            this.value = value;
            this.index = index;
        }

        public double[] getValues() {
            return value;
        }

        public int getIndex() {
            return index;
        }

        public int compareTo(Object o) {
            double[] x = ((DoubleArray) o).getValues();
            for (int i = 0; i < value.length; i++) {
                if (value[i] > x[i]) {
                    return 1;
                } else if (value[i] < x[i]) {
                    return -1;
                }
            }
            return 0;
        }
    }

    public boolean hasIdenticalTraits(Parameter trait, List<Integer> missingIndices, int dim) {
        int numTraits = trait.getDimension() / dim;

        List<DoubleArray> traitArray = new ArrayList<DoubleArray>();
        for (int i = 0; i < numTraits; i++) {

            if (!missingIndices.contains(i * dim)) { // TODO Assumes completely missing traits
                double[] x = new double[dim];
                for (int j = 0; j < dim; j++) {
                    x[j] = trait.getParameterValue(i * dim + j);
                }
                traitArray.add(new DoubleArray(x, i));
            }
        }

        DoubleArray[] sortedTraits = traitArray.toArray(new DoubleArray[0]);
        Arrays.sort(sortedTraits);
        // Mark duplicates
        for (int i = 1; i < sortedTraits.length; i++) {
            if (sortedTraits[i].compareTo(sortedTraits[i - 1]) == 0) {
                return true;
            }
        }
        return false;
    }


    public void jitter(Parameter trait, int dim, boolean[] missingIndicators, double[] window,
                       TaxonList taxonList,
                       boolean duplicates, boolean verbose) {
        int numTraits = trait.getDimension() / dim;
        boolean[] update = new boolean[numTraits];
        if (!duplicates) {
            Arrays.fill(update, true);
        } else {
            DoubleArray[] traitArray = new DoubleArray[numTraits];
            for (int i = 0; i < numTraits; i++) {
                double[] x = new double[dim];
                for (int j = 0; j < dim; j++) {
                    x[j] = trait.getParameterValue(i * dim + j);
                }
                traitArray[i] = new DoubleArray(x, i);
            }
            Arrays.sort(traitArray);
            // Mark duplicates
            for (int i = 1; i < numTraits; i++) {
                if (traitArray[i].compareTo(traitArray[i - 1]) == 0) {
                    update[traitArray[i - 1].getIndex()] = true;
                    update[traitArray[i].getIndex()] = true;
                }
            }
        }
        if (taxonList != null) {

            if (!(trait instanceof CompoundParameter)) {
                throw new IllegalArgumentException("Currently unable to match taxon names with a FastMatrixParameter");
            }

            Set<String> includedTaxonNames = new HashSet<>();
            for (int i = 0; i < taxonList.getTaxonCount(); ++i) {
                includedTaxonNames.add(taxonList.getTaxonId(i));
            }

            CompoundParameter cParameter = (CompoundParameter) trait;
            int offset = 0;
            for (int i = 0; i < cParameter.getParameterCount(); ++i) {
                String name = cParameter.getParameter(i).getParameterName();
                if (!includedTaxonNames.contains(name)) {
                    update[offset] = false;
                    if (verbose) {
                        Logger.getLogger("dr.evomodel.continuous").info(
                                "  Excluding taxon '" + name + "' from jitter.");
                    }
                }
                ++offset;
            }
        }
        for (int i = 0; i < numTraits; i++) {
            if (update[i]) {
                StringBuffer sb1 = null;
                StringBuffer sb2 = null;
                if (verbose) {
                    sb1 = new StringBuffer();
                    sb2 = new StringBuffer();
                }
                boolean hitAtLeastOneComponent = false;
                for (int j = 0; j < dim; j++) {
                    final double oldValue = trait.getParameterValue(i * dim + j);
                    final double newValue;
                    if (!missingIndicators[i * dim + j]) {
                        newValue = window[j % window.length] * (MathUtils.nextDouble() - 0.5) +
                                oldValue;
                        trait.setParameterValue(i * dim + j, newValue);
                        hitAtLeastOneComponent = true;
                    } else {
                        newValue = oldValue;
                    }
                    if (verbose) {
                        sb1.append(" ").append(oldValue);
                        sb2.append(" ").append(newValue);
                    }
                }
                if (verbose && hitAtLeastOneComponent) {
                    Logger.getLogger("dr.evomodel.continuous").info(
                            "  Replacing trait #" + (i + 1) + "  Old:" + sb1.toString() + " New: " + sb2.toString()
                    );
                }
            }
        }
    }

    public class TraitsAndMissingIndices {
        public CompoundParameter traitParameter;
        public boolean[] missingIndicators;
        public String traitName;
        public Parameter sampleMissingParameter;
        public boolean useMissingIndices;

        TraitsAndMissingIndices(CompoundParameter traitParameter, boolean[] missingIndicators, String traitName,
                                Parameter sampleMissingParameter, boolean useMissingIndices) {
            this.traitParameter = traitParameter;
            this.missingIndicators = missingIndicators;
            this.traitName = traitName;
            this.sampleMissingParameter = sampleMissingParameter;
            this.useMissingIndices = useMissingIndices;
        }

        public boolean[] getMissingIndicators() {
            return missingIndicators;
        }

        public List<Integer> getMissingIndices() { //TODO: deprecate
            return ContinuousTraitPartialsProvider.indicatorToIndices(missingIndicators);
        }
    }

    public static boolean isMissing(String oneValue) {
        return oneValue.equals("NA") || oneValue.equals("?");
    }

    public TraitsAndMissingIndices parseTraitsFromTaxonAttributes(
            XMLObject xo,
            Tree treeModel,
            boolean integrateOutInternalStates) throws XMLParseException {

        XMLObject xoc = xo.getChild(TRAIT_PARAMETER);
        Parameter parameter = (Parameter) xoc.getChild(Parameter.class);
        boolean existingTraitParameter = false;
        int randomSampleSizeFlag = xo.getAttribute(RANDOM_SAMPLE, -1);

        final String traitName;

        CompoundParameter traitParameter;
        boolean[] missingIndicators = null;
        int nMissing = 0;
        Parameter sampleMissingParameter = null;

        boolean isMatrixParameter = false;
        if (parameter instanceof MatrixParameter || parameter instanceof FastMatrixParameter) {
            traitParameter = (CompoundParameter) parameter;
            isMatrixParameter = true;
        } else if (parameter instanceof CompoundParameter) {
            // if we have been passed a CompoundParameter, this will be a leaf trait
            // parameter from a tree model so use this to allow for individual sampling
            // of leaf parameters.
            traitParameter = (CompoundParameter) parameter;
            existingTraitParameter = true;
        } else {
            // create a compound parameter of appropriate dimensions
            traitParameter = new CompoundParameter(parameter.getId());
            ParameterParser.replaceParameter(xoc, traitParameter);
        }

        if (xo.hasAttribute(TRAIT_NAME)) {

            Map<Integer, Integer> randomSample = null;
            traitName = xo.getStringAttribute(TRAIT_NAME);

            StringBuilder warnings = new StringBuilder();
            int warningLength = 0;
            final int maxWarnings = 10;

            // Fill in attributeValues
            int taxonCount = treeModel.getTaxonCount();
            for (int i = 0; i < taxonCount; i++) {
                String taxonName = treeModel.getTaxonId(i);

                // changed to just label the rows by the taxonName so it can be picked up elsewhere
                String paramName = taxonName;
                String altParamName = taxonName + "." + traitName;

                String object = (String) treeModel.getTaxonAttribute(i, traitName);
                if (object == null) {
                    throw new RuntimeException("Trait \"" + traitName + "\" not found for taxa \"" + taxonName + "\"");
                } else {
                    StringTokenizer st = new StringTokenizer(object);
                    int count = st.countTokens();

                    Parameter traitParam;
                    if (existingTraitParameter) {
                        traitParam = getTraitParameterByName(traitParameter, paramName);
                        if (traitParam == null) {
                            // try the alternative param name
                            traitParam = getTraitParameterByName(traitParameter, altParamName);
                            if (traitParam == null) {
                                throw new RuntimeException("Missing trait parameters for taxon, " + paramName);
                            }
                        }
                    } else {
                        if (isMatrixParameter) {
                            traitParam = traitParameter.getParameter(i);
                            traitParam.setId(paramName);
                        } else {
                            // Make multidimensional, in earlier revisions only first dimension was stored
                            traitParam = new Parameter.Default(paramName, count);
                            traitParameter.addParameter(traitParam);
                        }
                    }


                    int sampleSize = count;
                    if (randomSampleSizeFlag > 0) {
                        if (randomSample == null) {
                            randomSample = drawRandomSample(randomSampleSizeFlag, count);
                        }
                        sampleSize = randomSampleSizeFlag;
                    }
                    if (sampleSize != traitParam.getDimension()) {
                        if (existingTraitParameter) {
                            throw new RuntimeException("Trait length must match trait parameter dimension for taxon, " +
                                    taxonName + ": " +
                                    sampleSize + " != " + traitParam.getDimension());
                        } else {
                            traitParam.setDimension(sampleSize);
                        }
                    }
                    int index = 0;
                    for (int j = 0; j < count; j++) {
                        String oneValue = st.nextToken();
                        if (randomSampleSizeFlag == -1 || randomSample.containsKey(j)) {
                            double value = Double.NaN;
                            if (isMissing(oneValue)) {
//                            if (oneValue.equals("NA") || oneValue.equals("?") ) {
                                if (warningLength < maxWarnings) {
                                    warnings.append(
                                            "Warning: Missing value in tip for taxon " + taxonName +
                                                    " (filling with 0 as starting value when sampling only)\n"   // See comment below
                                    );
                                    ++warningLength;
                                }
                            } else {
                                try {
                                    value = new Double(oneValue);
                                    if (Double.isNaN(value)) {
                                        if (warningLength < maxWarnings) {
                                            warnings.append(
                                                    "Warning: Unrecognizable number " + oneValue + " for taxon " + taxonName + "\n"
                                            );
                                            ++warningLength;
                                        }
                                    }
                                } catch (NumberFormatException e) {
                                    throw new RuntimeException(e.getMessage());
                                }
                            }

                            int replicates = 1;
                            if (randomSampleSizeFlag != -1) {
                                // Count how many times to add this datum
                                replicates = randomSample.get(j);
                            }
                            for (int k = 0; k < replicates; k++) {
                                traitParam.setParameterValue(index, value);
                                index++;
                            }
                        }
                    }
                }
            }

            if (warningLength > 0) {
                Logger.getLogger("dr.evomodel.continuous").info(warnings.toString());
                if (warningLength == maxWarnings) {
                    Logger.getLogger("dr.evomodel.continuous").info("Warning: only first " + maxWarnings + " trait warnings were displayed\n");
                }
            }

            // Find missing values
            double[] allValues = traitParameter.getParameterValues();

            missingIndicators = new boolean[allValues.length];
            for (int i = 0; i < allValues.length; i++) {
                if ((new Double(allValues[i])).isNaN()) {
                    traitParameter.setParameterValue(i, 0); // Here, missings are set to zero
                    missingIndicators[i] = true;
                    nMissing++;
                }
            }

            // Standardize
            if (xo.getAttribute(STANDARDIZE, false) && traitParameter instanceof MatrixParameterInterface) {

                double targetSd = xo.getAttribute(TARGET_SD, 1.0);

                StandardizeTraits st = new StandardizeTraits((MatrixParameterInterface) traitParameter,
                        missingIndicators, targetSd);
                String message = st.doStandardization(false);

                Logger.getLogger("dr.evomodel.continous").info(message);
            }

            if (xo.hasChildNamed(MISSING)) {
                XMLObject cxo = xo.getChild(MISSING);

                Parameter missingParameter = new Parameter.Default(allValues.length, 0.0);
                for (int i = 0; i < missingIndicators.length; i++) {
                    if (missingIndicators[i]) {
                        missingParameter.setParameterValue(i, 1.0);
                    }
                }

                if (cxo.hasAttribute(LATENT_FROM) && cxo.hasAttribute(LATENT_TO)) {
                    int from = cxo.getIntegerAttribute(LATENT_FROM);
                    int to = cxo.getIntegerAttribute(LATENT_TO);

                    final int dimTrait = allValues.length / taxonCount;

                    if (from < 1 || to < 1 || from > dimTrait || to > dimTrait) {
                        throw new XMLParseException("Invalid latent dimension specification");
                    }

                    int index = 0;
                    for (int taxon = 0; taxon < taxonCount; ++taxon) {
                        for (int trait = from - 1; trait < to; ++trait) {
                            missingParameter.setParameterValue(index + trait, 1.0);
                        }
                        index += dimTrait;
                    }
                }

                missingParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, allValues.length));
                ParameterParser.replaceParameter(cxo, missingParameter);
                sampleMissingParameter = missingParameter;
            }

            // Give warnings if trait exist for internal and root nodes when integrating them out
            if (false) {
                int numTraits = traitParameter.getParameterCount();
                if (numTraits != treeModel.getExternalNodeCount()) {
                    throw new XMLParseException(
                            "Dimensionality of '" + traitParameter.getId() + "' (" + numTraits + ") is not equal to the number" +
                                    " of tree tips (" + treeModel.getExternalNodeCount() + ")");
                }

                for (int j = 0; j < numTraits; j++) {
                    String parameterName = traitParameter.getParameter(j).getId();
                    if (parameterName.startsWith("node") || parameterName.startsWith("root")) {
                        throw new XMLParseException(
                                "Internal/root node trait parameters are not allowed when " +
                                        "using the integrated observed data likelihood");
                    }
                }
            }
        } else {
            traitName = DEFAULT_TRAIT_NAME;
        }

        boolean useMissingIndices = true;
        if (xo.getAttribute(SAMPLE_MISSING_TRAITS, false) || xo.hasChildNamed(MISSING)) {
//            missingIndices = new ArrayList<Integer>(); // return empty
            useMissingIndices = false;

        }

        if (missingIndicators == null || nMissing == 0) {
            useMissingIndices = false;
        }

        return new TraitsAndMissingIndices(traitParameter, missingIndicators, traitName,
                sampleMissingParameter, useMissingIndices);
    }

    public static List<Integer> parseMissingIndices(Parameter traitParameter, double[] allValues) {
        List<Integer> missingIndices = new ArrayList<Integer>();
        for (int i = 0; i < allValues.length; i++) {
            if ((new Double(allValues[i])).isNaN()) {
                traitParameter.setParameterValue(i, 0); // Here, missings are set to zero
                missingIndices.add(i);
            }
        }
        return missingIndices;
    }


    private Parameter getTraitParameterByName(CompoundParameter traits, String name) {

        for (int i = 0; i < traits.getParameterCount(); i++) {
            Parameter found = traits.getParameter(i);
            if (found.getStatisticName().compareTo(name) == 0)
                return found;
        }
        return null;
    }

    private Map<Integer, Integer> drawRandomSample(int total, int length) {
        Map<Integer, Integer> thisMap = new HashMap<Integer, Integer>(total);
        for (int i = 0; i < total; i++) {
            int item = MathUtils.nextInt(length);
            if (thisMap.containsKey(item)) {
                thisMap.put(item, thisMap.get(item) + 1);
            } else {
                thisMap.put(item, 1);
            }
        }
        return thisMap;
    }

    public static TreeTrait getTreeTraitFromDataLikelihood(TreeDataLikelihood dataLikelihood) {
        return dataLikelihood.getTreeTrait(getTipTraitNameFromDataLikelihood(dataLikelihood));
    }

    public static String getTipTraitNameFromDataLikelihood(TreeDataLikelihood dataLikelihood) {
        ContinuousDataLikelihoodDelegate delegate =
                (ContinuousDataLikelihoodDelegate) dataLikelihood.getDataLikelihoodDelegate();

        ContinuousTraitPartialsProvider dataModel = delegate.getDataModel();
        String traitName = dataModel.getTipTraitName();
        return traitName;
    }
}
