/*
 * TreeTraitParserUtilities.java
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

package dr.evomodelxml.treelikelihood;

import dr.evolution.tree.MultivariateTraitTree;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;

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

    public void randomize(Parameter trait, double[] lower, double[] upper) {
        // Draws each dimension in each trait from U[lower, upper)
        for (int i = 0; i < trait.getDimension(); i++) {
            final int whichLower = i % lower.length;
            final int whichUpper = i % upper.length;
            final double newValue = MathUtils.uniform(lower[whichLower], upper[whichUpper]);
            trait.setParameterValue(i, newValue);
        }
    }

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

        }, optional);
    }

    public void jitter(XMLObject xo, int length, List<Integer> missingIndices) throws XMLParseException {
        XMLObject cxo = xo.getChild(TreeTraitParserUtilities.JITTER);
        Parameter traits = (Parameter) cxo.getChild(Parameter.class);
        double[] window = cxo.getDoubleArrayAttribute(TreeTraitParserUtilities.WINDOW); // Must be included, no default value
        boolean duplicates = cxo.getAttribute(TreeTraitParserUtilities.DUPLICATES, true); // default = true
        jitter(traits, length, missingIndices, window, duplicates, true);
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

    public boolean hasIdenticalTraits(Parameter trait, int dim) {
        int numTraits = trait.getDimension() / dim;

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
                return true;
            }
        }
        return false;
    }

    public void jitter(Parameter trait, int dim, List<Integer> missingIndices, double[] window, boolean duplicates, boolean verbose) {
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
                    if (!missingIndices.contains(i * dim + j)) {
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
        public List<Integer> missingIndices;
        public String traitName;

        TraitsAndMissingIndices(CompoundParameter traitParameter, List<Integer> missingIndices, String traitName) {
            this.traitParameter = traitParameter;
            this.missingIndices = missingIndices;
            this.traitName = traitName;
        }
    }

    public TraitsAndMissingIndices parseTraitsFromTaxonAttributes(
            XMLObject xo,
            String inTraitName,
            MultivariateTraitTree treeModel,
            boolean integrateOutInternalStates) throws XMLParseException {

        XMLObject xoc = xo.getChild(TRAIT_PARAMETER);
        Parameter parameter = (Parameter) xoc.getChild(Parameter.class);
        boolean existingTraitParameter = false;
        int randomSampleSizeFlag = xo.getAttribute(RANDOM_SAMPLE, -1);

        String traitName = inTraitName;

        CompoundParameter traitParameter;
        List<Integer> missingIndices = null;

        boolean isMatrixParameter = false;
        if (parameter instanceof MatrixParameter) {
            traitParameter = (CompoundParameter) parameter;
            isMatrixParameter = true;
        } else


        if (parameter instanceof CompoundParameter) {
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
                            if (oneValue.equals("NA") || oneValue.equals("?") ) {
                                Logger.getLogger("dr.evomodel.continuous").info(
                                        "Warning: Missing value in tip for taxon " + taxonName +
                                                " (filling with 0 as starting value when sampling only)"   // See comment below
                                );
                            } else {
                                try {
                                    value = new Double(oneValue);
                                    if (Double.isNaN(value)) {
                                        Logger.getLogger("dr.evomodel.continuous").info(
                                                "Warning: Unrecognizable number " + oneValue + " for taxon " + taxonName
                                        );
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

            // Find missing values
            double[] allValues = traitParameter.getParameterValues();
            missingIndices = new ArrayList<Integer>();
            for (int i = 0; i < allValues.length; i++) {
                if ((new Double(allValues[i])).isNaN()) {
                    traitParameter.setParameterValue(i, 0); // Here, missings are set to zero
                    missingIndices.add(i);
                }
            }

            if (xo.hasChildNamed(MISSING)) {
                XMLObject cxo = xo.getChild(MISSING);
                Parameter missingParameter = new Parameter.Default(allValues.length, 0.0);
                for (int i : missingIndices) {
                    missingParameter.setParameterValue(i, 1.0);
                }
                missingParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, allValues.length));
                ParameterParser.replaceParameter(cxo, missingParameter);
            }

            // Give warnings if trait exist for internal and root nodes when integrating them out
            if (integrateOutInternalStates) {
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
                                        "using the integrated observed data multivariateTraitLikelihoood");
                    }
                }
            }
        }
        return new TraitsAndMissingIndices(traitParameter, missingIndices, traitName);
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
}
