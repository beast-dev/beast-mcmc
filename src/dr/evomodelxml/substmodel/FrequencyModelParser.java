/*
 * FrequencyModelParser.java
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

package dr.evomodelxml.substmodel;

//import dr.app.bss.Utils;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.HiddenDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.CodonFromNucleotideFrequencyModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.*;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class FrequencyModelParser extends AbstractXMLObjectParser {

    public static final String FREQUENCIES = "frequencies";
    public static final String FREQUENCY_MODEL = "frequencyModel";
    public static final String NORMALIZE = "normalize";
    private static final String COMPRESS = "compress";

    private static final String COMPOSITION = "composition";
    private static final String FREQ_3x4 = "3x4";
    private static final String CODON_FROM_NUCLEOTIDE = "codonFromNucleotide";
    private static final String[] COMPOSITION_TYPES = new String[]{FREQ_3x4};

    public String getParserName() {
        return FREQUENCY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);

        Parameter freqsParam = (Parameter) xo.getElementFirstChild(FREQUENCIES);
        double[] frequencies = null;

        List<PatternList> multiplePatterns = new ArrayList<>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object obj = xo.getChild(i);
            if (obj instanceof PatternList) {

                PatternList patternList = (PatternList) obj;
                if (xo.getAttribute(COMPRESS, false)
                        && (patternList.getDataType() instanceof HiddenDataType)) {

                    double[] hiddenFrequencies = patternList
                            .getStateFrequencies();
                    int hiddenCount = ((HiddenDataType) patternList
                            .getDataType()).getHiddenClassCount();
                    int baseStateCount = hiddenFrequencies.length / hiddenCount;
                    frequencies = new double[baseStateCount];
                    for (int j = 0; j < baseStateCount; ++j) {
                        for (int k = 0; k < hiddenCount; ++k) {
                            frequencies[j] += hiddenFrequencies[j + k
                                    * baseStateCount];
                        }
                    }
                    break;
                } else {

                    // TODO
                    if (xo.hasAttribute(COMPOSITION)) {

                        String type = xo.getStringAttribute(COMPOSITION);
                        if (type.equalsIgnoreCase(FREQ_3x4)) {

                            frequencies = getEmpirical3x4Freqs(patternList);

                        }
                        break;
                    } else {

                        multiplePatterns.add(patternList);
                        frequencies = patternList.getStateFrequencies();

                    } // END: composition check

                }
            }// END: patternList check
        }

        if (multiplePatterns.size() > 1) {
            StringBuilder sb = new StringBuilder("\n Concatenating multiple pattern lists to calculate initial frequencies from: \n\t");
            for (PatternList patterns : multiplePatterns) {
                sb.append(" " + patterns.getId());
            }

            int totalWeight = 0;
            for (PatternList p : multiplePatterns) {
                totalWeight += p.getPatternCount();
            }

            Arrays.fill(frequencies, 0.0);

            for (PatternList p : multiplePatterns) {
                double[] increment = p.getStateFrequencies();
                double weight = (double) p.getPatternCount() / (double) totalWeight;

                for (int i = 0; i < frequencies.length; ++i) {
                    frequencies[i] += increment[i] * weight;
                }
            }

            Logger.getLogger("dr.evomodel").info(sb.toString());
        }

        StringBuilder sb = new StringBuilder(
                "\nCreating state frequencies model '"
                        + freqsParam.getParameterName() + "': ");
        if (frequencies != null) {

            if (freqsParam.getDimension() != frequencies.length) {
                throw new XMLParseException(
                        "dimension of frequency parameter and number of sequence states don't match.");
            }

            for (int j = 0; j < frequencies.length; j++) {
                freqsParam.setParameterValue(j, frequencies[j]);
            }

            sb.append("Using empirical frequencies from data ");

        } else {
            sb.append("Initial frequencies ");
        }
        sb.append("= {");

        if (xo.getAttribute(NORMALIZE, false)) {
            double sum = 0;
            for (int j = 0; j < freqsParam.getDimension(); j++)
                sum += freqsParam.getParameterValue(j);
            for (int j = 0; j < freqsParam.getDimension(); j++) {
                if (sum != 0)
                    freqsParam.setParameterValue(j,
                            freqsParam.getParameterValue(j) / sum);
                else
                    freqsParam.setParameterValue(j,
                            1.0 / freqsParam.getDimension());
            }
        }

        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(5);

        sb.append(format.format(freqsParam.getParameterValue(0)));
        for (int j = 1; j < freqsParam.getDimension(); j++) {
            sb.append(", ");
            sb.append(format.format(freqsParam.getParameterValue(j)));
        }
        sb.append("}");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        if (dataType instanceof Codons && xo.getAttribute(CODON_FROM_NUCLEOTIDE, false)) {

            FrequencyModel nucleotideFrequencyModel = new FrequencyModel(Nucleotides.INSTANCE, freqsParam);
            Parameter codonFrequencies = new Parameter.Default(dataType.getStateCount(), 1.0 / dataType.getStateCount());
            codonFrequencies.setId(xo.getId());

            return new CodonFromNucleotideFrequencyModel((Codons) dataType, nucleotideFrequencyModel, codonFrequencies);

        } else {
            return new FrequencyModel(dataType, freqsParam);
        }
    }// END: parseXMLObject

    private double[] getEmpirical3x4Freqs(PatternList patternList) {

        DataType nucleotideDataType = Nucleotides.INSTANCE;
        Codons codonDataType = Codons.UNIVERSAL;
        List<String> stopCodonsList = Arrays.asList(STOP_CODONS);

        int cStateCount = codonDataType.getStateCount();
        int nStateCount = nucleotideDataType.getStateCount();
        int nPosition = 3;

        double[] stopCodonFreqs = new double[STOP_CODONS.length];
        double counts[][] = new double[nStateCount][nPosition];
        int countsPos[] = new int[nPosition];

        int patternCount = patternList.getPatternCount();
        for (int i = 0; i < patternCount; i++) {

            int[] sitePatterns = patternList.getPattern(i);
            for (int codonState : sitePatterns) {

                int[] nucleotideStates = codonDataType.getTripletStates(codonState);

                String triplet = codonDataType.getTriplet(codonState);

//				if(triplet.equals("TAA")) {
//					System.err.println(triplet);
//				}

                if (stopCodonsList.contains(triplet)) {

                    int stopCodonIndex = stopCodonsList.indexOf(triplet);
                    stopCodonFreqs[stopCodonIndex]++;

                } //END: stopCodon check

                for (int pos = 0; pos < nPosition; pos++) {

                    int nucleotideState = nucleotideStates[pos];
                    counts[nucleotideState][pos]++;
                    countsPos[pos]++;

                } //END: nucleotide positions loop

            } //END: sitePatterns loop
        } // sites loop

        int total = 0;
        for (int pos = 0; pos < nPosition; pos++) {

            int totalPos = countsPos[pos];
            for (int s = 0; s < nStateCount; s++) {
                counts[s][pos] = counts[s][pos] / totalPos;
            }//END: nucleotide states loop

            total += totalPos;
        }//END: nucleotide positions loop

//		Utils.printArray(stopCodonFreqs);
//		System.out.println(stopCodonFreqs.length);

        // add stop codon frequencies
        double pi = 0.0;
        for (double stopCodonFreq : stopCodonFreqs) {
            double freq = stopCodonFreq / total;
            pi += freq;
        }

//		System.out.println(pi);

        double[] freqs = new double[cStateCount];
        Arrays.fill(freqs, 1.0);
        for (int codonState = 0; codonState < cStateCount; codonState++) {

            int[] nucleotideStates = codonDataType.getTripletStates(codonState);
            for (int pos = 0; pos < nPosition; pos++) {

                int nucleotide = nucleotideStates[pos];
                freqs[codonState] *= counts[nucleotide][pos];

            }// END: nucleotide positions loop

            // TODO: stop codons freqs
            freqs[codonState] = freqs[codonState] / (1 - pi);

        }//END: codon states loop


        return freqs;
    }// END: getEmpirical3x4Freqs


    public String getParserDescription() {
        return "A model of equilibrium base frequencies.";
    }

    public Class getReturnType() {
        return FrequencyModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

            new StringAttributeRule(COMPOSITION, "Composition type",
                    COMPOSITION_TYPES, true),

            AttributeRule.newBooleanRule(NORMALIZE, true),
            AttributeRule.newBooleanRule(COMPRESS, true),

            new ElementRule(PatternList.class, "Initial value", 0, 1),

            new XORRule(new StringAttributeRule(DataType.DATA_TYPE,
                    "The type of sequence data",
                    DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)),

            new ElementRule(FREQUENCIES, new XMLSyntaxRule[]{new ElementRule(
                    Parameter.class)}),
    };

    public static final String STOP_CODONS[] = new String[] { "TAA", "TAG", "TGA" };
}