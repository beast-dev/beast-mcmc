/*
 * SimpleAlignment.java
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

package dr.evolution.alignment;

import dr.app.bss.XMLExporter;
import dr.app.tools.NexusExporter;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.sequence.UncertainSequence;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.NumberFormatter;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * A simple alignment class that implements gaps by characters in the sequences.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SimpleAlignment.java,v 1.46 2005/06/21 16:25:15 beth Exp $
 */
@SuppressWarnings("serial")
public class SimpleAlignment extends Sequences implements Alignment, dr.util.XHTMLable {

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private OutputType outputType = OutputType.FASTA;
    private DataType dataType = null;
    private int siteCount = 0;
    private boolean siteCountKnown = false;
    private boolean countStatistics = !(dataType instanceof Codons) && !(dataType instanceof GeneralDataType);

    // **************************************************************
    // SimpleAlignment METHODS
    // **************************************************************

    /**
     * parameterless constructor.
     */
    public SimpleAlignment() {
    }

    /**
     * Constructs a sub alignment based on the provided taxa.
     *
     * @param a
     * @param taxa
     */
    public SimpleAlignment(Alignment a, TaxonList taxa) {

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);

            Sequence sequence = a.getSequence(a.getTaxonIndex(taxon));
            addSequence(sequence);
        }
    }

    public void setOutputType(OutputType out) {
        outputType = out;
    }

    public List<Sequence> getSequences() {
        return Collections.unmodifiableList(sequences);
    }

    /**
     * Calculates the siteCount by finding the longest sequence.
     */
    public void updateSiteCount() {
        siteCount = 0;
        int i, len, n = getSequenceCount();

        for (i = 0; i < n; i++) {
            len = getSequence(i).getLength();
            if (len > siteCount)
                siteCount = len;
        }

        siteCountKnown = true;
    }

    // **************************************************************
    // Alignment IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the dataType of this alignment. This should be the same as
     * the sequences.
     */
    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * @return number of sites
     */
    public int getSiteCount(DataType dataType) {
        return getSiteCount();
    }

    /**
     * sequence character at (sequence, site)
     */
    public char getChar(int sequenceIndex, int siteIndex) {
        return getSequence(sequenceIndex).getChar(siteIndex);
    }

    /**
     * Returns string representation of single sequence in
     * alignment with gap characters included.
     */
    public String getAlignedSequenceString(int sequenceIndex) {
        return getSequence(sequenceIndex).getSequenceString();
    }


    /**
     * Returns string representation of single sequence in
     * alignment with gap characters excluded.
     */
    public String getUnalignedSequenceString(int sequenceIndex) {

        StringBuffer unaligned = new StringBuffer();
        for (int i = 0, n = getSiteCount(); i < n; i++) {

            int state = getState(sequenceIndex, i);
            if (!dataType.isGapState(state)) {
                unaligned.append(dataType.getChar(state));
            }
        }

        return unaligned.toString();
    }

    // **************************************************************
    // Sequences METHODS
    // **************************************************************

    /**
     * Add a sequence to the sequence list
     */
    public void addSequence(Sequence sequence) {
        if (dataType == null) {
            if (sequence.getDataType() == null) {
                dataType = sequence.guessDataType();
                sequence.setDataType(dataType);
            } else {
                setDataType(sequence.getDataType());
            }
        } else if (sequence.getDataType() == null) {
            sequence.setDataType(dataType);
        } else if (dataType != sequence.getDataType()) {
            throw new IllegalArgumentException("Sequence's dataType does not match the alignment's");
        }
        
        int invalidCharAt = sequence.getInvalidChar();
        if (invalidCharAt >= 0)
            throw new IllegalArgumentException("Sequence of " + sequence.getTaxon().getId()
                    + " contains invalid char \'" + sequence.getChar(invalidCharAt) + "\' at index " + invalidCharAt);

        super.addSequence(sequence);
        updateSiteCount();
    }

    /**
     * Insert a sequence to the sequence list at position
     */
    public void insertSequence(int position, Sequence sequence) {
        if (dataType == null) {
            if (sequence.getDataType() == null) {
                dataType = sequence.guessDataType();
                sequence.setDataType(dataType);
            } else {
                setDataType(sequence.getDataType());
            }
        } else if (sequence.getDataType() == null) {
            sequence.setDataType(dataType);
        } else if (dataType != sequence.getDataType()) {
            throw new IllegalArgumentException("Sequence's dataType does not match the alignment's");
        }

        int invalidCharAt = sequence.getInvalidChar();
        if (invalidCharAt >= 0)
            throw new IllegalArgumentException("Sequence of " + sequence.getTaxon().getId()
                    + " contains invalid char \'" + sequence.getChar(invalidCharAt) + "\' at index " + invalidCharAt);

        super.insertSequence(position, sequence);
    }
    
    // **************************************************************
    // SiteList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of sites
     */
    public int getSiteCount() {
        if (!siteCountKnown)
            updateSiteCount();
        return siteCount;
    }

    /**
     * Gets the pattern of site as an array of state numbers (one per sequence)
     *
     * @return the site pattern at siteIndex
     */
    public int[] getSitePattern(int siteIndex) {
        Sequence seq;
        int i, n = getSequenceCount();

        int[] pattern = new int[n];

        for (i = 0; i < n; i++) {
            seq = getSequence(i);

            if (siteIndex >= seq.getLength())
                pattern[i] = dataType.getGapState();
            else
                pattern[i] = seq.getState(siteIndex);
        }

        return pattern;
    }

    @Override
    public double[][] getUncertainSitePattern(int siteIndex) {
        if (areUncertain())   {

            double[][] pattern = new double[getSequenceCount()][];
            for (int i = 0; i < getSequenceCount(); ++i) {

                Sequence seq = getSequence(i);
                if (siteIndex > seq.getLength()) {
                    pattern[i] = new double[dataType.getStateCount()];
                    Arrays.fill(pattern[i], 1.0);
                } else {
                    if (seq instanceof UncertainSequence) {
                        pattern[i] = ((UncertainSequence) seq).getUncertainPattern(siteIndex);
                    } else {
                        pattern[i] = new double[dataType.getStateCount()];
                        int[] states = dataType.getStates(seq.getState(siteIndex));
                        for (int state : states) {
                            pattern[i][state] = 1.0;
                        }
                    }
                }
            }

            return pattern;
        } else {
            throw new UnsupportedOperationException("getUncertainSitePattern not implemented yet");
        }
    }

    /**
     * Gets the pattern index at a particular site
     *
     * @return the patternIndex
     */
    public int getPatternIndex(int siteIndex) {
        return siteIndex;
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public int getState(int taxonIndex, int siteIndex) {

        Sequence seq = getSequence(taxonIndex);

        if (siteIndex >= seq.getLength()) {
            return dataType.getGapState();
        }

        return seq.getState(siteIndex);
    }

    @Override
    public double[] getUncertainState(int taxonIndex, int siteIndex) {
        throw new UnsupportedOperationException("getUncertainState not implemented yet");
    }

    /**
     */
    public void setState(int taxonIndex, int siteIndex, int state) {

        Sequence seq = getSequence(taxonIndex);

        if (siteIndex >= seq.getLength()) {
            throw new IllegalArgumentException();
        }

        seq.setState(siteIndex, state);
    }

    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of patterns
     */
    public int getPatternCount() {
        return getSiteCount();
    }

    /**
     * @return number of invariant sites
     */
    public int getInvariantCount() {
        int invariantSites = 0;
        for (int i = 0; i < getSiteCount(); i++) {
            int[] pattern = getSitePattern(i);
            if (Patterns.isInvariant(pattern)) {
                invariantSites++;
            }
        }
        return invariantSites;
    }

    public int getUniquePatternCount() {
        Patterns patterns = new Patterns(this);
        return patterns.getPatternCount();
    }

    public int getInformativeCount() {
        Patterns patterns = new Patterns(this);
        int informativeCount = 0;

        for (int i = 0; i < patterns.getPatternCount(); i++) {
            int[] pattern = patterns.getPattern(i);

            if (isInformative(pattern)) {
                informativeCount += patterns.getPatternWeight(i);
            }
        }

        return informativeCount;
    }

    public int getSingletonCount() {
        Patterns patterns = new Patterns(this);
        int singletonCount = 0;

        for (int i = 0; i < patterns.getPatternCount(); i++) {
            int[] pattern = patterns.getPattern(i);

            if (!Patterns.isInvariant(pattern) && !isInformative(pattern)) {
                singletonCount += patterns.getPatternWeight(i);
            }
        }

        return singletonCount;
    }

    private boolean isInformative(int[] pattern) {
        int[] stateCounts = new int[getStateCount()];
        for (int j = 0; j < pattern.length; j++) {
            stateCounts[pattern[j]]++;
        }

        boolean oneStateGreaterThanOne = false;
        boolean secondStateGreaterThanOne = false;
        for (int j = 0; j < stateCounts.length; j++) {
            if (stateCounts[j] > 1) {
                if (!oneStateGreaterThanOne) {
                    oneStateGreaterThanOne = true;
                } else {
                    secondStateGreaterThanOne = true;
                }

            }
        }

        return secondStateGreaterThanOne;
    }

    /**
     * @return number of states for this siteList
     */
    public int getStateCount() {
        return getDataType().getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    public int getPatternLength() {
        return getSequenceCount();
    }

    /**
     * Gets the pattern as an array of state numbers (one per sequence)
     *
     * @return the pattern at patternIndex
     */
    public int[] getPattern(int patternIndex) {
        return getSitePattern(patternIndex);
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        throw new UnsupportedOperationException("getUncertainPattern not implemented yet");
    }

    /**
     * @return state at (taxonIndex, patternIndex)
     */
    public int getPatternState(int taxonIndex, int patternIndex) {
        return getState(taxonIndex, patternIndex);
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        throw new UnsupportedOperationException("getUncertainPatternState not implemented yet");
    }

    /**
     * Gets the weight of a site pattern (always 1.0)
     */
    public double getPatternWeight(int patternIndex) {
        return 1.0;
    }

    /**
     * @return the array of pattern weights
     */
    public double[] getPatternWeights() {
        double[] weights = new double[siteCount];
        for (int i = 0; i < siteCount; i++)
            weights[i] = 1.0;
        return weights;
    }

    /**
     * @return the DataType of this siteList
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * @return the frequency of each state
     */
    public double[] getStateFrequencies() {
        return PatternList.Utils.empiricalStateFrequencies(this);
    }

    @Override
    public boolean areUnique() {
        return false;
    }

    @Override
    public boolean areUncertain() {

        for (Sequence seq : sequences) {
            if (seq instanceof UncertainSequence) {
                return true;
            }
        }
        return false;
    }

    public void setReportCountStatistics(boolean report) {
        countStatistics = report;
    }

    public String toString() {
        return outputType.makeOutputString(this);   // generic delegation to ease extensibility
    }// END: toString

    public String toXHTML() {
        String xhtml = "<p><em>Alignment</em> data type = ";
        xhtml += getDataType().getDescription();
        xhtml += ", no. taxa = ";
        xhtml += getTaxonCount();
        xhtml += ", no. sites = ";
        xhtml += getSiteCount();
        xhtml += "</p>";

        xhtml += "<pre>";

        int length, maxLength = 0;
        for (int i = 0; i < getTaxonCount(); i++) {
            length = getTaxonId(i).length();
            if (length > maxLength)
                maxLength = length;
        }

        for (int i = 0; i < getTaxonCount(); i++) {
            length = getTaxonId(i).length();
            xhtml += getTaxonId(i);
            for (int j = length; j <= maxLength; j++)
                xhtml += " ";
            xhtml += getAlignedSequenceString(i) + "\n";
        }
        xhtml += "</pre>";
        return xhtml;
    }

    public enum OutputType {
        FASTA("fasta", "fsa") {
            @Override
            public String makeOutputString(SimpleAlignment alignment) {

                NumberFormatter formatter = new NumberFormatter(6);
                StringBuffer buffer = new StringBuffer();

                if (alignment.countStatistics) {
                    buffer.append("Site count = ").append(alignment.getSiteCount()).append("\n");
                    buffer.append("Invariant sites = ").append(alignment.getInvariantCount()).append("\n");
                    buffer.append("Singleton sites = ").append(alignment.getSingletonCount()).append("\n");
                    buffer.append("Parsimony informative sites = ").append(alignment.getInformativeCount()).append("\n");
                    buffer.append("Unique site patterns = ").append(alignment.getUniquePatternCount()).append("\n\n");
                }

                for (int i = 0; i < alignment.getSequenceCount(); i++) {
                    String name = formatter.formatToFieldWidth(alignment.getTaxonId(i), 10);
                    buffer.append(">" + name + "\n");
                    buffer.append(alignment.getAlignedSequenceString(i) + "\n");
                }

                return buffer.toString();
            }
        },

        NEXUS("nexus", "nxs") {
            @Override
            public String makeOutputString(SimpleAlignment alignment) {

                StringBuffer buffer = new StringBuffer();

                try {

                    File tmp = File.createTempFile("tempfile", ".tmp");
                    PrintStream ps = new PrintStream(tmp);

                    NexusExporter nexusExporter = new NexusExporter(ps);
                    buffer.append(nexusExporter.exportAlignment(alignment));

                } catch (IllegalArgumentException e) {
                    e.printStackTrace();
                } catch (IOException e) {
                    e.printStackTrace();
                }

                return buffer.toString();
            }// END: makeOutputString
        },
        XML("xml", "xml") {
            @Override
            public String makeOutputString(SimpleAlignment alignment) {

                StringBuffer buffer = new StringBuffer();

                try {

                    XMLExporter xmlExporter = new XMLExporter();
                    buffer.append(xmlExporter.exportAlignment(alignment));

                } catch (IllegalArgumentException e) {
                    e.printStackTrace();
                } catch (IOException e) {
					e.printStackTrace();
				}

                return buffer.toString();
            }// END: makeOutputString
        };

//        public static OutputType getValue(String str) {
//			if (FASTA.name().equalsIgnoreCase(str)) {
//				return FASTA;
//			} else if (NEXUS.name().equalsIgnoreCase(str)) {
//				return NEXUS;
//			} else if (XML.name().equalsIgnoreCase(str)) {
//				return XML;
//			}
//			return null;
//		}// END: getValue
//        
//        public static Enum[] getValues() {
//        	
//        	Enum values[] = new Enum[values().length];
//        	
//        	int i = 0;
//        	for(Enum value : OutputType.values()) {
//        		
//        		values[i] = getValue(value.toString());// value;
//        		i++;
//        		
//        	}
//        	
//        	return values;
//        }
        
        private final String text;
        private final String extension;

        private OutputType(String text, String extension) {
            this.text = text;
            this.extension = extension;
        }

        public String getText() {
            return text;
        }

        public String getExtension() {
            return extension;
        }

        public abstract String makeOutputString(SimpleAlignment alignment);

        public static OutputType parseFromString(String text) {
            for (OutputType type : OutputType.values()) {
                if (type.getText().compareToIgnoreCase(text) == 0) {
                    return type;
                }
            }
            return null;
        }

        public static OutputType parseFromExtension(String extension) {
            for (OutputType type : OutputType.values()) {
                if (type.getExtension().compareToIgnoreCase(extension) == 0) {
                    return type;
                }
            }
            return null;
        }
    }
}// END: class
