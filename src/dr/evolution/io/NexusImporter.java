/*
 * NexusImporter.java
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

package dr.evolution.io;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.TwoStates;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.SequenceList;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.*;
import dr.util.Attributable;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class for importing NEXUS file format
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: NexusImporter.java,v 1.30 2006/04/25 14:39:37 rambaut Exp $
 */
public class NexusImporter extends Importer implements SequenceImporter, TreeImporter {

    public static final NexusBlock UNKNOWN_BLOCK = new NexusBlock("unknown");
    public static final NexusBlock TAXA_BLOCK = new NexusBlock("TAXA");
    public static final NexusBlock CHARACTERS_BLOCK = new NexusBlock("CHARACTERS");
    public static final NexusBlock DATA_BLOCK = new NexusBlock("DATA");
    public static final NexusBlock UNALIGNED_BLOCK = new NexusBlock("UNALIGNED");
    public static final NexusBlock DISTANCES_BLOCK = new NexusBlock("DISTANCES");
    public static final NexusBlock TREES_BLOCK = new NexusBlock("TREES");
    public static final NexusBlock CALIBRATION_BLOCK = new NexusBlock("CALIBRATION");

    public static boolean suppressWarnings = false;

    private final boolean ignoreMetaComments;

    public static void setSuppressWarnings(boolean sw) {
        suppressWarnings = sw;
    }

    // NEXUS specific ImportException classes
    public static class MissingBlockException extends ImportException {
        /**
         *
         */
        private static final long serialVersionUID = -6287423449717453999L;

        public MissingBlockException() {
            super();
        }

        public MissingBlockException(String message) {
            super(message);
        }
    }

    public static class NexusBlock {
        private final String name;

        public NexusBlock(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }
    }

    /**
     * Constructor
     */
    public NexusImporter(Reader reader) {
        super(reader);
        this.ignoreMetaComments = false;
        setCommentDelimiters('[', ']', '\0', '\0', '&');
    }

    public NexusImporter(Reader reader, boolean ignoreMetaComments) {
        super(reader);
        this.ignoreMetaComments = ignoreMetaComments;
        setCommentDelimiters('[', ']', '\0', '!', '&');
    }

    public NexusImporter(Reader reader, Writer commentWriter) {
        super(reader, commentWriter);
        this.ignoreMetaComments = false;
        setCommentDelimiters('[', ']', '\0', '!', '&');
    }

    /**
     * This function returns an integer to specify what the
     * next block in the file is. The internal variable nextBlock is also set to this
     * value. This should be overridden to provide support for other blocks. Once
     * the block is read in, nextBlock is automatically set to UNKNOWN_BLOCK by
     * findEndBlock.
     */
    public NexusBlock findNextBlock() throws IOException {
        findToken("BEGIN", true);
        String blockName = readToken(";");
        return findBlockName(blockName);
    }

    /**
     * This function returns an enum class to specify what the
     * block given by blockName is.
     */
    public NexusBlock findBlockName(String blockName) {
        if (blockName.equalsIgnoreCase(TAXA_BLOCK.toString())) {
            nextBlock = TAXA_BLOCK;
        } else if (blockName.equalsIgnoreCase(CHARACTERS_BLOCK.toString())) {
            nextBlock = CHARACTERS_BLOCK;
        } else if (blockName.equalsIgnoreCase(DATA_BLOCK.toString())) {
            nextBlock = DATA_BLOCK;
        } else if (blockName.equalsIgnoreCase(UNALIGNED_BLOCK.toString())) {
            nextBlock = UNALIGNED_BLOCK;
        } else if (blockName.equalsIgnoreCase(DISTANCES_BLOCK.toString())) {
            nextBlock = DISTANCES_BLOCK;
        } else if (blockName.equalsIgnoreCase(TREES_BLOCK.toString())) {
            nextBlock = TREES_BLOCK;
        } else if (blockName.equalsIgnoreCase(CALIBRATION_BLOCK.toString())) {
            nextBlock = CALIBRATION_BLOCK;
        }

        return nextBlock;
    }

    /**
     * Parses a 'TAXA' block.
     */
    public TaxonList parseTaxaBlock() throws ImportException, IOException {
        return readTaxaBlock();
    }

    /**
     * Parses a 'CHARACTERS' block.
     */
    public Alignment parseCharactersBlock(TaxonList taxonList) throws ImportException, IOException {
        return readCharactersBlock(taxonList);
    }

    /**
     * Parses a 'DATA' block.
     */
    public Alignment parseDataBlock(TaxonList taxonList) throws ImportException, IOException {
        return readDataBlock(/*taxonList*/);
    }

    /**
     * Parses a 'TREES' block.
     */
    public Tree[] parseTreesBlock(TaxonList taxonList) throws ImportException, IOException {
        return readTreesBlock(taxonList, false);
    }

    /**
     * Parses a 'CALIBRATION' block putting the dates into the appropriate taxa
     */
    public dr.evolution.util.Date[] parseCalibrationBlock(TaxonList taxonList) throws ImportException, IOException {
        return readCalibrationBlock(taxonList);
    }

    // **************************************************************
    // SequenceImporter IMPLEMENTATION
    // **************************************************************

    /**
     * importAlignment.
     */
    public Alignment importAlignment() throws IOException, Importer.ImportException {
        boolean done = false;

        TaxonList taxonList = null;
        Alignment alignment = null;

        while (!done) {
            try {

                NexusImporter.NexusBlock block = findNextBlock();

                if (block == NexusImporter.TAXA_BLOCK) {

                    taxonList = readTaxaBlock();

                } else if (block == NexusImporter.CALIBRATION_BLOCK) {
                    if (taxonList == null) {
                        throw new MissingBlockException("TAXA block is missing");
                    }

                    readCalibrationBlock(taxonList);

                } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                    if (taxonList == null) {
                        throw new MissingBlockException("TAXA block is missing");
                    }

                    alignment = readCharactersBlock(taxonList);
                    done = true;

                } else if (block == NexusImporter.DATA_BLOCK) {

                    // A data block doesn't need a taxon block before it
                    // but if one exists then it will use it.
                    alignment = readDataBlock(/*taxonList*/);
                    done = true;

                } else {
                    // Ignore the block..
                }

            } catch (EOFException ex) {
                done = true;
            }
        }

        if (alignment == null) {
            throw new MissingBlockException("DATA or CHARACTERS block is missing");
        }

        return alignment;
    }

    /**
     * importSequences.
     */
    public SequenceList importSequences() throws IOException, ImportException {
        return importAlignment();
    }

    // **************************************************************
    // TreeImporter IMPLEMENTATION
    // **************************************************************

    private boolean isReadingTreesBlock = false;
    private HashMap<String, Taxon> translationList = null;
    private Tree nextTree = null;
    private final String[] lastToken = new String[1];

    /**
     * import a single tree.
     */
    public Tree importTree(TaxonList taxonList) throws IOException, ImportException {
        return importTree(taxonList, false);
    }

    public Tree importTree(TaxonList taxonList, boolean useTaxonListNumbering) throws IOException, ImportException {
        isReadingTreesBlock = false;
        TaxonList[] aTaxonList = new TaxonList[1];
        aTaxonList[0] = taxonList;
        if (!startReadingTrees(aTaxonList)) {
            throw new MissingBlockException("TREES block is missing");
        }
        translationList = readTranslationList(aTaxonList[0], lastToken);
        return readNextTree(translationList, lastToken, useTaxonListNumbering ? taxonList : null);
    }

    /**
     * import an array of all trees.
     */
    public Tree[] importTrees(TaxonList taxonList) throws IOException, ImportException {
        return importTrees(taxonList, false);
    }

    public Tree[] importTrees(TaxonList taxonList, boolean useTaxonListNumbering) throws IOException, ImportException {
        isReadingTreesBlock = false;
        TaxonList[] aTaxonList = new TaxonList[1];
        aTaxonList[0] = taxonList;
        if (!startReadingTrees(aTaxonList)) {
            throw new MissingBlockException("TREES block is missing");
        }
        return readTreesBlock(aTaxonList[0], useTaxonListNumbering);
    }

    /**
     * return whether another tree is available.
     */
    public boolean hasTree() throws IOException, ImportException {
        if (!isReadingTreesBlock) {
            TaxonList[] taxonList = new TaxonList[1];
            taxonList[0] = null;
            isReadingTreesBlock = startReadingTrees(taxonList);
            if (!isReadingTreesBlock) return false;

            translationList = readTranslationList(taxonList[0], lastToken);
        }

        if (nextTree == null) {
            nextTree = readNextTree(translationList, lastToken, null);
        }

        return (nextTree != null);
    }


    /**
     * import the next tree.
     * return the tree or null if no more trees are available
     */
    public Tree importNextTree() throws IOException, ImportException {
        // call hasTree to do the hard work...
        if (!hasTree()) {
            isReadingTreesBlock = false;
            return null;
        }

        Tree tree = nextTree;
        nextTree = null;

        return tree;
    }

    public boolean startReadingTrees(TaxonList[] taxonList) throws IOException, ImportException {
        boolean done = false;

        while (!done) {
            try {

                NexusImporter.NexusBlock block = findNextBlock();

                if (block == NexusImporter.TAXA_BLOCK && taxonList[0] == null) {
                    // only read the taxon list if one hasn't been set already...
                    taxonList[0] = readTaxaBlock();

                } else if (block == NexusImporter.TREES_BLOCK) {

                    return true;
                } else {
                    // Ignore the block..
                }

            } catch (EOFException ex) {
                done = true;
            }
        }

        return false;
    }

    /**
     * Finds the end of the current block.
     */
    private void findToken(String query, boolean ignoreCase) throws IOException {
        String token;
        boolean found = false;

        do {

            token = readToken();

            if ((ignoreCase && token.equalsIgnoreCase(query)) || token.equals(query)) {
                found = true;
            }
        } while (!found);
    }

    /**
     * Finds the end of the current block.
     */
    public void findEndBlock() throws IOException {
        try {
            String token;

            do {
                token = readToken(";");
            } while (!token.equalsIgnoreCase("END") && !token.equalsIgnoreCase("ENDBLOCK"));
        } catch (EOFException e) {
            // Doesn't matter if the End is missing
        }

        nextBlock = UNKNOWN_BLOCK;
    }

    /**
     * Reads the header information for a 'DATA', 'CHARACTERS' or 'TAXA' block.
     */
    private void readDataBlockHeader(String tokenToLookFor, NexusBlock block) throws ImportException, IOException {

        boolean dim = false, ttl = false, fmt = false;
        String token;

        do {
            token = readToken();

            if (token.equalsIgnoreCase("TITLE")) {
                if (ttl) {
                    throw new DuplicateFieldException("TITLE");
                }

                ttl = true;
            } else if (token.equalsIgnoreCase("DIMENSIONS")) {

                if (dim) {
                    throw new DuplicateFieldException("DIMENSIONS");
                }

                boolean nchar = (block == TAXA_BLOCK);
                boolean ntax = (block == CHARACTERS_BLOCK);

                do {
                    String token2 = readToken("=;");

                    if (getLastDelimiter() != '=') {
                        throw new BadFormatException("Unknown subcommand, '" + token2 + "', or missing '=' in DIMENSIONS command");
                    }

                    if (token2.equalsIgnoreCase("NTAX")) {

                        if (block == CHARACTERS_BLOCK) {
                            throw new BadFormatException("NTAX subcommand in CHARACTERS block");
                        }

                        taxonCount = readInteger(";");
                        ntax = true;

                    } else if (token2.equalsIgnoreCase("NCHAR")) {

                        if (block == TAXA_BLOCK) {
                            throw new BadFormatException("NCHAR subcommand in TAXA block");
                        }

                        siteCount = readInteger(";");
                        nchar = true;

                    } else {
                        throw new BadFormatException("Unknown subcommand, '" + token2 + "', in DIMENSIONS command");
                    }

                } while (getLastDelimiter() != ';');

                if (!ntax) {
                    throw new BadFormatException("NTAX subcommand missing from DIMENSIONS command");
                }
                if (!nchar) {
                    throw new BadFormatException("NCHAR subcommand missing from DIMENSIONS command");
                }
                dim = true;

            } else if (token.equalsIgnoreCase("FORMAT")) {

                if (fmt) {
                    throw new DuplicateFieldException("FORMAT");
                }

                dataType = null;

                do {
                    String token2 = readToken("=;");

                    if (token2.equalsIgnoreCase("GAP")) {

                        if (getLastDelimiter() != '=') {
                            throw new BadFormatException("Expecting '=' after GAP subcommand in FORMAT command");
                        }

                        gapCharacters = readToken(";");

                    } else if (token2.equalsIgnoreCase("MISSING")) {

                        if (getLastDelimiter() != '=') {
                            throw new BadFormatException("Expecting '=' after MISSING subcommand in FORMAT command");
                        }

                        missingCharacters = readToken(";");

                    } else if (token2.equalsIgnoreCase("MATCHCHAR")) {

                        if (getLastDelimiter() != '=') {
                            throw new BadFormatException("Expecting '=' after MATCHCHAR subcommand in FORMAT command");
                        }

                        matchCharacters = readToken(";");

                    } else if (token2.equalsIgnoreCase("DATATYPE")) {

                        if (getLastDelimiter() != '=') {
                            throw new BadFormatException("Expecting '=' after DATATYPE subcommand in FORMAT command");
                        }

                        String token3 = readToken(";");
                        if (token3.equalsIgnoreCase("NUCLEOTIDE") ||
                                token3.equalsIgnoreCase("DNA") ||
                                token3.equalsIgnoreCase("RNA")) {

                            dataType = Nucleotides.INSTANCE;

                        } else if (token3.equalsIgnoreCase("STANDARD") || token3.equalsIgnoreCase("BINARY")) {

                            dataType = TwoStates.INSTANCE;

                        } else if (token3.equalsIgnoreCase("PROTEIN")) {

                            dataType = AminoAcids.INSTANCE;

                        } else if (token3.equalsIgnoreCase("CONTINUOUS")) {

                            throw new UnparsableDataException("Continuous data cannot be parsed at present");

                        }
                    } else if (token2.equalsIgnoreCase("INTERLEAVE")) {
                        isInterleaved = true;
                    }

                } while (getLastDelimiter() != ';');

                fmt = true;
            }
        } while (!token.equalsIgnoreCase(tokenToLookFor));

        if (!dim) {
            throw new MissingFieldException("DIMENSIONS");
        }
        if (block != TAXA_BLOCK && dataType == null) {
            throw new MissingFieldException("DATATYPE");
        }
    }

    /**
     * Reads sequences in a 'DATA' or 'CHARACTERS' block.
     */
    private void readSequenceData(Sequences sequences, TaxonList taxonList) throws ImportException, IOException {
        int n, i;
        String firstSequence = null;

        if (isInterleaved) {
            boolean firstLoop = true;

            int readCount = 0;
            while (readCount < siteCount) {

                n = -1;

                for (i = 0; i < taxonCount; i++) {

                    String token = readToken().trim();

                    Sequence sequence;

                    if (firstLoop) {

                        sequence = new Sequence();
                        sequence.setDataType(dataType);
                        sequences.addSequence(sequence);

                        Taxon taxon;

                        if (taxonList != null) {
                            int index = taxonList.getTaxonIndex(token.trim());
                            if (index == -1) {
                                // taxon not found in taxon list...
                                // ...perhaps it is a numerical taxon reference?
                                throw new UnknownTaxonException(token);
                            } else {
                                taxon = taxonList.getTaxon(index);
                            }
                        } else {
                            taxon = new Taxon(token.trim());
                        }

                        sequence.setTaxon(taxon);

                    } else {

                        sequence = sequences.getSequence(i);
                        Taxon taxon = sequence.getTaxon();
                        if (!taxon.getId().equals(token)) {
                            throw new UnknownTaxonException("Unknown taxon label: expecting '" +
                                    taxon.getId() + "', found '" + token + "'");
                        }
                    }

                    StringBuffer buffer = new StringBuffer();
                    readSequenceLine(buffer, dataType, ";", gapCharacters, missingCharacters,
                            matchCharacters, firstSequence);
                    String seqString = buffer.toString();
                    sequence.appendSequenceString(seqString);
                    if (i == 0) {
                        firstSequence = seqString;
                    }

                    if (getLastDelimiter() == ';') {
                        if (i < taxonCount - 1) {
                            throw new TooFewTaxaException();
                        }
                        if (readCount + n < siteCount) {
                            throw new ShortSequenceException(sequence.getTaxon().getId());
                        }
                    }

                    if (n == -1) {
                        n = seqString.length();
                    }

                    if (n != seqString.length()) {
                        throw new ShortSequenceException(sequence.getTaxon().getId());
                    }
                }

                firstLoop = false;
                readCount += n;

            }
            if (getLastDelimiter() != ';') {
                throw new BadFormatException("Expecting ';' after sequences data");
            }

        } else {

            for (i = 0; i < taxonCount; i++) {
                String token = readToken().trim();

                Sequence sequence = new Sequence();
                sequence.setDataType(dataType);
                sequences.addSequence(sequence);

                Taxon taxon;

                if (taxonList != null) {
                    int index = taxonList.getTaxonIndex(token);
                    if (index == -1) {
                        // taxon not found in taxon list...
                        // ...perhaps it is a numerical taxon reference?
                        throw new UnknownTaxonException(token);
                    } else {
                        taxon = taxonList.getTaxon(index);
                    }
                } else {
                    taxon = new Taxon(token);
                }

                sequence.setTaxon(taxon);

                StringBuffer buffer = new StringBuffer();
                readSequence(buffer, dataType, ";", siteCount, gapCharacters,
                        missingCharacters, matchCharacters, firstSequence);
                String seqString = buffer.toString();
                if (seqString.length() != siteCount) {
                    throw new ShortSequenceException(sequence.getTaxon().getId());
                }

                sequence.appendSequenceString(seqString);
                if (i == 0) {
                    firstSequence = seqString;
                }

                if (getLastDelimiter() == ';' && i < taxonCount - 1) {
                    throw new TooFewTaxaException();
                }

            }

            if (getLastDelimiter() != ';') {
                throw new BadFormatException("Expecting ';' after sequences data, has '"
                        + (char) getLastDelimiter() + "' in line " + getLineNumber());
            }

        }
    }


    /**
     * Reads a 'TAXA' block.
     */
    private TaxonList readTaxaBlock() throws ImportException, IOException, IllegalArgumentException {

        taxonCount = 0;

        readDataBlockHeader("TAXLABELS", TAXA_BLOCK);

        if (taxonCount == 0) {
            throw new MissingFieldException("NTAXA");
        }

        Taxa taxa = new Taxa();

        do {
            String name = readToken(";").trim();
            if (name.length() > 0) {
                Taxon taxon = new Taxon(name);
                taxa.addTaxon(taxon);
            }
        } while (getLastDelimiter() != ';');

        if (taxa.getTaxonCount() != taxonCount) {
            throw new BadFormatException("Number of taxa doesn't match NTAXA field");
        }

        findEndBlock();

        int duplicateTaxon = TaxonList.Utils.findDuplicateTaxon(taxa);
        if (duplicateTaxon >= 0)
            throw new IllegalArgumentException("Tree contains duplicate taxon name: " + taxa.getTaxon(duplicateTaxon).getId() +
                    "!\nAll taxon names should be unique.");

        return taxa;
    }

    /**
     * Reads a 'CHARACTERS' block.
     */
    private Alignment readCharactersBlock(TaxonList taxonList) throws ImportException, IOException {

        siteCount = 0;
        dataType = null;

        readDataBlockHeader("MATRIX", CHARACTERS_BLOCK);

        SimpleAlignment alignment = new SimpleAlignment();
        readSequenceData(alignment, taxonList);
        alignment.updateSiteCount();

        findEndBlock();

        return alignment;
    }

    /**
     * Reads a 'DATA' block.
     */
    private Alignment readDataBlock(/*TaxonList taxonList*/) throws ImportException, IOException {

        taxonCount = 0;
        siteCount = 0;
        dataType = null;

        readDataBlockHeader("MATRIX", DATA_BLOCK);

        SimpleAlignment alignment = new SimpleAlignment();
        readSequenceData(alignment, null);
        alignment.updateSiteCount();

        findEndBlock();

        return alignment;
    }


    /**
     * Reads a 'TREES' block.
     */
    private Tree[] readTreesBlock(TaxonList taxonList, boolean useTaxonListNumbering) throws ImportException, IOException {
        ArrayList<Tree> trees = new ArrayList<Tree>();

        String[] lastToken = new String[1];
        HashMap<String, Taxon> translationList = readTranslationList(taxonList, lastToken);

        boolean done = false;
        do {

            Tree tree = readNextTree(translationList, lastToken, useTaxonListNumbering ? taxonList : null);

            if (tree != null) {
                trees.add(tree);
            } else {
                done = true;
            }
        } while (!done);

        if (trees.size() == 0) {
            throw new BadFormatException("No trees defined in TREES block");
        }

        Tree[] treeArray = new Tree[trees.size()];
        trees.toArray(treeArray);

        nextBlock = UNKNOWN_BLOCK;

        return treeArray;
    }

    private HashMap<String, Taxon> readTranslationList(TaxonList taxonList, String[] lastToken) throws ImportException, IOException {
        HashMap<String, Taxon> translationList = new HashMap<String, Taxon>();

        String token = readToken(";");

        if (token.equalsIgnoreCase("TRANSLATE")) {

            do {
                String token2 = readToken(",;");

                if (getLastDelimiter() == ',' || getLastDelimiter() == ';') {
                    throw new BadFormatException("Missing taxon label in TRANSLATE command of TREES block");
                }

                String token3 = readToken(",;");
                Taxon taxon;

                if (getLastDelimiter() != ',' && getLastDelimiter() != ';') {
                    throw new BadFormatException("Expecting ',' or ';' after taxon label in TRANSLATE command of TREES block");
                }

                if (taxonList != null) {
                    int index = taxonList.getTaxonIndex(token3);
                    if (index == -1) {
                        // taxon not found in taxon list...
                        // ...perhaps it is a numerical taxon reference?
                        throw new UnknownTaxonException(token3);
                    } else {
                        taxon = taxonList.getTaxon(index);
                    }
                } else {
                    taxon = new Taxon(token3);
                }

                if (translationList.containsKey(token2)) {
                    throw new BadFormatException("Translation list uses the key, " + token2 + ", more than once.");
                }
                translationList.put(token2, taxon);

            } while (getLastDelimiter() != ';');

            token = readToken(";");
        } else if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                translationList.put(taxon.getId(), taxon);
            }
        }

        lastToken[0] = token;

        return translationList;
    }

    private Tree readNextTree(HashMap<String, Taxon> translationList, String[] lastToken, TaxonList taxonList) throws ImportException, IOException {
        try {
            Tree tree = null;
            String token = lastToken[0];

            if (token.equalsIgnoreCase("UTREE") || token.equalsIgnoreCase("TREE")) {

                if (nextCharacter() == '*') {
                    // Star is used to specify a default tree - ignore it
                    readCharacter();
                }

                String token2 = readToken("=;");
                // Save tree comment and attach it later
                final String comment = getLastMetaComment();
                clearLastMetaComment();

                if (getLastDelimiter() != '=') {
                    throw new BadFormatException("Missing label for tree'" + token2 + "' or missing '=' in TREE command of TREES block");
                }

                try {
                    if (nextCharacter() != '(') {
                        throw new BadFormatException("Missing tree definition in TREE command of TREES block");
                    }

                    // tree special comments
                    final String scomment = getLastMetaComment();
                    clearLastMetaComment();

                    FlexibleNode root = readInternalNode(translationList);

                    if (translationList != null) {
                        // this ensures that if a translation list is used, the external node numbers
                        // of the trees correspond as well.

                        Map<Taxon, Integer> taxonNumberMap = new HashMap<Taxon, Integer>();
                        int count = 0;
                        for (String label : translationList.keySet()) {
                            Taxon taxon = translationList.get(label);
                            int number;

                            if (taxonList != null) { // Map back to original numbering from TaxonList
                                number =  taxonList.getTaxonIndex(taxon);
                            } else { // Old functionality
                                try {
                                    number = Integer.parseInt(label) - 1;
                                } catch (NumberFormatException nfe) {
                                    number = count;
                                }
                            }

                            taxonNumberMap.put(taxon, number);
                            count++;
                        }

                        tree = new FlexibleTree(root, false, true, taxonNumberMap);
                    } else {
                        tree = new FlexibleTree(root, false, true, null);
                    }

                    tree.setId(token2);

                    if (getLastDelimiter() == ':') {
                        // in case the root has a branch length, skip it
                        readToken(";");

                        if (getLastMetaComment() != null) {
                            // There was a meta-comment which should be in the form:
                            // \[&label[=value][,label[=value]>[,/..]]\]
                            try {
                                parseMetaCommentPairs(getLastMetaComment(), root);
                            } catch (BadFormatException bfe) {
                                // ignore it
                            }
                            clearLastMetaComment();
                        }
                    }

                    if (getLastDelimiter() != ';') {
                        throw new BadFormatException("Expecting ';' after tree, '" + token2 + "', TREE command of TREES block");
                    }

                    if (scomment != null) {
                        // below is correct only if [&W] appears on it own
                        String c = scomment;
                        while (c.length() > 0) {
                            final char ch = c.charAt(0);
                            if (ch == ';') {
                                c = c.substring(1);
                                continue;
                            }
                            if (ch == 'R') {
                                // we only have rooted trees anyway
                                c = c.substring(1);
                            } else if (ch == 'W') {
                                int e = c.indexOf(';');
                                if (e < 0) e = c.length();

                                try {
                                    final Float value = new Float(c.substring(2, e));
                                    tree.setAttribute("weight", value);
                                } catch (NumberFormatException ex) {
                                    // don't fail, ignore
                                }
                                c = c.substring(e);
                            } else {
                                c = c.substring(1);
                            }
                        }
                    }

                    if (comment != null) {
                        try {
                            parseMetaCommentPairs(comment, tree);
                        } catch (Importer.BadFormatException e) {
                            // set generic comment attribute
                            tree.setAttribute("comment", comment);
                        }
                    }

                } catch (EOFException e) {
                    // If we reach EOF we may as well return what we have?
                    return tree;
                }

                token = readToken(";");
            } else if (token.equalsIgnoreCase("ENDBLOCK") || token.equalsIgnoreCase("END")) {
                return null;
            } else {
                throw new BadFormatException("Unknown command '" + token + "' in TREES block");
            }

            //added this to escape readNextTree loop correctly -- AJD
            lastToken[0] = token;

            return tree;

        } catch (EOFException e) {
            return null;
        }
    }

    /**
     * Reads a branch in. This could be a node or a tip (calls readNode or readTip
     * accordingly). It then reads the branch length and SimpleNode that will
     * point at the new node or tip.
     */
    FlexibleNode readBranch(HashMap<String, Taxon> translationList) throws IOException, ImportException {
        double length = 0.0;
        FlexibleNode branch;

        clearLastMetaComment();

        if (nextCharacter() == '(') {
            // is an internal node
            branch = readInternalNode(translationList);

        } else {
            // is an external node
            branch = readExternalNode(translationList);
        }

        if (getLastDelimiter() != ':' && getLastDelimiter() != ',' && getLastDelimiter() != ')') {
            String label = readToken(",():;");
            if (label.length() > 0) {
                branch.setAttribute("label", label);
            }
        }

        if (getLastDelimiter() == ':') {
            length = readDouble(",():;");

            if (getLastMetaComment() != null) {
                if (!ignoreMetaComments) {
                    // There was a meta-comment which should be in the form:
                    // \[&label[=value][,label[=value]>[,/..]]\]
                    try {
                        parseMetaCommentPairs(getLastMetaComment(), branch);
                    } catch (BadFormatException bfe) {
                        // ignore it
                    }
                }
                clearLastMetaComment();
            }

        }

        branch.setLength(length);

        return branch;
    }

    /**
     * Reads a node in. This could be a polytomy. Calls readBranch on each branch
     * in the node.
     */
    FlexibleNode readInternalNode(HashMap<String, Taxon> translationList) throws IOException, ImportException {
        FlexibleNode node = new FlexibleNode();

        // read the opening '('
        readCharacter();

        // read the first child
        node.addChild(readBranch(translationList));

        if (getLastDelimiter() != ',' && !suppressWarnings) {
            java.util.logging.Logger.getLogger("dr.evolution.io").warning("Internal node only has a single child.");
        }

        // this allows one or more children
        while(getLastDelimiter()==',') {
            node.addChild(readBranch(translationList));
        }

        // should have had a closing ')'
        if (getLastDelimiter() != ')') {
            throw new BadFormatException("Missing closing ')' in tree in TREES block");
        }

        readToken(":(),;");

        if (getLastMetaComment() != null) {
            if (!ignoreMetaComments) {
                // There was a meta-comment which should be in the form:
                // \[&label[=value][,label[=value]>[,/..]]\]
                try {
                    parseMetaCommentPairs(getLastMetaComment(), node);
                } catch (BadFormatException bfe) {
                    // ignore it
                }
            }
            clearLastMetaComment();
        }

        // find the next delimiter
        return node;
    }

//	private void labelNode(FlexibleNode node, String label, String value) {
//		// Attempt to format the value as a number
//		Number number = null;
//		try {
//			number = Integer.valueOf(value);
//		} catch (NumberFormatException nfe1) {
//			try {
//				number = Double.valueOf(value);
//			} catch (NumberFormatException nfe2) {
//				//
//			}
//		}
//		if (number != null) {
//			node.setAttribute(label, number);
//		} else {
//			node.setAttribute(label, value);
//		}
//	}

    /**
     * Reads an external node in.
     */
    FlexibleNode readExternalNode(HashMap<String, Taxon> translationList) throws ImportException, IOException {
        FlexibleNode node = new FlexibleNode();

        String label = readToken(":(),;");

        Taxon taxon;

        if (translationList.size() > 0) {
            taxon = translationList.get(label);

            if (taxon == null) {
                // taxon not found in taxon list...
                throw new UnknownTaxonException("Taxon in tree, '" + label + "' is unknown");
            }
        } else {
            taxon = new Taxon(label);
        }

        if (getLastMetaComment() != null) {
            if (!ignoreMetaComments) {
                // There was a meta-comment which should be in the form:
                // \[&label[=value][,label[=value]>[,/..]]\]
                try {
                    parseMetaCommentPairs(getLastMetaComment(), node);
                } catch (BadFormatException bfe) {
                    // ignore it
                }
            }
            clearLastMetaComment();
        }

        node.setTaxon(taxon);
        return node;
    }

    /**
     * Reads a 'CALIBRATION' block.
     */
    private dr.evolution.util.Date[] readCalibrationBlock(TaxonList taxonList) throws ImportException, IOException {
        double origin = 0.0;
        boolean isBackwards = false;
        Units.Type units = Units.Type.YEARS;
        ArrayList<Date> dates = new ArrayList<Date>();

        String token;

        boolean done = false;
        do {
            token = readToken(";");

            if (token.equalsIgnoreCase("OPTIONS")) {

                do {
                    String token2 = readToken("=;");

                    if (getLastDelimiter() != '=') {
                        throw new BadFormatException("Unknown subcommand, '" + token2 + "', or missing '=' in OPTIONS command of CALIBRATION block");
                    }

                    if (token2.equalsIgnoreCase("SCALE")) {

                        String token3 = readToken(";");
                        if (token3.equalsIgnoreCase("DAYS")) {

                            units = Units.Type.DAYS;

                        } else if (token3.equalsIgnoreCase("MONTHS")) {

                            units = Units.Type.MONTHS;

                        } else if (token3.equalsIgnoreCase("YEARS")) {

                            units = Units.Type.YEARS;

                        } else {
                            throw new BadFormatException("SCALE in OPTIONS command of CALIBRATION block must be one of DAYS, MONTHS or YEARS");
                        }

                    } else if (token2.equalsIgnoreCase("ORIGIN")) {

                        origin = readDouble(";");

                    } else if (token2.equalsIgnoreCase("DIRECTION")) {

                        String token3 = readToken(";");
                        if (token3.equalsIgnoreCase("FORWARDS")) {

                            isBackwards = false;

                        } else if (token3.equalsIgnoreCase("BACKWARDS")) {

                            isBackwards = true;

                        } else {
                            throw new BadFormatException("DIRECTION in OPTIONS command of CALIBRATION block must be either FORWARDS or BACKWARDS");
                        }


                    } else {
                        throw new BadFormatException("Unknown subcommand, '" + token2 + "', in OPTIONS command of CALIBRATION block");
                    }

                } while (getLastDelimiter() != ';');

            } else if (token.equalsIgnoreCase("TIPCALIBRATION")) {

                do {
                    String token2 = readToken("=;");

                    if (getLastDelimiter() != '=') {
                        throw new BadFormatException("Missing date for label '" + token2 + "' or missing '=' in TIPCALIBRATION command of CALIBRATION block");
                    }

                    double value = readDouble(":;");

                    if (getLastDelimiter() != ':') {
                        throw new BadFormatException("Missing taxon list for label '" + token2 + "' or missing ':' in TIPCALIBRATION command of CALIBRATION block");
                    }

                    dr.evolution.util.Date date;
                    if (isBackwards) {
                        date = dr.evolution.util.Date.createTimeAgoFromOrigin(value, units, origin);
                    } else {
                        date = dr.evolution.util.Date.createTimeSinceOrigin(value, units, origin);
                    }

                    dates.add(date);

                    do {
                        String token3 = readToken(",;");
                        Taxon taxon;

                        int index = taxonList.getTaxonIndex(token3);
                        if (index == -1) {
                            // taxon not found in taxon list...
                            throw new UnknownTaxonException("Unknown taxon '" + token3 + "' for label '" + token2 + "' in TIPCALIBRATION command of CALIBRATION block");
                        } else {
                            taxon = taxonList.getTaxon(index);
                        }

                        taxon.setAttribute("date", date);

                    } while (getLastDelimiter() != ',' && getLastDelimiter() != ';');

                } while (getLastDelimiter() == ',');

            } else if (token.equalsIgnoreCase("NODECALIBRATION")) {
                throw new BadFormatException("NODECALIBRATION not suppored in CALIBRATION block");
            } else if (token.equalsIgnoreCase("ENDBLOCK") || token.equalsIgnoreCase("END")) {
                done = true;
            } else {
                throw new BadFormatException("Unknown command '" + token + "' in CALIBRATION block");
            }
        } while (!done);

        dr.evolution.util.Date[] dateArray = new dr.evolution.util.Date[dates.size()];
        dates.toArray(dateArray);

        nextBlock = UNKNOWN_BLOCK;

        return dateArray;
    }

    static void parseMetaCommentPairs(String meta, Attributable item) throws Importer.BadFormatException {
        if (meta.startsWith("B ")) {
            // a MrBayes annotation
            String[] parts = meta.split(" ");
            if (parts.length == 3 && parts[1].length() > 0 && parts[2].length() > 0) {
                item.setAttribute(parts[1], parseValue(parts[2]));
            } else if (parts.length == 2 && parts[1].length() > 0) {
                item.setAttribute(parts[1], Boolean.TRUE);
            } else {
                throw new Importer.BadFormatException("Badly formatted attribute: '" + meta + "'");
            }
            return;
        }

        // This regex should match key=value pairs, separated by commas
        // This can match the following types of meta comment pairs:
        // value=number, value="string", value={item1, item2, item3}
        // (label must be quoted if it contains spaces (i.e. "my label"=label)

        // TODO MAS Minor change in line below for nested arrays may cause other unforeseen bugs
        Pattern pattern = Pattern.compile("(\"[^\"]*\"+|[^,=\\s]+)\\s*(=\\s*(\\{[^=]*\\}|\"[^\"]*\"+|[^,]+))?");
        Matcher matcher = pattern.matcher(meta);

        while (matcher.find()) {
            String label = matcher.group(1);
            if (label.charAt(0) == '\"') {
                label = label.substring(1, label.length() - 1);
            }
            if (label == null || label.trim().length() == 0) {
                throw new Importer.BadFormatException("Badly formatted attribute: '" + matcher.group() + "'");
            }
            final String value = matcher.group(2);
            if (value != null && value.trim().length() > 0) {
                // there is a specified value so try to parse it
                item.setAttribute(label, parseValue(value.substring(1)));
            } else {
                item.setAttribute(label, Boolean.TRUE);
            }
        }
    }

    public static boolean isInt(String str)
    {
        if (str == null) {
            return false;
        }
        int length = str.length();
        if (length == 0) {
            return false;
        }
        int i = 0;
        if (str.charAt(0) == '-') {
            if (length == 1) {
                return false;
            }
            i = 1;
        }
        for (; i < length; i++) {
            char c = str.charAt(i);
            if (c <= '/' || c >= ':') {
                return false;
            }
        }
        return true;
    }

    /**
     * This method takes a string and tries to decode it returning the object
     * that best fits the data. It will recognize comma delimited lists enclosed
     * in {..} and call parseValue() on each element. It will also recognize Boolean,
     * Integer and Double. If the value starts with a # then it will attempt to decode
     * the following integer as an RGB colour - see Color.decode(). If nothing else fits
     * then the value will be returned as a string but trimmed of leading and trailing
     * white space.
     *
     * @param value the string
     * @return the object
     */
    public static Serializable parseValue(String value) {

        value = value.trim();

        if (value.startsWith("{")) {
            // the value is a list so recursively parse the elements
            // and return an array
            String inside = value.substring(1, value.length() - 1);

            if (inside.length() == 0) {
                return null;
            }

             // Determine depth of further nesting
            int depth = 0;
            while (inside.charAt(depth) == '{') {
                depth++;
            }

            StringBuilder split;
            if (depth == 0) {
                split = new StringBuilder(",");
            } else {
                StringBuilder rightBookEnd = new StringBuilder("(?<=");
                StringBuilder leftBookEnd = new StringBuilder("(?=");
                for (int i = 0; i < depth; ++i) {
                    rightBookEnd.append("\\}");
                    leftBookEnd.append("\\{");
                }
                leftBookEnd.append(")");
                rightBookEnd.append(")");
                split = rightBookEnd.append(",").append(leftBookEnd);
            }

            // Non-destructive split
            String[] elements = inside.split(split.toString());
            Object[] values = new Object[elements.length];
            for (int i = 0; i < elements.length; i++) {
                values[i] = parseValue(elements[i]);
            }
            return values;
        }

        if (value.startsWith("#")) {
            // I am not sure whether this is a good idea but
            // I am going to assume that a # denotes an RGB colour
            try {
                return Color.decode(value.substring(1));
            } catch (NumberFormatException nfe1) {
                // not a colour
            }
        }

        if (value.equalsIgnoreCase("TRUE") || value.equalsIgnoreCase("FALSE")) {
            return Boolean.valueOf(value);
        }

        // Attempt to format the value as an integer
//        try {
//            return new Integer(value);
//        } catch (NumberFormatException nfe1) {
//            // not an integer
//        }
        // throwing exception to test for an integer is slow
        if (isInt(value)) {
            return Integer.valueOf(value);
        }

        // Attempt to format the value as a double
        try {
            return new Double(value);
        } catch (NumberFormatException nfe2) {
            // not a double
        }

        // return the trimmed string
        return value;
    }

    // private stuff

    private NexusBlock nextBlock = null;

    private int taxonCount = 0, siteCount = 0;
    private DataType dataType = null;
    private String gapCharacters = "-";
    private String matchCharacters = ".";
    private String missingCharacters = "?";
    private boolean isInterleaved = false;

    public static void main(String[] args) throws IOException, ImportException {

        if (args.length > 3) {

            int sampleFrequency = Integer.parseInt(args[1]);
            boolean includeBranchLengths = Boolean.getBoolean(args[2]);
            boolean isNexus = Boolean.getBoolean(args[3]);

            NexusImporter nexusImporter = null;

            BufferedReader reader = null;

            if (isNexus) {
                nexusImporter = new NexusImporter(new FileReader(args[0]));
            } else {
                reader = new BufferedReader(new FileReader(args[0]));
            }


            int index = 0;
            int count = 0;
            String line = null;
            if (!isNexus) line = reader.readLine();
            while (line != null || (isNexus && nexusImporter.hasTree())) {

                Tree tree;
                if (isNexus) {
                    tree = nexusImporter.importNextTree();
                } else {
                    String treeString = line.substring(line.indexOf('(')).trim();

                    java.io.Reader stringReader = new java.io.StringReader(treeString);
                    NewickImporter importer = new NewickImporter(stringReader);
                    tree = importer.importNextTree();
                }


                if (index % sampleFrequency == 0) {
                    if (includeBranchLengths) {
                        System.out.println(TreeUtils.newick(tree));
                    } else {
                        System.out.println(TreeUtils.newickNoLengths(tree));
                        count += 1;
                    }
                }
                index += 1;
                if (!isNexus) line = reader.readLine();
            }
            System.out.println(count + " trees");

        } else {
            System.err.println("usage: filterTrees <tree-file-name> <sample-frequency> <include-branch-lengths>");
        }


    }

}
