/*
 * BEAUTiImporter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.util;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.datapanel.BadTraitFormatDialog;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.options.*;
import dr.app.util.Utils;
import dr.evolution.alignment.*;
import dr.evolution.datatype.*;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.io.NexusImporter.NexusBlock;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.util.DataTable;
import org.jdom.JDOMException;

import javax.swing.*;
import java.io.*;
import java.util.*;
import java.util.List;


/**
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class BEAUTiImporter {

    private static final long MAX_FILE_SIZE = 1024 * 1024 * 128; // 0.1GB
    private final String DEFAULT_NAME = "default";

    private final BeautiOptions options;
    private final BeautiFrame frame;

    private PartitionNameDialog partitionNameDialog = null;

    public BEAUTiImporter(BeautiFrame frame, BeautiOptions options) {
        this.frame = frame;
        this.options = options;
    }

    public void importFromFile(File file) throws IOException, ImportException, JDOMException {
        long fileSize = file.length();
        if (fileSize > MAX_FILE_SIZE) {
            String size = String.format("%.2f", ((double)fileSize) / (1024 * 1024 * 1024));
            int result = JOptionPane.showConfirmDialog(this.frame, "This file is " + size +"GB in size.\nAre you sure you want to import it?",
                    "Importing Data", JOptionPane.OK_CANCEL_OPTION);
            if (result == JOptionPane.CANCEL_OPTION) {
                return;
            }
        }

        try {
            String line = findFirstLine(file).trim();

            if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
                // is a NEXUS file
                importNexusFile(file);
            } else if ((line != null && line.startsWith("" + FastaImporter.FASTA_FIRST_CHAR))) {
                // is a FASTA file
                importFastaFile(file);
            } else if ((line != null && (line.toUpperCase().contains("<?XML") || line.toUpperCase().contains("<BEAST")))) {
                // assume it is a BEAST XML file and see if that works...
                importBEASTFile(file);
            } else if ((line != null && line.startsWith("("))) {
                importNewickFile(file);
            } else {
                // assume it is a tab-delimited traits file and see if that works...
                importTraits(file);
//            } else {
//                throw new ImportException("Unrecognized format for imported file.");
            }

        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    public String findFirstLine(File file) throws IOException {
        Reader reader = new FileReader(file);

        BufferedReader bufferedReader = new BufferedReader(reader);
        String line = bufferedReader.readLine();
        while (line != null && line.isEmpty()) {
            line = bufferedReader.readLine();
        }

        bufferedReader.close();
        return line;
    }

    // xml
    private void importBEASTFile(File file) throws IOException, ImportException, JDOMException {
        try {
            FileReader reader = new FileReader(file);

            BeastImporter importer = new BeastImporter(reader);

            List<TaxonList> taxonLists = new ArrayList<TaxonList>();
            List<Alignment> alignments = new ArrayList<Alignment>();

            importer.importBEAST(taxonLists, alignments);

            TaxonList taxa = taxonLists.get(0);

            int count = 1;
            for (Alignment alignment : alignments) {
                String name = file.getName();
                if (alignment.getId() != null && alignment.getId().length() > 0) {
                    name = alignment.getId();
                } else {
                    if (alignments.size() > 1) {
                        name += count;
                    }
                }
                setData(name, taxa, alignment, null, null, null, null, null, 0, false);

                count++;
            }

            // assume that any additional taxon lists are taxon sets...
            for (int i = 1; i < taxonLists.size(); i++) {
                Taxa taxonSet = (Taxa) taxonLists.get(i);

                options.taxonSets.add(taxonSet);
                options.taxonSetsMono.put(taxonSet, false);
                options.taxonSetsIncludeStem.put(taxonSet, false);
                options.taxonSetsTreeModel.put(taxonSet, options.getPartitionTreeModels().get(0));
                options.updateTMRCAStatistic(taxonSet); // force the creation of a tmrca statistic
            }
            Collections.sort(options.taxonSets);

            reader.close();
        } catch (JDOMException e) {
            throw new JDOMException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }

    }

    public void importNexusFile(File file) throws IOException, ImportException {
        importNexusFile(file, false);
    }

    // nexus
    public void importNexusFile(File file, Boolean allowEmpty) throws IOException, ImportException {
        TaxonList taxa = null;
        SimpleAlignment alignment = null;
        List<Tree> trees = new ArrayList<Tree>();
        PartitionSubstitutionModel model = null;
        List<CharSet> charSets = new ArrayList<CharSet>();
        List<NexusApplicationImporter.TaxSet> taxSets = new ArrayList<NexusApplicationImporter.TaxSet>();

        int treeCount = 0;

        try {
            NexusApplicationImporter importer = new NexusApplicationImporter(new FileReader(file));

            boolean done = false;
            while (!done) {
                try {
                    NexusBlock block = importer.findNextBlock();
                    if (block == NexusImporter.TREES_BLOCK) {
                        treeCount += importer.countTrees();
                    }
                } catch (EOFException ex) {
                    done = true;
                }
            }

            boolean importAllTrees = true;
            if (treeCount > 10) {
                int result = JOptionPane.showConfirmDialog(this.frame,
                        "This file contains " + treeCount + " trees.\n" +
                                "Do you want to just import the first tree in the file (recommended)?\n" +
                                "Choose 'No' to import all the trees.",
                        "Importing Data", JOptionPane.YES_NO_CANCEL_OPTION);
                if (result == JOptionPane.CANCEL_OPTION) {
                    return;
                }
                importAllTrees = result == JOptionPane.NO_OPTION;
            }

            // reset file for reading
            importer = new NexusApplicationImporter(new FileReader(file));

            done = false;
            while (!done) {
                try {

                    NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxa != null) {
                            throw new MissingBlockException("TAXA block already defined");
                        }

                        taxa = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.CALIBRATION_BLOCK) {
                        if (taxa == null) {
                            throw new MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
                        }

                        importer.parseCalibrationBlock(taxa);

                    } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                        if (taxa == null) {
                            throw new MissingBlockException("TAXA block must be defined before a CHARACTERS block");
                        }

                        if (alignment != null) {
                            throw new MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        alignment = (SimpleAlignment) importer.parseCharactersBlock(taxa);

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        if (alignment != null) {
                            throw new MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = (SimpleAlignment) importer.parseDataBlock(taxa);
                        if (taxa == null) {
                            taxa = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {
                        if (!importAllTrees) {
                            if (trees.isEmpty()) {
                                trees.addAll(importer.parseTreesBlock(taxa, 1));
                            }
                        } else {
                            trees.addAll(importer.parseTreesBlock(taxa));
                        }
                        if (taxa == null && !trees.isEmpty()) {
                            taxa = trees.get(0);
                        }
                    } else if (block == NexusApplicationImporter.PAUP_BLOCK) {

                        model = importer.parsePAUPBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

                        model = importer.parseMrBayesBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.ASSUMPTIONS_BLOCK || block == NexusApplicationImporter.SETS_BLOCK) {

                        importer.parseAssumptionsBlock(charSets, taxSets);

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (alignment == null && taxa == null) {
                throw new MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (IOException e) {
            throw new IOException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
//        } catch (Exception e) {
//            throw new Exception(e.getMessage());
        }

        setData(file.getName(), taxa, alignment, charSets, taxSets, model, null, trees, treeCount, allowEmpty);
    }

    // FASTA

    private void importFastaFile(File file) throws IOException, ImportException {
        try {
            FileReader reader = new FileReader(file);

            FastaImporter importer = new FastaImporter(reader, Nucleotides.INSTANCE);

            Alignment alignment = importer.importAlignment();

            reader.close();

            setData(file.getName(), alignment, alignment, null, null, null, null, null, 0, false);
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    private boolean isMissingValue(String value) {
        return (value.equals("?") || value.equals("NA") || value.length() == 0);
    }

    public void importTraits(final File file) throws IOException {
        List<TraitData> importedTraits = new ArrayList<TraitData>();

        if (options.taxonList == null) {
            options.taxonList = new Taxa();
        }
        Taxa taxa = options.taxonList;

        DataTable<String[]> dataTable = DataTable.Text.parse(new FileReader(file));

        String[] traitNames = dataTable.getColumnLabels();
        String[] taxonNames = dataTable.getRowLabels();

        for (int i = 0; i < dataTable.getColumnCount(); i++) {
            boolean warningGiven = false;

            String traitName = traitNames[i];

            String[] values = dataTable.getColumn(i);
            Class c = null;
            if (!isMissingValue(values[0])) {
                c = Utils.detectType(values[0]);
            }
            for (int j = 1; j < values.length; j++) {
                if (!isMissingValue(values[j])) {
                    if (c == null) {
                        c = Utils.detectType(values[j]);
                    } else {
                        Class c1 = Utils.detectType(values[j]);
                        if (c == Integer.class && c1 == Double.class) {
                            // change the type to double
                            c = Double.class;
                        }

                        if (c1 != c &&
                                !(c == Double.class && c1 == Integer.class) &&
                                !warningGiven) {
                            JOptionPane.showMessageDialog(frame, "Not all values of same type for trait" + traitName,
                                    "Incompatible values", JOptionPane.WARNING_MESSAGE);
                            warningGiven = true;
                        }
                    }
                }
            }

            TraitData.TraitType t = (c == Boolean.class || c == String.class || c == null) ? TraitData.TraitType.DISCRETE :
                    (c == Integer.class) ? TraitData.TraitType.INTEGER : TraitData.TraitType.CONTINUOUS;
            TraitData newTrait = new TraitData(options, traitName, file.getName(), t);

            if (frame.validateTraitName(traitName)) {
                importedTraits.add(newTrait);

                int j = 0;
                for (final String taxonName : taxonNames) {
                    final int index = taxa.getTaxonIndex(taxonName);
                    Taxon taxon;
                    if (index >= 0) {
                        taxon = taxa.getTaxon(index);
                    } else {
                        taxon = new Taxon(taxonName);
                        taxa.addTaxon(taxon);
                    }
                    if (!isMissingValue(values[j])) {
                        taxon.setAttribute(traitName, Utils.constructFromString(c, values[j]));
                    } else {
                        // AR - merge rather than replace existing trait values
                        if (taxon.getAttribute(traitName) == null) {
                            taxon.setAttribute(traitName, "?");
                        }
                    }
                    j++;
                }
            }
        }
        try {
            setData(file.getName(), taxa, null, null, null, null, importedTraits, null, 0, true);
        } catch (ImportException ie) {
            BadTraitFormatDialog dialog = new BadTraitFormatDialog(frame);
            dialog.showDialog();
        }
    }

    public void importTaxaFromTraits(final File file) throws ImportException, IOException {

        DataTable<String[]> dataTable = DataTable.Text.parse(new FileReader(file));

        String[] taxonNames = dataTable.getRowLabels();

        Taxa taxa = new Taxa();
        for (int i = 0; i < taxonNames.length; i++) {
            taxa.addTaxon(new Taxon(taxonNames[i]));
        }

        addTaxonList(taxa);

        SimpleAlignment dummyAlignment = new SimpleAlignment();
        dummyAlignment.setDataType(new DummyDataType());

        setData(file.getName(), taxa, dummyAlignment, null, null, null, null, null, 0,true);
    }

    public void importNewickFile(final File file) throws ImportException, IOException {
        FileReader reader = new FileReader(file);
        NewickImporter importer = new NewickImporter(reader);
        int treeCount = importer.countTrees();

        boolean importAllTrees = true;
        if (treeCount > 10) {
            int result = JOptionPane.showConfirmDialog(this.frame,
                    "This file contains " + treeCount + " trees.\n" +
                            "Do you want to just import the first tree in the file (recommended)?\n" +
                            "Choose 'No' to import all the trees.",
                    "Importing Data", JOptionPane.YES_NO_CANCEL_OPTION);
            if (result == JOptionPane.CANCEL_OPTION) {
                return;
            }
            importAllTrees = result == JOptionPane.NO_OPTION;
        }

        // reset file after counting trees
        importer = new NewickImporter(new FileReader(file));

        List<Tree> trees;
        if (importAllTrees) {
            trees = importer.importTrees(options.taxonList);
        } else {
            trees = Collections.singletonList(importer.importTree(options.taxonList));
        }
        setData(file.getName(), trees.get(0), null, null, null, null, null, trees, treeCount, true);
    }

    public boolean importPredictors(final File file, final TraitData trait) throws ImportException, IOException {
        List<Predictor> importedPredictors = new ArrayList<Predictor>();

        // This is the order that states will be in...
        Set<String> states = trait.getStatesOfTrait();

        String fileName = file.getName();
        String extension = "";
        if (fileName.lastIndexOf(".") != -1) {
            extension = fileName.substring(fileName.lastIndexOf(".") + 1);
            fileName = fileName.substring(0, fileName.lastIndexOf("."));
        }
        boolean isCSV = extension.toUpperCase().equals("CSV");

        DataTable<String[]> dataTable = DataTable.Text.parse(new FileReader(file), isCSV);

        String[] stateNamesCol = dataTable.getColumnLabels();
        String[] stateNamesRow = dataTable.getRowLabels();

        if (dataTable.getRowCount() != states.size()) {
            JOptionPane.showMessageDialog(frame, "Predictor data does not have the same number of rows as trait states",
                    "Mismatched states", JOptionPane.ERROR_MESSAGE);
            return false;
        }

        for (String name : stateNamesRow) {
            if (!states.contains(name)) {
                JOptionPane.showMessageDialog(frame, "Predictor row label contains unrecognized state '" + name + "'",
                        "Unrecognized state", JOptionPane.ERROR_MESSAGE);
                return false;
            }
        }

        boolean isMatrix = false;

        if (dataTable.getColumnCount() == states.size()) {
            // likely this is a matrix
            int matches = states.size();
            for (String name : stateNamesCol) {
                if (!states.contains(name)) {
                    matches -= 1;
                }
            }
            if (((double) matches) / states.size() > 0.5) { // arbitrary cut off but if 50% matching then presumably they all should
                for (String name : stateNamesCol) {
                    if (!states.contains(name)) {
                        JOptionPane.showMessageDialog(frame, "Predictor row label contains unrecognized state '" + name + "'",
                                "Unrecognized state", JOptionPane.ERROR_MESSAGE);
                        return false;
                    }
                }

                isMatrix = true;
            }
        }

        if (isMatrix) {
            double[][] data = new double[dataTable.getRowCount()][dataTable.getRowCount()];

            String name = stateNamesCol[0];
            for (int i = 0; i < dataTable.getRowCount(); i++) {
                String[] row = dataTable.getRow(i);
                for (int j = 0; j < row.length; j++) {
                    try {
                        data[i][j] = Double.parseDouble(row[j]);
                    } catch (NumberFormatException nfe) {
                        JOptionPane.showMessageDialog(frame, "Predictor '" + name + "' has a bad value at row " + (i + 1),
                                "Missing value", JOptionPane.ERROR_MESSAGE);
                    }
                }


            }
            Predictor predictor = new Predictor(options, fileName, trait, stateNamesRow, data, Predictor.Type.MATRIX);
            importedPredictors.add(predictor);

        } else {
            // a series of vectors
            for (int i = 0; i < dataTable.getColumnCount(); i++) {
                String name = stateNamesCol[i];

                double[][] data = new double[dataTable.getRowCount()][1];

                String[] values = dataTable.getColumn(i);
                for (int j = 0; j < values.length; j++) {
                    try {
                        data[j][0] = Double.parseDouble(values[j]);
                    } catch (NumberFormatException nfe) {
                        JOptionPane.showMessageDialog(frame, "Predictor '" + name + "' has a bad value at position " + (j + 1),
                                "Missing value", JOptionPane.ERROR_MESSAGE);
                    }
                }

                Predictor predictor = new Predictor(options, name, trait, stateNamesRow, data, Predictor.Type.BOTH_VECTOR);
                importedPredictors.add(predictor);
            }
        }

        for (Predictor predictor : importedPredictors) {
            trait.addPredictor(predictor);
        }

        return true;
    }

    // for Alignment
    private void setData(String fileName, TaxonList taxonList, Alignment alignment,
                         List<CharSet> charSets,
                         List<NexusApplicationImporter.TaxSet> taxSets,
                         PartitionSubstitutionModel model,
                         List<TraitData> traits, List<Tree> trees, int treeCount,
                         boolean allowEmpty) throws ImportException, IllegalArgumentException {

        String fileNameStem = Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "FA", "FAS", "FASTA", "TRE", "TREE", "XML", "TXT", "TSV"});
        if (options.fileNameStem == null || options.fileNameStem.equals(MCMCPanel.DEFAULT_FILE_NAME_STEM)) {
            options.fileNameStem = fileNameStem;
        }

        // Importing some data - see what it is and ask what to do with it...

        boolean createTreePartition = false;
        if (alignment == null && traits == null && trees != null) {
            // If only importing trees then ask if a tree-data partition should be created.
            // Doing this first in case the cancel button is pressed
            int result = JOptionPane.showConfirmDialog(this.frame, "File for importing contains trees:\nCreate 'tree-as-data' partition?",
                    "Importing Trees", JOptionPane.YES_NO_CANCEL_OPTION);
            if (result == JOptionPane.CANCEL_OPTION) {
                return;
            }
            createTreePartition = result == JOptionPane.YES_OPTION;
        }

        boolean createTraitPartition = false;
        if (alignment == null && trees == null && traits != null && !traits.isEmpty()) {
            // If only importing traits then ask if a trait-data partition should be created.
            // Doing this first in case the cancel button is pressed
            int result = JOptionPane.showConfirmDialog(this.frame, "File for importing contains traits:\nCreate a traits partition?",
                    "Importing Traits", JOptionPane.YES_NO_CANCEL_OPTION);
            if (result == JOptionPane.CANCEL_OPTION) {
                return;
            }
            createTraitPartition = result == JOptionPane.YES_OPTION;
        }

        boolean compressPatterns = true;
        if (alignment != null) {
            // check the alignment before adding it...
            for (Sequence seq : alignment.getSequences()) {
                if (seq.getLength() != alignment.getSiteCount()) {
                    // sequences are different lengths
                    throw new ImportException("The sequences in the alignment file are of different lengths - BEAST requires aligned sequences");
                }
            }

            if (alignment.getSiteCount() == 0) {
                if (allowEmpty) {
                    int result = JOptionPane.showConfirmDialog(this.frame, "File contains zero length sequences:\nCreate empty partition?",
                            "Importing zero length sequences", JOptionPane.OK_CANCEL_OPTION);
                    if (result == JOptionPane.CANCEL_OPTION) {
                        return;
                    }
                } else {
                    throw new ImportException("This alignment is of zero length");
                }
            }

            if (charSets == null && alignment.getSiteCount() > 50000) {
//                int result = JOptionPane.showConfirmDialog(this.frame,
//                        "The alignment contains long sequences (" + alignment.getSiteCount() + " sites):\n" +
//                                "Do you want to compress the alignment into site-patterns? This will\n" +
//                                "take some time and will preclude some analyses but will make the\n" +
//                                "BEAST XML file size smaller and faster to load.",
//                        "Importing long sequences", JOptionPane.YES_NO_CANCEL_OPTION);
//                if (result == JOptionPane.CANCEL_OPTION) {
//                    return;
//                }
//                compressPatterns = result == JOptionPane.YES_OPTION;
                // for now, if the sequences are long, we just don't do a pattern compression (this just
                // means the table doesn't show the pattern count. Future options may be to compress
                // the data to only variable sites and a count of constant sites.
                compressPatterns = false;
            }

        }

        List<String> attributeNames = new ArrayList<>();

        // first add them in the order they come in the imported trait file (the taxon attributes
        // may be out of order by using a Set).
        if (traits != null) {
            for (TraitData trait : traits) {
                attributeNames.add(trait.getName());
            }
        }

        for (Taxon taxon : taxonList) {
            if (taxon.getDate() != null) {
                options.useTipDates = true;
            }
            for (Iterator<String> it = taxon.getAttributeNames(); it.hasNext(); ) {
                String name = it.next();
                if (!name.equalsIgnoreCase("Date") && !name.equals("_height") &&
                        !attributeNames.contains(name)) {
                    // "_height" is a magic value used by BEAUti
                    attributeNames.add(name);
                }
            }
        }

        // attempt to work out the type of the trait
        for (String name : attributeNames) {
            if (!options.traitExists(name)) {
                TraitData.TraitType type = null;
                for (Taxon taxon : taxonList) {
                    if (taxon.getAttribute(name) != null) {
                        String value = taxon.getAttribute(name).toString();
                        if (value.equals("NA") || value.equals("?")) { // need to check "?" to avoid 'else' block
                            taxon.setAttribute(name, "?");
                        } else {
                            try {
                                Integer.parseInt(value);
                                if (type == null || type == TraitData.TraitType.INTEGER) {
                                    type = TraitData.TraitType.INTEGER;
                                }
                            } catch (NumberFormatException e) {
                                try {
                                    Double.parseDouble(value);
                                    if (type == null || type == TraitData.TraitType.INTEGER || type == TraitData.TraitType.CONTINUOUS) {
                                        type = TraitData.TraitType.CONTINUOUS;
                                    }
                                } catch (NumberFormatException e1) {
                                    type = TraitData.TraitType.DISCRETE;
                                }
                            }
                        }
                    }
                }
                if (type == null) {
                    type = TraitData.TraitType.DISCRETE;
                }

                options.traits.add(new TraitData(options, name, fileName, type));
            }
        }

        addTaxonList(taxonList);

        if (alignment != null) {
            addAlignment(alignment, charSets, compressPatterns, model, fileName, fileNameStem);
        }

        if (taxSets != null) {
            addTaxonSets(taxonList, taxSets);
        }

        if (traits != null) {
            addTraits(traits);
        }

        TreeHolder treeHolder = addTrees(fileNameStem, fileName, trees, treeCount);

        if (createTreePartition) {
            options.createPartitionForTree(treeHolder, fileNameStem);
        }

        if (createTraitPartition) {
            frame.getDataPanel().createPartitionFromTraits(null, fileNameStem, frame);
        }
    }

    private void addTaxonList(TaxonList taxonList) throws ImportException {
        checkTaxonList(taxonList);
        if (options.taxonList == null) {
            // This is the first partition to be loaded...
            options.taxonList = new Taxa(taxonList);
        } else {
            // otherwise just add the new ones...
            for (Taxon taxon : taxonList) {
                if (!options.taxonList.contains(taxon)) {
                    options.taxonList.addTaxon(taxon);
                }
            }
        }
    }

    private void checkTaxonList(TaxonList taxonList) throws ImportException {

        // check the taxon names for invalid characters
        boolean foundAmp = false;
        for (Taxon taxon : taxonList) {
            String name = taxon.getId();
            if (name.indexOf('&') >= 0) {
                foundAmp = true;
            }
        }
        if (foundAmp) {
            throw new ImportException("One or more taxon names include an illegal character ('&').\n" +
                    "These characters will prevent BEAST from reading the resulting XML file.\n\n" +
                    "Please edit the taxon name(s) before reloading the data file.");
        }

        // make sure they all have dates...
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            if (taxonList.getTaxonAttribute(i, "date") == null) {
                Date origin = new Date(0);

                dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                taxonList.getTaxon(i).setAttribute("date", date);
            }
        }
    }

    private void addTaxonSets(TaxonList taxonList, List<NexusApplicationImporter.TaxSet> taxSets) throws ImportException {
        assert taxSets != null;
        for (NexusApplicationImporter.TaxSet taxSet : taxSets) {
            Taxa taxa = new Taxa(taxSet.getName());
            for (CharSetBlock block : taxSet.getBlocks()) {
                for (int i = block.getFromSite(); i <= block.getToSite(); i++) {
                    taxa.addTaxon(taxonList.getTaxon(i - 1));
                }
            }
            options.taxonSets.add(taxa);
            options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
            options.taxonSetsMono.put(taxa, false);
            options.taxonSetsIncludeStem.put(taxa, false);
            options.taxonSetsHeights.put(taxa, 0.0);
        }
    }

    private void addAlignment(Alignment alignment,
                              List<CharSet> charSets,
                              boolean compressPatterns,
                              PartitionSubstitutionModel model,
                              String fileName, String fileNameStem) {
        assert alignment != null;

        List<AbstractPartitionData> partitions = new ArrayList<AbstractPartitionData>();
        if (charSets != null && !charSets.isEmpty()) {
            for (CharSet charSet : charSets) {
                PartitionData partitionData = new PartitionData(options, charSet.name, fileName,
                        charSet.constructCharSetAlignment(alignment));
                if (compressPatterns) {
                    partitionData.compressPatterns(SitePatterns.CompressionType.UNIQUE_ONLY);
                }
                partitions.add(partitionData);
            }
        } else {
            PartitionData partitionData = new PartitionData(options, fileNameStem, fileName, alignment);
            if (compressPatterns) {
                partitionData.compressPatterns(SitePatterns.CompressionType.UNIQUE_ONLY);
            }
            partitions.add(partitionData);
        }
        createPartitionFramework(model, partitions);
    }


    private void createPartitionFramework(PartitionSubstitutionModel model, List<AbstractPartitionData> partitions) {
        for (AbstractPartitionData partition : partitions) {
            String name = partition.getName();

            while (name.isEmpty() || options.hasPartitionData(name)) {
                String text;
                if (options.hasPartitionData(name)) {
                    text = "<html>" +
                            "A partition named, " + name + ", already exists.<br>" +
                            "Please provide a unique name for this partition." +
                            "</html>";
                } else {
                    text = "<html>" +
                            "Invalid partition name. Please provide a unique<br>" +
                            "name for this partition." +
                            "</html>";
                }
                if (partitionNameDialog == null) {
                    partitionNameDialog = new PartitionNameDialog(frame);
                }

                partitionNameDialog.setDescription(text);

                int result = partitionNameDialog.showDialog();

                if (result == -1 || result == JOptionPane.CANCEL_OPTION) {
                    return;
                }

                name = partitionNameDialog.getName();
            }
            partition.setName(name);

            options.dataPartitions.add(partition);

            if (model != null) {
                setSubstModel(partition, model);

                options.setClockAndTree(partition);//TODO Cannot load Clock Model and Tree Model from BEAST file yet

            } else {// only this works
                if (options.getPartitionSubstitutionModels(partition.getDataType()).isEmpty()) {// use same substitution model in beginning
                    // PartitionSubstitutionModel based on PartitionData
                    PartitionSubstitutionModel psm = new PartitionSubstitutionModel(options, DEFAULT_NAME, partition);
                    partition.setPartitionSubstitutionModel(psm);
                } else { //if (options.getPartitionSubstitutionModels() != null) {
//                        && options.getPartitionSubstitutionModels().size() == 1) {
                    PartitionSubstitutionModel psm = options.getPartitionSubstitutionModels(partition.getDataType()).get(0);
                    setSubstModel(partition, psm);
                }

                options.setClockAndTree(partition);
            }
        }
    }

    private void setSubstModel(AbstractPartitionData partition, PartitionSubstitutionModel psm) {
        partition.setPartitionSubstitutionModel(psm);

        if (psm.getDataType().getType() != partition.getDataType().getType())
            throw new IllegalArgumentException("Partition " + partition.getName()
                    + "\ncannot assign to Substitution Model\n" + psm.getName()
                    + "\nwith different data type.");
    }

    private void addTraits(List<TraitData> traits) {
        assert traits != null;
        for (TraitData trait : traits) {
            options.addTrait(trait);
        }
    }

    /**
     * Add a new tree/tree set
     * @param fileNameStem
     * @param fileName
     * @param trees
     * @param treeCount the number of trees - if > than the size of trees then is a set of trees which haven't
     *                  all been loaded.
     * @return the TreeHolder
     */
    private TreeHolder addTrees(String fileNameStem, String fileName, List<Tree> trees, int treeCount) {
        TreeHolder treeHolder = null;
        if (trees != null && !trees.isEmpty()) {
            int i = 1;
            for (Tree tree : trees) {
                String id = tree.getId();
                if (id == null || id.trim().isEmpty()) {
                    tree.setId("tree_" + i);
                }
                i++;
            }
            treeHolder = new TreeHolder(trees, treeCount, fileNameStem, fileName);
            options.userTrees.put(fileNameStem, treeHolder);
        }

        return treeHolder;
    }

}
