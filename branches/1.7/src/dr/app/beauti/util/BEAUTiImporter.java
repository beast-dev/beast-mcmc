/*
 * BEAUTiImporter.java
 *
 * Copyright (C) 2002-2011 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.util;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.mcmcpanel.MCMCPanel;
import dr.app.beauti.options.*;
import dr.app.util.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.MicroSatImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.io.NexusImporter.NexusBlock;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.util.DataTable;
import org.jdom.JDOMException;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;


/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class BEAUTiImporter {

    private final BeautiOptions options;
    private final BeautiFrame frame;

    private PartitionNameDialog partitionNameDialog = null;

    public BEAUTiImporter(BeautiFrame frame, BeautiOptions options) {
        this.frame = frame;
        this.options = options;
    }

    public void importFromFile(File file) throws IOException, ImportException, JDOMException {
        try {
            Reader reader = new FileReader(file);

            BufferedReader bufferedReader = new BufferedReader(reader);
            String line = bufferedReader.readLine();
            while (line != null && line.length() == 0) {
                line = bufferedReader.readLine();
            }

            if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
                // is a NEXUS file
                importNexusFile(file);
            } else if ((line != null && line.trim().startsWith("" + FastaImporter.FASTA_FIRST_CHAR))) {
                // is a FASTA file
                importFastaFile(file);
            } else if ((line != null && (line.toUpperCase().contains("<?XML") || line.toUpperCase().contains("<BEAST")))) {
                // assume it is a BEAST XML file and see if that works...
                importBEASTFile(file);
//            } else {
//                // assume it is a tab-delimited traits file and see if that works...
//                importTraits(file);
            } else if ((line != null && line.toUpperCase().contains("#MICROSAT"))) {
                // MicroSatellite
                importMicroSatFile(file);
            } else {
                throw new ImportException("Unrecognized format for imported file.");
            }

            bufferedReader.close();
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    // micro-sat
    private void importMicroSatFile(File file) throws IOException, ImportException {
        try {
            Reader reader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(reader);

            MicroSatImporter importer = new MicroSatImporter(bufferedReader);

            List<Patterns> microsatPatList = importer.importPatterns();
            Taxa unionSetTaxonList = importer.getUnionSetTaxonList();
            Microsatellite microsatellite = importer.getMicrosatellite();
//            options.allowDifferentTaxa = importer.isHasDifferentTaxon();

            bufferedReader.close();

            PartitionSubstitutionModel substModel = new PartitionSubstitutionModel(options, microsatPatList.get(0).getId());
            substModel.setMicrosatellite(microsatellite);

            for (Patterns patterns : microsatPatList) {
                setData(file.getName(), unionSetTaxonList, patterns, substModel, null);
            }
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
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
                setData(name, taxa, alignment, null, null, null, null);

                count++;
            }

            // assume that any additional taxon lists are taxon sets...
            for (int i = 1; i < taxonLists.size(); i++) {
                Taxa taxonSet = (Taxa) taxonLists.get(i);

                options.taxonSets.add(taxonSet);
                options.taxonSetsMono.put(taxonSet, false);
                options.taxonSetsIncludeStem.put(taxonSet, false);
                options.taxonSetsTreeModel.put(taxonSet, options.getPartitionTreeModels().get(0));
            }

            reader.close();
        } catch (JDOMException e) {
            throw new JDOMException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }

    }

    // nexus
    private void importNexusFile(File file) throws IOException, ImportException {
        TaxonList taxa = null;
        SimpleAlignment alignment = null;
        List<Tree> trees = new ArrayList<Tree>();
        PartitionSubstitutionModel model = null;
        List<NexusApplicationImporter.CharSet> charSets = new ArrayList<NexusApplicationImporter.CharSet>();

        try {
            FileReader reader = new FileReader(file);

            NexusApplicationImporter importer = new NexusApplicationImporter(reader);

            boolean done = false;

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

                        // I guess there is no reason not to allow multiple trees blocks
//                        if (trees.size() > 0) {
//                            throw new MissingBlockException("TREES block already defined");
//                        }

                        Tree[] treeArray = importer.parseTreesBlock(taxa);
                        trees.addAll(Arrays.asList(treeArray));

                        if (taxa == null && trees.size() > 0) {
                            taxa = trees.get(0);
                        }


                    } else if (block == NexusApplicationImporter.PAUP_BLOCK) {

                        model = importer.parsePAUPBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

                        model = importer.parseMrBayesBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.ASSUMPTIONS_BLOCK) {

                        importer.parseAssumptionsBlock(charSets);

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            reader.close();

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

        setData(file.getName(), taxa, alignment, charSets, model, null, trees);
    }

    // FASTA

    private void importFastaFile(File file) throws IOException, ImportException {
        try {
            FileReader reader = new FileReader(file);

            FastaImporter importer = new FastaImporter(reader, Nucleotides.INSTANCE);

            Alignment alignment = importer.importAlignment();

            reader.close();

            setData(file.getName(), alignment, alignment, null, null, null, null);
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    private boolean isMissingValue(String value) {
        return (value.equals("?") || value.equals("NA") || value.length() == 0);
    }

    public void importTraits(final File file) throws Exception {
        List<TraitData> importedTraits = new ArrayList<TraitData>();
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
                                !warningGiven ) {
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

            if (validateTraitName(traitName)) {
                importedTraits.add(newTrait);
            }

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
                    taxon.setAttribute(traitName, "?");
                }
                j++;
            }
        }
        setData(file.getName(), taxa, null, null, null, importedTraits, null);
    }

//    public void oldImportTraits(final File file) throws Exception {
//        List<TraitData> importedTraits = new ArrayList<TraitData>();
//        Taxa taxa = new Taxa();
//
//        try {
//            Map<String, List<String[]>> traits = Utils.importTraitsFromFile(file, "\t");
//
//
//            for (Map.Entry<String, List<String[]>> e : traits.entrySet()) {
//                final Class c = Utils.detectType(e.getValue().get(0)[1]);
//                final String traitName = e.getKey();
//
//                Boolean warningGiven = false;
//                for (String[] v : e.getValue()) {
//                    final Class c1 = Utils.detectType(v[1]);
//                    if (c != c1 && !warningGiven) {
//                        JOptionPane.showMessageDialog(frame, "Not all values of same type in column" + traitName,
//                                "Incompatible values", JOptionPane.WARNING_MESSAGE);
//                        warningGiven = true;
//                        // TODO Error - not all values of same type
//                    }
//                }
//
//                TraitData.TraitType t = (c == Boolean.class || c == String.class) ? TraitData.TraitType.DISCRETE :
//                        (c == Integer.class) ? TraitData.TraitType.INTEGER : TraitData.TraitType.CONTINUOUS;
//                TraitData newTrait = new TraitData(options, traitName, file.getName(), t);
//
////                TraitData newTrait = new TraitData(options, traitName, file.getName(), TraitData.TraitType.DISCRETE);
//
//                if (validateTraitName(traitName)) {
//                    importedTraits.add(newTrait);
//                }
//
//                for (final String[] v : e.getValue()) {
//                    if (v[0].equalsIgnoreCase(v[1])) {
//                        throw new Arguments.ArgumentException("Trait (" + traitName + ") value (" + v[1]
//                                + ")\n cannot be same as taxon name (" + v[0] + ") !");
//                    }
//
//                    final int index = options.taxonList.getTaxonIndex(v[0]);
//                    Taxon taxon;
//                    if (index >= 0) {
//                        taxon = options.taxonList.getTaxon(index);
////                        taxon.setAttribute(traitName, Utils.constructFromString(c, v[1]));
//                    } else {
//                        taxon = new Taxon(v[0]);
//                    }
//                    taxon.setAttribute(traitName, v[1]);
//                    taxa.addTaxon(taxon);
//                }
//            }
//        } catch (Arguments.ArgumentException e) {
//            JOptionPane.showMessageDialog(frame, "Error in loading traits file " + file.getName() + " :\n" + e.getMessage(),
//                    "Error Loading file", JOptionPane.ERROR_MESSAGE);
//            // AR: this will remove all the existing traits including those loaded previously:
////            traitsPanel.traitsTable.selectAll();
////            traitsPanel.removeTrait();
//        }
//
//        setData(file.getName(), taxa, null, null, null, importedTraits, null);
//    }

    public boolean validateTraitName(String traitName) {
        // check that the name is valid
        if (traitName.trim().length() == 0) {
            Toolkit.getDefaultToolkit().beep();
            return false;
        }

        // disallow a trait called 'date'
        if (traitName.equalsIgnoreCase("date")) {
            JOptionPane.showMessageDialog(frame,
                    "This trait name has a special meaning. Use the 'Tip Date' panel\n" +
                            " to set dates for taxa.",
                    "Reserved trait name",
                    JOptionPane.WARNING_MESSAGE);

            return false;
        }

        // check that the trait name doesn't exist
        if (options.traitExists(traitName)) {
            int option = JOptionPane.showConfirmDialog(frame,
                    "A trait of this name already exists. Do you wish to replace\n" +
                            "it with this new trait? This may result in the loss or change\n" +
                            "in trait values for the taxa.",
                    "Overwrite trait?",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE);

            if (option == JOptionPane.NO_OPTION) {
                return false;
            }
        }

        return true;
    }

    // for Alignment
    private void setData(String fileName, TaxonList taxonList, Alignment alignment,
                         List<NexusApplicationImporter.CharSet> charSets, PartitionSubstitutionModel model,
                         List<TraitData> traits, List<Tree> trees) throws ImportException, IllegalArgumentException {
        String fileNameStem = Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "TRE", "TREE", "XML", "TXT"});
        if (options.fileNameStem == null || options.fileNameStem.equals(MCMCPanel.DEFAULT_FILE_NAME_STEM)) {
            options.fileNameStem = fileNameStem;
        }

        addTaxonList(taxonList);

        addAlignment(alignment, charSets, model, fileName, fileNameStem);

        addTraits(traits);

        addTrees(trees);
    }

    // for Patterns
    private void setData(String fileName, Taxa taxonList, Patterns patterns,
                         PartitionSubstitutionModel model, List<TraitData> traits
    ) throws ImportException, IllegalArgumentException {
        String fileNameStem = Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "TRE", "TREE", "XML", "TXT"});
        if (options.fileNameStem == null || options.fileNameStem.equals(MCMCPanel.DEFAULT_FILE_NAME_STEM)) {
            options.fileNameStem = fileNameStem;
        }

        addTaxonList(taxonList);

        addPatterns(patterns, model, fileName);

        addTraits(traits);
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

    private void addAlignment(Alignment alignment, List<NexusApplicationImporter.CharSet> charSets,
                              PartitionSubstitutionModel model,
                              String fileName, String fileNameStem) {
        if (alignment != null) {
            List<AbstractPartitionData> partitions = new ArrayList<AbstractPartitionData>();
            if (charSets != null && charSets.size() > 0) {
                for (NexusApplicationImporter.CharSet charSet : charSets) {
                    partitions.add(new PartitionData(options, charSet.name, fileName,
                            charSet.constructCharSetAlignment(alignment)));
                }
            } else {
                partitions.add(new PartitionData(options, fileNameStem, fileName, alignment));
            }
            createPartitionFramework(model, partitions);
        }
    }

    private void addPatterns(Patterns patterns, PartitionSubstitutionModel model, String fileName) {
        if (patterns != null) {
            List<AbstractPartitionData> partitions = new ArrayList<AbstractPartitionData>();
            partitions.add(new PartitionPattern(options, patterns.getId(), fileName, patterns));

            createPartitionFramework(model, partitions);
        }
    }

    private void createPartitionFramework(PartitionSubstitutionModel model, List<AbstractPartitionData> partitions) {
        for (AbstractPartitionData partition : partitions) {
            String name = partition.getName();

            while (name.length() == 0 || options.hasPartitionData(name)) {
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

                setClockAndTree(partition);//TODO Cannot load Clock Model and Tree Model from BEAST file yet

            } else {// only this works
                if (options.getPartitionSubstitutionModels(partition.getDataType()).size() < 1) {// use same substitution model in beginning
                    // PartitionSubstitutionModel based on PartitionData
                    PartitionSubstitutionModel psm = new PartitionSubstitutionModel(options, partition);
                    partition.setPartitionSubstitutionModel(psm);
                } else { //if (options.getPartitionSubstitutionModels() != null) {
//                        && options.getPartitionSubstitutionModels().size() == 1) {
                    PartitionSubstitutionModel psm = options.getPartitionSubstitutionModels(partition.getDataType()).get(0);
                    setSubstModel(partition, psm);
                }

                setClockAndTree(partition);
            }
        }

        options.updatePartitionAllLinks();
        options.clockModelOptions.initClockModelGroup();
    }

    private void setClockAndTree(AbstractPartitionData partition) {
        // use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
        if (options.getPartitionClockModels(partition.getDataType()).size() < 1) {
            // PartitionClockModel based on PartitionData
            PartitionClockModel pcm = new PartitionClockModel(options, partition);
            partition.setPartitionClockModel(pcm);
        } else { //if (options.getPartitionClockModels() != null) {
//                        && options.getPartitionClockModels().size() == 1) {
            PartitionClockModel pcm;
            if (partition.getDataType().getType() == DataType.MICRO_SAT) {
                pcm = new PartitionClockModel(options, partition);
            } else {
                // make sure in the same data type
                pcm = options.getPartitionClockModels(partition.getDataType()).get(0);
            }
            partition.setPartitionClockModel(pcm);
        }

        // use same tree model and same tree prior in beginning
        if (options.getPartitionTreeModels().size() < 1) {
            // PartitionTreeModel based on PartitionData
            PartitionTreeModel ptm = new PartitionTreeModel(options, partition);
            partition.setPartitionTreeModel(ptm);

            // PartitionTreePrior always based on PartitionTreeModel
            PartitionTreePrior ptp = new PartitionTreePrior(options, ptm);
            ptm.setPartitionTreePrior(ptp);
        } else { //if (options.getPartitionTreeModels() != null) {
//                        && options.getPartitionTreeModels().size() == 1) {
            PartitionTreeModel ptm;
            if (partition.getDataType().getType() == DataType.MICRO_SAT) {
                ptm = new PartitionTreeModel(options, partition); // different tree model,
                PartitionTreePrior ptp = options.getPartitionTreePriors().get(0); // but same tree prior
                ptm.setPartitionTreePrior(ptp);
            } else {
                ptm = options.getPartitionTreeModels().get(0); // same tree model,
            }
            partition.setPartitionTreeModel(ptm); // if same tree model, therefore same prior
        }
    }

    private void setSubstModel(AbstractPartitionData partition, PartitionSubstitutionModel psm) {
        partition.setPartitionSubstitutionModel(psm);

        if (psm.getDataType() != partition.getDataType())
            throw new IllegalArgumentException("Partition " + partition.getName()
                    + "\ncannot assign to Substitution Model\n" + psm.getName()
                    + "\nwith different data type.");
    }

    private void addTraits(List<TraitData> traits) {
        if (traits != null) {
            for (TraitData trait : traits) {
                options.addTrait(trait);
            }

            options.updatePartitionAllLinks();
        }
    }

    private void addTrees(List<Tree> trees) {
        if (trees != null && trees.size() > 0) {
            for (Tree tree : trees) {
                String id = tree.getId();
                if (id == null || id.trim().length() == 0) {
                    tree.setId("tree_" + (options.userTrees.size() + 1));
                } else {
                    String newId = id;
                    int count = 1;
                    for (Tree tree1 : options.userTrees) {
                        if (tree1.getId().equals(newId)) {
                            newId = id + "_" + count;
                            count++;
                        }
                    }
                    tree.setId(newId);
                }
                options.userTrees.add(tree);
            }
        }
    }

}
