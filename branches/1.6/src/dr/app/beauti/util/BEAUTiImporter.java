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
import dr.app.beauti.options.*;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.io.NexusImporter.NexusBlock;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import org.jdom.JDOMException;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;


/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class BEAUTiImporter {

    private final BeautiOptions options;
    private final BeautiFrame frame;

    public BEAUTiImporter(BeautiFrame frame, BeautiOptions options) {
        this.frame = frame;
        this.options = options;
    }

    public void importFromFile(File file) throws Exception {
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
            } else {
                throw new ImportException("Unrecognized format for imported file.");
            }
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    // xml

    private void importBEASTFile(File file) throws IOException, JDOMException, ImportException {
        try {
            FileReader reader = new FileReader(file);

            BeastImporter importer = new BeastImporter(reader);

            List<TaxonList> taxonLists = new ArrayList<TaxonList>();
            List<Alignment> alignments = new ArrayList<Alignment>();

            importer.importBEAST(taxonLists, alignments);

            TaxonList taxa = taxonLists.get(0);

            for (Alignment alignment : alignments) {
                setData(taxa, alignment, null, null, null, null, file.getName());
            }
        } catch (JDOMException e) {
            throw new JDOMException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }

    }

    // nexus

    private void importNexusFile(File file) throws Exception {
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

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (alignment == null && taxa == null) {
                throw new MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (IOException e) {
            throw new IOException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (Exception e) {
            throw new Exception(e.getMessage());
        }

        setData(taxa, alignment, trees, null, model, charSets, file.getName());
    }

    // FASTA

    private void importFastaFile(File file) throws Exception {
        try {
            FileReader reader = new FileReader(file);

            FastaImporter importer = new FastaImporter(reader, Nucleotides.INSTANCE);

            Alignment alignment = importer.importAlignment();

            setData(alignment, alignment, null, null, null, null, file.getName());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    public void importTraits(final File file) throws Exception {
        List<TraitData> importedTraits = new ArrayList<TraitData>();
        Taxa taxa = new Taxa();

        try {
            Map<String, List<String[]>> traits = Utils.importTraitsFromFile(file, "\t");


            for (Map.Entry<String, List<String[]>> e : traits.entrySet()) {
                final Class c = Utils.detectTYpe(e.getValue().get(0)[1]);
                final String traitName = e.getKey();

                Boolean warningGiven = false;
                for (String[] v : e.getValue()) {
                    final Class c1 = Utils.detectTYpe(v[1]);
                    if (c != c1 && !warningGiven) {
                        JOptionPane.showMessageDialog(frame, "Not all values of same type in column" + traitName,
                                "Incompatible values", JOptionPane.WARNING_MESSAGE);
                        warningGiven = true;
                        // TODO Error - not all values of same type
                    }
                }

                TraitData.TraitType t = (c == Boolean.class || c == String.class) ? TraitData.TraitType.DISCRETE :
                        (c == Integer.class) ? TraitData.TraitType.INTEGER : TraitData.TraitType.CONTINUOUS;
                TraitData newTrait = new TraitData(options, traitName, file.getName(), t);

//                TraitData newTrait = new TraitData(options, traitName, file.getName(), TraitData.TraitType.DISCRETE);

                if (validateTraitName(traitName)) {
                    importedTraits.add(newTrait);
                }

                for (final String[] v : e.getValue()) {
                    if (v[0].equalsIgnoreCase(v[1])) {
                        throw new Arguments.ArgumentException("Trait (" + traitName + ") value (" + v[1]
                                + ")\n cannot be same as taxon name (" + v[0] + ") !");
                    }

                    final int index = options.taxonList.getTaxonIndex(v[0]);
                    Taxon taxon;
                    if (index >= 0) {
                        taxon = options.taxonList.getTaxon(index);
//                        taxon.setAttribute(traitName, Utils.constructFromString(c, v[1]));
                    } else {
                        taxon = new Taxon(v[0]);
                    }
                    taxon.setAttribute(traitName, v[1]);
                    taxa.addTaxon(taxon);
                }
            }
        } catch (Arguments.ArgumentException e) {
            JOptionPane.showMessageDialog(frame, "Error in loading traits file " + file.getName() + " :\n" + e.getMessage(),
                    "Error Loading file", JOptionPane.ERROR_MESSAGE);
            // AR: this will remove all the existing traits including those loaded previously:
//            traitsPanel.traitsTable.selectAll();
//            traitsPanel.removeTrait();
        }

        setData(taxa, null, null, importedTraits, null, null, file.getName());
    }

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
        if (options.containTrait(traitName)) {
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


    //TODO need refactory to simplify

    private void setData(TaxonList taxa, Alignment alignment, List<Tree> trees, List<TraitData> traits,
                         PartitionSubstitutionModel model,
                         List<NexusApplicationImporter.CharSet> charSets, String fileName) throws ImportException {
        String fileNameStem = dr.app.util.Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "TRE", "TREE", "XML"});

        if (options.taxonList == null) {
            // This is the first partition to be loaded...

            options.taxonList = new Taxa(taxa);

            // check the taxon names for invalid characters
            boolean foundAmp = false;
            for (Taxon taxon : taxa) {
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
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                if (taxa.getTaxonAttribute(i, "date") == null) {
                    Date origin = new Date(0);

                    dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                    taxa.getTaxon(i).setAttribute("date", date);
                }
            }

            options.fileNameStem = fileNameStem;

        } else {
            // AR - removed this distinction. I think we should always allow different taxa
            // for different partitions but give a warning if they are different

            // This is an additional partition so check it uses the same taxa
            if (!options.allowDifferentTaxa) { // not allow Different Taxa
                List<String> oldTaxa = new ArrayList<String>();
                for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                    oldTaxa.add(options.taxonList.getTaxon(i).getId());
                }
                List<String> newTaxa = new ArrayList<String>();
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    newTaxa.add(taxa.getTaxon(i).getId());
                }

                if (!(oldTaxa.containsAll(newTaxa) && oldTaxa.size() == newTaxa.size())) {
                    // set to Allow Different Taxa
                    options.allowDifferentTaxa = true;
                    //changeTabs();// can be added, if required in future

                    List<String> prevTaxa = new ArrayList<String>();
                    for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                        prevTaxa.add(options.taxonList.getTaxon(i).getId());
                    }
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        if (!prevTaxa.contains(taxa.getTaxon(i).getId())) {
                            options.taxonList.addTaxon(taxa.getTaxon(i));
                        }
                    }
                }
            } else { // allow Different Taxa
                // AR - it will be much simpler just to consider options.taxonList
                // to be the union set of all taxa. Each data partition has an alignment
                // which is a taxon list containing the taxa specific to that partition

                // add the new diff taxa
                List<String> prevTaxa = new ArrayList<String>();
                for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                    prevTaxa.add(options.taxonList.getTaxon(i).getId());
                }
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    if (!prevTaxa.contains(taxa.getTaxon(i).getId())) {
                        options.taxonList.addTaxon(taxa.getTaxon(i));
                    }
                }

            }
        }

        addAlignment(alignment, charSets, model, fileName, fileNameStem);

        addTraits(traits, fileName, fileNameStem);

        addTrees(trees);
    }

    private void addAlignment(Alignment alignment, List<NexusApplicationImporter.CharSet> charSets,
                              PartitionSubstitutionModel model,
                              String fileName, String fileNameStem) {
        if (alignment != null) {
            List<PartitionData> partitions = new ArrayList<PartitionData>();
            if (charSets != null && charSets.size() > 0) {
                for (NexusApplicationImporter.CharSet charSet : charSets) {
                    partitions.add(new PartitionData(options, charSet.getName(), fileName,
                            charSet.constructCharSetAlignment(alignment)));
                }
            } else {
                partitions.add(new PartitionData(options, fileNameStem, fileName, alignment));
            }
            for (PartitionData partition : partitions) {
                options.dataPartitions.add(partition);

                if (model != null) {//TODO Cannot load Clock Model and Tree Model from BEAST file yet
                    partition.setPartitionSubstitutionModel(model);
//                    model.addPartitionData(partition);

                    // use same tree model and same tree prior in beginning
                    if (options.getPartitionTreeModels() != null
                            && options.getPartitionTreeModels().size() == 1) {
                        PartitionTreeModel ptm = options.getPartitionTreeModels().get(0);
                        partition.setPartitionTreeModel(ptm); // same tree model, therefore same prior
                    }
                    if (partition.getPartitionTreeModel() == null) {
                        // PartitionTreeModel based on PartitionData
                        PartitionTreeModel ptm = new PartitionTreeModel(options, partition);
                        partition.setPartitionTreeModel(ptm);

                        // PartitionTreePrior always based on PartitionTreeModel
                        PartitionTreePrior ptp = new PartitionTreePrior(options, ptm);
                        ptm.setPartitionTreePrior(ptp);

//                        options.addPartitionTreeModel(ptm);
//                        options.shareSameTreePrior = true;
                    }

                    // use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
                    if (options.getPartitionClockModels() != null
                            && options.getPartitionClockModels().size() == 1) {
                        PartitionClockModel pcm = options.getPartitionClockModels().get(0);
                        partition.setPartitionClockModel(pcm);
                    }
                    if (partition.getPartitionClockModel() == null) {
                        // PartitionClockModel based on PartitionData
                        PartitionClockModel pcm = new PartitionClockModel(options, partition);
                        partition.setPartitionClockModel(pcm);
//                        options.addPartitionClockModel(pcm);
                    }
//                    options.clockModelOptions.fixRateOfFirstClockPartition();

                } else {// only this works
                    if (options.getPartitionSubstitutionModels() != null
                            && options.getPartitionSubstitutionModels().size() == 1) { // use same substitution model in beginning
                        PartitionSubstitutionModel psm = options.getPartitionSubstitutionModels().get(0);
                        if (psm.getDataType() == alignment.getDataType()) {
                            partition.setPartitionSubstitutionModel(psm);
                        } else {
                            //TODO exception?
                        }
                    }
                    if (partition.getPartitionSubstitutionModel() == null) {
                        // PartitionSubstitutionModel based on PartitionData
                        PartitionSubstitutionModel psm = new PartitionSubstitutionModel(options, partition);
                        partition.setPartitionSubstitutionModel(psm);
//                        options.addPartitionSubstitutionModel(psm);
                    }

                    // use same tree model and same tree prior in beginning
                    if (options.getPartitionTreeModels() != null
                            && options.getPartitionTreeModels().size() == 1) {
                        PartitionTreeModel ptm = options.getPartitionTreeModels().get(0);
                        partition.setPartitionTreeModel(ptm); // same tree model, therefore same prior
                    }
                    if (partition.getPartitionTreeModel() == null) {
                        // PartitionTreeModel based on PartitionData
                        PartitionTreeModel ptm = new PartitionTreeModel(options, partition);
                        partition.setPartitionTreeModel(ptm);

                        // PartitionTreePrior always based on PartitionTreeModel
                        PartitionTreePrior ptp = new PartitionTreePrior(options, ptm);
                        ptm.setPartitionTreePrior(ptp);

//                        options.addPartitionTreeModel(ptm);
//                        options.shareSameTreePrior = true;
                    }

                    // use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
                    if (options.getPartitionClockModels() != null
                            && options.getPartitionClockModels().size() == 1) {
                        PartitionClockModel pcm = options.getPartitionClockModels().get(0);
                        partition.setPartitionClockModel(pcm);
                    }
                    if (partition.getPartitionClockModel() == null) {
                        // PartitionClockModel based on PartitionData
                        PartitionClockModel pcm = new PartitionClockModel(options, partition);
                        partition.setPartitionClockModel(pcm);
//                        options.addPartitionClockModel(pcm);
                    }
                }
            }

//            options.updateLinksBetweenPDPCMPSMPTMPTPP();
            options.updatePartitionAllLinks();
            options.clockModelOptions.fixRateOfFirstClockPartition();
        }
    }

    private void addTraits(List<TraitData> traits,
                           String fileName, String fileNameStem) {
        if (traits != null) {
            for (TraitData trait : traits) {
                options.addTrait(trait);
            }

//            options.updateLinksBetweenPDPCMPSMPTMPTPP();
            options.updatePartitionAllLinks();
            options.clockModelOptions.fixRateOfFirstClockPartition();
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
