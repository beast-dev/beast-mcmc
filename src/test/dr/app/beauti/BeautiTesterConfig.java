/*
 * BeautiTesterConfig.java
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

package test.dr.app.beauti;

import dr.evomodel.substmodel.aminoacid.AminoAcidModelType;
import dr.evomodel.substmodel.nucleotide.NucModelType;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.util.NexusApplicationImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.CharSetAlignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.*;
import dr.evoxml.util.DateUnitsType;

import java.io.*;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiTester.java,v 1.2 2005/07/11 14:07:25 rambaut Exp $
 */
public class BeautiTesterConfig {

    private PrintWriter scriptWriter;
    BeastGenerator generator;

    private String COMMEND = "beast ";

    public void setCOMMEND(String commend) {
        COMMEND = commend;
    }

    public BeautiTesterConfig() {
    }

    public BeautiOptions createOptions() {
        BeautiOptions beautiOptions = new BeautiOptions();

        beautiOptions.fileNameStem = "";
        beautiOptions.treeFileName.clear();
        beautiOptions.substTreeLog = false;
        beautiOptions.substTreeFileName.clear();

        // MCMC options
        beautiOptions.chainLength = 100;
        beautiOptions.logEvery = 100;
        beautiOptions.echoEvery = 100;
        beautiOptions.burnIn = 10;
        beautiOptions.fileName = null;
        beautiOptions.autoOptimize = true;

        // Data options
        beautiOptions.taxonList = null;

        beautiOptions.datesUnits = DateUnitsType.YEARS;
        beautiOptions.datesDirection = DateUnitsType.FORWARDS;

//        beautiOptions.startingTreeType = StartingTreeType.RANDOM;
//        beautiOptions.fixedTree = false;

        beautiOptions.performTraceAnalysis = false;
        beautiOptions.generateCSV = true;  // until beuati button

        beautiOptions.units = Units.Type.SUBSTITUTIONS;
        beautiOptions.maximumTipHeight = 0.0;

//        beautiOptions.meanSubstitutionRate = 1.0;
//        beautiOptions.fixedSubstitutionRate = true;

        generator = new BeastGenerator(beautiOptions, null);

        return beautiOptions;
    }

    public void closeScriptWriter() {
        scriptWriter.close();
    }

    public void createScriptWriter(String pathAndFile) {
        try {
            scriptWriter = new PrintWriter(new FileWriter(pathAndFile));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void printlnScriptWriter(String line) {
        scriptWriter.println(line);
    }


    public void buildNucModels(String key, BeautiOptions beautiOptions) {
        PartitionSubstitutionModel model = beautiOptions.getPartitionSubstitutionModels().get(0);

        model.setNucSubstitutionModel(NucModelType.HKY);
        buildCodonModels(key + "HKY", beautiOptions);
        model.setNucSubstitutionModel(NucModelType.GTR);
        buildCodonModels(key + "GTR", beautiOptions);
    }

    public void buildCodonModels(String key, BeautiOptions beautiOptions) {
        PartitionSubstitutionModel model = beautiOptions.getPartitionSubstitutionModels().get(0);

        model.setCodonHeteroPattern(null);
        model.setUnlinkedSubstitutionModel(false);
        model.setUnlinkedHeterogeneityModel(false);
        buildHeteroModels(key + "", beautiOptions);

        model.setCodonHeteroPattern("123");
        buildHeteroModels(key + "+C123", beautiOptions);

        model.setUnlinkedSubstitutionModel(true);
        model.setUnlinkedHeterogeneityModel(false);
        buildHeteroModels(key + "+C123^S", beautiOptions);

        model.setUnlinkedSubstitutionModel(false);
        model.setUnlinkedHeterogeneityModel(true);
        buildHeteroModels(key + "+C123^H", beautiOptions);

        model.setUnlinkedSubstitutionModel(true);
        model.setUnlinkedHeterogeneityModel(true);
        buildHeteroModels(key + "+C123^SH", beautiOptions);

        model.setCodonHeteroPattern("112");
        buildHeteroModels(key + "+C112", beautiOptions);

        model.setUnlinkedSubstitutionModel(true);
        model.setUnlinkedHeterogeneityModel(false);
        buildHeteroModels(key + "+C112^S", beautiOptions);

        model.setUnlinkedSubstitutionModel(false);
        model.setUnlinkedHeterogeneityModel(true);
        buildHeteroModels(key + "+C112^H", beautiOptions);

        model.setUnlinkedSubstitutionModel(true);
        model.setUnlinkedHeterogeneityModel(true);
        buildHeteroModels(key + "+C112^SH", beautiOptions);

    }

    public void buildHeteroModels(String key, BeautiOptions beautiOptions) {
        PartitionSubstitutionModel model = beautiOptions.getPartitionSubstitutionModels().get(0);

        model.setGammaHetero(false);
        model.setGammaCategories(4);
        model.setInvarHetero(false);
        buildTreePriorModels(key + "", beautiOptions);

        model.setGammaHetero(true);
        model.setInvarHetero(false);
        buildTreePriorModels(key + "+G", beautiOptions);

        model.setGammaHetero(false);
        model.setInvarHetero(true);
        buildTreePriorModels(key + "+I", beautiOptions);

        model.setGammaHetero(true);
        model.setInvarHetero(true);
        buildTreePriorModels(key + "+GI", beautiOptions);
    }

    public void buildAAModels(String key, BeautiOptions beautiOptions) {
        PartitionSubstitutionModel model = beautiOptions.getPartitionSubstitutionModels().get(0);

        model.setAaSubstitutionModel(AminoAcidModelType.BLOSUM_62);
        buildHeteroModels(key + "BLOSUM62", beautiOptions);

        model.setAaSubstitutionModel(AminoAcidModelType.CP_REV_45);
        buildHeteroModels(key + "CPREV45", beautiOptions);

        model.setAaSubstitutionModel(AminoAcidModelType.DAYHOFF);
        buildHeteroModels(key + "DAYHOFF", beautiOptions);

        model.setAaSubstitutionModel(AminoAcidModelType.JTT);
        buildHeteroModels(key + "JTT", beautiOptions);

        model.setAaSubstitutionModel(AminoAcidModelType.MT_REV_24);
        buildHeteroModels(key + "MTREV24", beautiOptions);

        model.setAaSubstitutionModel(AminoAcidModelType.WAG);
        buildHeteroModels(key + "WAG", beautiOptions);
        
        model.setAaSubstitutionModel(AminoAcidModelType.LG);
        buildHeteroModels(key + "LG", beautiOptions);
    }

    public void buildTreePriorModels(String key, BeautiOptions options) {

//        options.nodeHeightPrior = TreePrior.CONSTANT;
//        buildClockModels(key + "+CP", options);
//
//        options.nodeHeightPrior = TreePrior.EXPONENTIAL;
//        options.parameterization = BeautiOptions.GROWTH_RATE;
//        buildClockModels(key + "+EG", options);
//
//        options.nodeHeightPrior = TreePrior.LOGISTIC;
//        options.parameterization = BeautiOptions.GROWTH_RATE;
//        buildClockModels(key + "+LG", options);
//
//        options.nodeHeightPrior = TreePrior.EXPANSION;
//        options.parameterization = BeautiOptions.GROWTH_RATE;
//        buildClockModels(key + "+XG", options);
//
//        options.nodeHeightPrior = TreePrior.SKYLINE;
//        options.skylineGroupCount = 3;
//        options.skylineModel = BeautiOptions.CONSTANT_SKYLINE;
//        buildClockModels(key + "+SKC", options);
//
//        options.skylineModel = BeautiOptions.LINEAR_SKYLINE;
//        buildClockModels(key + "+SKL", options);

    }

    public void buildClockModels(String key, BeautiOptions options) {
//        options.clockType = ClockType.STRICT_CLOCK;
//        generate(key + "+CLOC", options);
//        options.clockType = ClockType.UNCORRELATED_EXPONENTIAL;
//        generate(key + "+UCED", options);
//        options.clockType = ClockType.UNCORRELATED_LOGNORMAL;
//        generate(key + "+UCLD", options);
    }

    public void generate(String name, BeautiOptions options) {
        options.logFileName = name + ".log";
//        options.treeFileName = name + ".trees";

        System.out.println("Generating: " + name);
        String fileName = name + ".xml";
        try {
            generator.generateXML(new File(fileName));
        } catch (Exception e) {
            e.printStackTrace();
        }

//        scriptWriter.println("beast " + fileName);
        printlnScriptWriter(COMMEND + fileName);
    }

    public void importFromFile(String fileName, BeautiOptions options, boolean translate) {

        TaxonList taxa = null;
        Alignment alignment = null;
        Tree tree = null;
        PartitionSubstitutionModel model = null;
        java.util.List<NexusApplicationImporter.CharSet> charSets = new ArrayList<NexusApplicationImporter.CharSet>();

        try {
            FileReader reader = new FileReader(fileName);

            NexusApplicationImporter importer = new NexusApplicationImporter(reader);

            boolean done = false;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxa != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        taxa = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.CALIBRATION_BLOCK) {
                        if (taxa == null) {
                            throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
                        }

                        importer.parseCalibrationBlock(taxa);

                    } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                        if (taxa == null) {
                            throw new NexusImporter.MissingBlockException("TAXA block must be defined before a CHARACTERS block");
                        }

                        if (alignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        alignment = importer.parseCharactersBlock(options.taxonList);

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        if (alignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = importer.parseDataBlock(options.taxonList);
                        if (taxa == null) {
                            taxa = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        if (taxa == null) {
                            throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a TREES block");
                        }

                        if (tree != null) {
                            throw new NexusImporter.MissingBlockException("TREES block already defined");
                        }

                        Tree[] trees = importer.parseTreesBlock(taxa);
                        if (trees.length > 0) {
                            tree = trees[0];
                        }

                    } else if (block == NexusApplicationImporter.PAUP_BLOCK) {

                        model = importer.parsePAUPBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

                        model = importer.parseMrBayesBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.ASSUMPTIONS_BLOCK) {

                        importer.parseAssumptionsBlock(charSets, null);

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (alignment == null && taxa == null) {
                throw new NexusImporter.MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (FileNotFoundException fnfe) {
            System.err.println("File not found: " + fnfe);
            System.exit(1);

        } catch (Importer.ImportException ime) {
            System.err.println("Error parsing imported file: " + ime);
            System.exit(1);
        } catch (IOException ioex) {
            System.err.println("File I/O Error: " + ioex);
            System.exit(1);
        } catch (Exception ex) {
            System.err.println("Fatal exception: " + ex);
            System.exit(1);
        }

        if (options.taxonList == null) {
            // This is the first partition to be loaded...

            options.taxonList = new Taxa(taxa);

            // check the taxon names for invalid characters
            boolean foundAmp = false;
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                String name = taxa.getTaxon(i).getId();
                if (name.indexOf('&') >= 0) {
                    foundAmp = true;
                }
            }
            if (foundAmp) {
                System.err.println("One or more taxon names include an illegal character ('&').");
                return;
            }

            // make sure they all have dates...
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                if (taxa.getTaxonAttribute(i, "date") == null) {
                    java.util.Date origin = new java.util.Date(0);

                    dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                    taxa.getTaxon(i).setAttribute("date", date);
                }
            }

        } else {
            // This is an additional partition so check it uses the same taxa

            java.util.List<String> oldTaxa = new ArrayList<String>();
            for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                oldTaxa.add(options.taxonList.getTaxon(i).getId());
            }
            java.util.List<String> newTaxa = new ArrayList<String>();
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                newTaxa.add(taxa.getTaxon(i).getId());
            }

            if (!(oldTaxa.containsAll(newTaxa) && oldTaxa.size() == newTaxa.size())) {
                System.err.println("This file contains different taxa from the previously loaded");
                return;
            }
        }

        String fileNameStem = dr.app.util.Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "TRE", "TREE"});

        if (alignment != null) {
            if (translate) {
                alignment = new ConvertAlignment(AminoAcids.INSTANCE, GeneticCode.UNIVERSAL, alignment);
            }

            java.util.List<PartitionData> partitions = new ArrayList<PartitionData>();
            if (charSets != null && charSets.size() > 0) {
                for (NexusApplicationImporter.CharSet charSet : charSets) {
                    partitions.add(new PartitionData(options, charSet.getName(), fileName,
                            new CharSetAlignment(charSet, alignment)));
                }
            } else {
                partitions.add(new PartitionData(options, fileNameStem, fileName, alignment));
            }
            for (PartitionData partition : partitions) {
                if (model != null) {
                    partition.setPartitionSubstitutionModel(model);
//                    options.addPartitionSubstitutionModel(model);
                } else {
                    for (PartitionSubstitutionModel pm : options.getPartitionSubstitutionModels()) {
                        if (pm.getDataType() == alignment.getDataType()) {
                            partition.setPartitionSubstitutionModel(pm);
                        }
                    }
                    if (partition.getPartitionSubstitutionModel() == null) {
                        PartitionSubstitutionModel pm = new PartitionSubstitutionModel(options, "default", partition);
                        partition.setPartitionSubstitutionModel(pm);
//                        options.addPartitionSubstitutionModel(pm);
                    }
                }
                options.dataPartitions.add(partition);
            }
        }

        calculateHeights(options);
    }

    private void calculateHeights(BeautiOptions options) {

        options.maximumTipHeight = 0.0;
        if (!options.hasData()) return;

        dr.evolution.util.Date mostRecent = null;
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Date date = options.taxonList.getTaxon(i).getDate();
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                Date date = options.taxonList.getTaxon(i).getDate();
                if (date != null) {
                    double height = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    if (height > options.maximumTipHeight) options.maximumTipHeight = height;
                }
            }
        }
    }


}

