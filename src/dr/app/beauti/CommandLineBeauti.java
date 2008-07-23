/*
 * BeautiMacFileMenuFactory.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
package dr.app.beauti;

import dr.app.beauti.options.DataPartition;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import org.jdom.Document;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import java.io.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
 */
public class CommandLineBeauti {
    private BeastGenerator beautiOptions = new BeastGenerator();

    public CommandLineBeauti(String inputFileName, String templateFileName, String outputFileName) {

        try {
            if (!importFromFile(new File(inputFileName))) {
                return;
            }
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: Input file not found");
            return;
        } catch (IOException ioe) {
            System.err.println("Error reading input file: " + ioe.getMessage());
            return;
        }

        try {
            if (!readFromFile(new File(templateFileName))) {
                return;
            }
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: Template file not found");
            return;
        } catch (IOException ioe) {
            System.err.println("Error reading template file: " + ioe.getMessage());
            return;
        }

        beautiOptions.guessDates();

        try {
            generate(new File(outputFileName));

        } catch (IOException ioe) {
            System.err.println("Unable to generate file: " + ioe.getMessage());
            return;
        }
    }

    private boolean readFromFile(File file) throws FileNotFoundException, IOException {
        try {
            SAXBuilder parser = new SAXBuilder();
            Document doc = parser.build(file);
            beautiOptions.parse(doc);

        } catch (dr.xml.XMLParseException xpe) {
            System.err.println("Error reading file: This may not be a BEAUti Template file");
            System.err.println(xpe.getMessage());
            return false;
        } catch (JDOMException e) {
            System.err.println("Unable to open file: This may not be a BEAUti Template file");
            System.err.println(e.getMessage());
            return false;
        }
        return true;
    }

    private boolean importFromFile(File file) throws FileNotFoundException, IOException {

        Alignment alignment = null;
        Tree tree = null;
        TaxonList taxa = null;

        try {
            FileReader reader = new FileReader(file);

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

                        importer.parseCalibrationBlock(beautiOptions.taxonList);

                    } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                        if (taxa == null) {
                            throw new NexusImporter.MissingBlockException("TAXA block must be defined before a CHARACTERS block");
                        }

                        if (alignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        alignment = (SimpleAlignment) importer.parseCharactersBlock(beautiOptions.taxonList);

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        if (alignment != null) {
                            throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = (SimpleAlignment) importer.parseDataBlock(beautiOptions.taxonList);
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

/*					} else if (block == NexusApplicationImporter.PAUP_BLOCK) {

						importer.parsePAUPBlock(beautiOptions);

					} else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

						importer.parseMrBayesBlock(beautiOptions);

					} else if (block == NexusApplicationImporter.RHINO_BLOCK) {

						importer.parseRhinoBlock(beautiOptions);
*/
                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            if (alignment == null && taxa == null) {
                throw new NexusImporter.MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (Importer.ImportException ime) {
            System.err.println("Error parsing imported file: " + ime);
            return false;
        } catch (IOException ioex) {
            System.err.println("File I/O Error: " + ioex);
            return false;
        } catch (Exception ex) {
            System.err.println("Fatal exception: " + ex);
            return false;
        }

        if (beautiOptions.taxonList == null) {
            // This is the first partition to be loaded...

            beautiOptions.taxonList = taxa;

            // check the taxon names for invalid characters
            boolean foundAmp = false;
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                String name = taxa.getTaxon(i).getId();
                if (name.indexOf('&') >= 0) {
                    foundAmp = true;
                }
            }
            if (foundAmp) {
                System.err.println("One or more taxon names include an illegal character ('&').\n" +
                        "These characters will prevent BEAST from reading the resulting XML file.\n\n" +
                        "Please edit the taxon name(s) before reloading the data file.");
                return false;
            }

            // make sure they all have dates...
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                if (taxa.getTaxonAttribute(i, "date") == null) {
                    java.util.Date origin = new java.util.Date(0);

                    dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                    taxa.getTaxon(i).setAttribute("date", date);
                }
            }

            beautiOptions.fileNameStem = dr.app.util.Utils.trimExtensions(file.getName(),
                    new String[]{"nex", "NEX", "tre", "TRE", "nexus", "NEXUS"});

            if (alignment != null) {
                DataPartition partition = new DataPartition(beautiOptions.fileNameStem, file.getName(), alignment);
                beautiOptions.dataPartitions.add(partition);
                beautiOptions.dataType = alignment.getDataType();

                Patterns patterns = new Patterns(alignment);
                DistanceMatrix distances = new JukesCantorDistanceMatrix(patterns);
                beautiOptions.meanDistance = distances.getMeanDistance();

            } else {
                beautiOptions.meanDistance = 0.0;
            }
        } else {
            // This is an additional partition so check it uses the same taxa
        }

        return true;
    }

    private void generate(File file) throws IOException {
        FileWriter fw = new FileWriter(file);
        beautiOptions.generateXML(fw);
        fw.close();
    }


}
