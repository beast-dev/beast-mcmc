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
package dr.app.oldbeauti;

import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import org.jdom.Document;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import java.io.*;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: BeautiFrame.java,v 1.22 2006/09/09 16:07:06 rambaut Exp $
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

		try {
			FileReader reader = new FileReader(file);

			NexusApplicationImporter importer = new NexusApplicationImporter(reader);

			boolean done = false;

			beautiOptions.originalAlignment = null;
			beautiOptions.alignment = null;
			beautiOptions.tree = null;
			beautiOptions.taxonList = null;

			while (!done) {
				try {

					NexusImporter.NexusBlock block = importer.findNextBlock();

					if (block == NexusImporter.TAXA_BLOCK) {

						if (beautiOptions.taxonList != null) {
							throw new NexusImporter.MissingBlockException("TAXA block already defined");
						}

						beautiOptions.taxonList = importer.parseTaxaBlock();

					} else if (block == NexusImporter.CALIBRATION_BLOCK) {
						if (beautiOptions.taxonList == null) {
							throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
						}

						importer.parseCalibrationBlock(beautiOptions.taxonList);

					} else if (block == NexusImporter.CHARACTERS_BLOCK) {

						if (beautiOptions.taxonList == null) {
							throw new NexusImporter.MissingBlockException("TAXA block must be defined before a CHARACTERS block");
						}

						if (beautiOptions.originalAlignment != null) {
							throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
						}

						beautiOptions.originalAlignment = (SimpleAlignment)importer.parseCharactersBlock(beautiOptions.taxonList);

					} else if (block == NexusImporter.DATA_BLOCK) {

						if (beautiOptions.originalAlignment != null) {
							throw new NexusImporter.MissingBlockException("CHARACTERS or DATA block already defined");
						}

						// A data block doesn't need a taxon block before it
						// but if one exists then it will use it.
						beautiOptions.originalAlignment = (SimpleAlignment)importer.parseDataBlock(beautiOptions.taxonList);
						if (beautiOptions.taxonList == null) {
							beautiOptions.taxonList = beautiOptions.originalAlignment;
						}

					} else if (block == NexusImporter.TREES_BLOCK) {

						if (beautiOptions.taxonList == null) {
							throw new NexusImporter.MissingBlockException("TAXA or DATA block must be defined before a TREES block");
						}

						if (beautiOptions.tree != null) {
							throw new NexusImporter.MissingBlockException("TREES block already defined");
						}

						Tree[] trees = importer.parseTreesBlock(beautiOptions.taxonList);
						if (trees.length > 0) {
							beautiOptions.tree = trees[0];
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

			if (beautiOptions.originalAlignment == null) {
				throw new NexusImporter.MissingBlockException("DATA or CHARACTERS block is missing");
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

		// check the taxon names for invalid characters
		boolean foundAmp = false;
		for (int i = 0; i < beautiOptions.originalAlignment.getTaxonCount(); i++) {
			String name = beautiOptions.originalAlignment.getTaxon(i).getId();
			if (name.indexOf('&') >= 0) {
				foundAmp = true;
			}
		}
		if (foundAmp) {
			System.err.println("One or more taxon names include an illegal character ('&').\n" +
					"These characters will prevent BEAST from reading the resulting XML file.\n\n" +
					"Please edit the taxon name(s) before generating the BEAST file.");
		}


		// make sure they all have dates...
		for (int i = 0; i < beautiOptions.originalAlignment.getTaxonCount(); i++) {
			if (beautiOptions.originalAlignment.getTaxonAttribute(i, "date") == null) {
				java.util.Date origin = new java.util.Date(0);

				dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
				beautiOptions.originalAlignment.getTaxon(i).setAttribute("date", date);
			}
		}

		beautiOptions.fileNameStem = dr.app.util.Utils.trimExtensions(file.getName(),
				new String[] {"nex", "NEX", "tre", "TRE", "nexus", "NEXUS"});

		beautiOptions.alignment = beautiOptions.originalAlignment;
		beautiOptions.alignmentReset = true;
		if (beautiOptions.alignment != null) {
			Patterns patterns = new Patterns(beautiOptions.alignment);
			DistanceMatrix distances = new JukesCantorDistanceMatrix(patterns);
			beautiOptions.meanDistance = distances.getMeanDistance();
		}

		return true;
	}

	private void generate(File file) throws IOException {
		FileWriter fw = new FileWriter(file);
		beautiOptions.generateXML(fw);
		fw.close();
	}



}
