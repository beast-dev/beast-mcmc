/*
 * BeautiTester.java
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

import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.TimeScale;
import dr.evolution.util.Units;

import java.io.*;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: BeautiTester.java,v 1.2 2005/07/11 14:07:25 rambaut Exp $
 */
public class BeautiTester {

    PrintWriter scriptWriter;

    public BeautiTester() {
        BeastGenerator beautiOptions = createOptions();

        try {
            scriptWriter = new PrintWriter(new FileWriter("tests/run_script.sh"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        importFromFile("examples/Primates.nex", beautiOptions);

        buildNucModels("tests/pri_", beautiOptions);
        buildAAModels("tests/pri_", beautiOptions);

        importFromFile("examples/Dengue4.env.nex", beautiOptions);
        beautiOptions.fixedSubstitutionRate = false;

        buildNucModels("tests/den_", beautiOptions);
        buildAAModels("tests/den_", beautiOptions);

        scriptWriter.close();
    }

    public BeastGenerator createOptions() {

        BeastGenerator beautiOptions = new BeastGenerator();

        beautiOptions.fileNameStem = "";
        beautiOptions.substTreeLog = false;
        beautiOptions.substTreeFileName = null;

        // MCMC options
        beautiOptions.chainLength = 100;
        beautiOptions.logEvery = 100;
        beautiOptions.echoEvery = 100;
        beautiOptions.burnIn = 10;
        beautiOptions.fileName = null;
        beautiOptions.autoOptimize = true;

        // Data options
        beautiOptions.taxonList = null;
        beautiOptions.originalAlignment = null;
        beautiOptions.alignment = null;
        beautiOptions.tree = null;

        beautiOptions.datesUnits = BeautiOptions.YEARS;
        beautiOptions.datesDirection = BeautiOptions.FORWARDS;

        beautiOptions.userTree = false;
        beautiOptions.fixedTree = false;

        beautiOptions.performTraceAnalysis = false;
        beautiOptions.generateCSV = true;  // until beuati button

        beautiOptions.units = Units.Type.SUBSTITUTIONS;
        beautiOptions.maximumTipHeight = 0.0;

        beautiOptions.meanSubstitutionRate = 1.0;
        beautiOptions.fixedSubstitutionRate = true;
        return beautiOptions;
    }

    public void buildNucModels(String key, BeastGenerator beautiOptions) {
        beautiOptions.alignment = beautiOptions.originalAlignment;

        beautiOptions.nucSubstitutionModel = BeautiOptions.HKY;
        buildCodonModels(key+"HKY", beautiOptions);
        beautiOptions.nucSubstitutionModel = BeautiOptions.GTR;
        buildCodonModels(key+"GTR", beautiOptions);
    }

    public void buildCodonModels(String key, BeastGenerator beautiOptions) {
        beautiOptions.codonHeteroPattern = null;
        beautiOptions.unlinkedSubstitutionModel = false;
        beautiOptions.unlinkedHeterogeneityModel = false;
        buildHeteroModels(key+"", beautiOptions);

        beautiOptions.codonHeteroPattern = "123";
        buildHeteroModels(key+"+C123", beautiOptions);

        beautiOptions.unlinkedSubstitutionModel = true;
        beautiOptions.unlinkedHeterogeneityModel = false;
        buildHeteroModels(key+"+C123^S", beautiOptions);

        beautiOptions.unlinkedSubstitutionModel = false;
        beautiOptions.unlinkedHeterogeneityModel = true;
        buildHeteroModels(key+"+C123^H", beautiOptions);

        beautiOptions.unlinkedSubstitutionModel = true;
        beautiOptions.unlinkedHeterogeneityModel = true;
        buildHeteroModels(key+"+C123^SH", beautiOptions);

        beautiOptions.codonHeteroPattern = "112";
	    buildHeteroModels(key+"+C112", beautiOptions);

	    beautiOptions.unlinkedSubstitutionModel = true;
	    beautiOptions.unlinkedHeterogeneityModel = false;
	    buildHeteroModels(key+"+C112^S", beautiOptions);

	    beautiOptions.unlinkedSubstitutionModel = false;
	    beautiOptions.unlinkedHeterogeneityModel = true;
	    buildHeteroModels(key+"+C112^H", beautiOptions);

	    beautiOptions.unlinkedSubstitutionModel = true;
	    beautiOptions.unlinkedHeterogeneityModel = true;
	    buildHeteroModels(key+"+C112^SH", beautiOptions);

    }

    public void buildHeteroModels(String key, BeastGenerator beautiOptions) {

        beautiOptions.gammaHetero = false;
        beautiOptions.gammaCategories = 4;
        beautiOptions.invarHetero = false;
        buildTreePriorModels(key+"", beautiOptions);

        beautiOptions.gammaHetero = true;
        beautiOptions.invarHetero = false;
        buildTreePriorModels(key+"+G", beautiOptions);

        beautiOptions.gammaHetero = false;
        beautiOptions.invarHetero = true;
        buildTreePriorModels(key+"+I", beautiOptions);

        beautiOptions.gammaHetero = true;
        beautiOptions.invarHetero = true;
        buildTreePriorModels(key+"+GI", beautiOptions);
    }

    public void buildAAModels(String key, BeastGenerator beautiOptions) {

        beautiOptions.alignment = new ConvertAlignment(AminoAcids.INSTANCE, GeneticCode.UNIVERSAL, beautiOptions.originalAlignment);
        /*
        beautiOptions.aaSubstitutionModel = BeautiOptions.BLOSUM_62;
        buildHeteroModels(key+"BLOSUM62", beautiOptions);

        beautiOptions.aaSubstitutionModel = BeautiOptions.CP_REV_45;
        buildHeteroModels(key+"CPREV45", beautiOptions);

        beautiOptions.aaSubstitutionModel = BeautiOptions.DAYHOFF;
        buildHeteroModels(key+"DAYHOFF", beautiOptions);

        beautiOptions.aaSubstitutionModel = BeautiOptions.JTT;
        buildHeteroModels(key+"JTT", beautiOptions);

        beautiOptions.aaSubstitutionModel = BeautiOptions.MT_REV_24;
        buildHeteroModels(key+"MTREV24", beautiOptions);
        */
        beautiOptions.aaSubstitutionModel = BeautiOptions.WAG;
        buildHeteroModels(key+"WAG", beautiOptions);
    }

    public void buildTreePriorModels(String key, BeastGenerator beautiOptions) {

        beautiOptions.nodeHeightPrior = BeautiOptions.CONSTANT;
        buildClockModels(key+"+CP", beautiOptions);

        beautiOptions.nodeHeightPrior = BeautiOptions.EXPONENTIAL;
        beautiOptions.parameterization = BeautiOptions.GROWTH_RATE;
        buildClockModels(key+"+EG", beautiOptions);

        beautiOptions.nodeHeightPrior = BeautiOptions.LOGISTIC;
        beautiOptions.parameterization = BeautiOptions.GROWTH_RATE;
        buildClockModels(key+"+LG", beautiOptions);

        beautiOptions.nodeHeightPrior = BeautiOptions.EXPANSION;
        beautiOptions.parameterization = BeautiOptions.GROWTH_RATE;
        buildClockModels(key+"+XG", beautiOptions);

        beautiOptions.nodeHeightPrior = BeautiOptions.SKYLINE;
        beautiOptions.skylineGroupCount = 3;
        beautiOptions.skylineModel = BeautiOptions.CONSTANT_SKYLINE;
        buildClockModels(key+"+SKC", beautiOptions);

        beautiOptions.skylineModel = BeautiOptions.LINEAR_SKYLINE;
        buildClockModels(key+"+SKL", beautiOptions);

    }

    public void buildClockModels(String key, BeastGenerator beautiOptions) {
        beautiOptions.clockModel = BeautiOptions.STRICT_CLOCK;
        generate(key+"+CLOC", beautiOptions);
        beautiOptions.clockModel = BeautiOptions.UNCORRELATED_EXPONENTIAL;
        generate(key+"+UCED", beautiOptions);
        beautiOptions.clockModel = BeautiOptions.UNCORRELATED_LOGNORMAL;
        generate(key+"+UCLD", beautiOptions);
    }

    public void generate(String name, BeastGenerator beautiOptions) {
        beautiOptions.logFileName = name + ".log";
        beautiOptions.treeFileName = name + ".trees";

        System.out.println("Generating: " + name);
        String fileName = name + ".xml";
        try {
            FileWriter fw = new FileWriter(fileName);
            beautiOptions.generateXML(fw);
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        scriptWriter.println("beast " + fileName);
    }

    protected void importFromFile(String fileName, BeastGenerator beautiOptions) {

        try {
            FileReader reader = new FileReader(fileName);

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

        // make sure they all have dates...
        for (int i = 0; i < beautiOptions.originalAlignment.getTaxonCount(); i++) {
            if (beautiOptions.originalAlignment.getTaxonAttribute(i, "date") == null) {
                java.util.Date origin = new java.util.Date(0);

                dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                beautiOptions.originalAlignment.getTaxon(i).setAttribute("date", date);
            }
        }

        beautiOptions.alignment = beautiOptions.originalAlignment;
        beautiOptions.taxonList = beautiOptions.originalAlignment;

        calculateHeights(beautiOptions);
    }

    private void calculateHeights(BeautiOptions options) {

        options.maximumTipHeight = 0.0;
        if (options.alignment == null) return;

        dr.evolution.util.Date mostRecent = null;
        for (int i = 0; i < options.alignment.getSequenceCount(); i++) {
            Date date = options.alignment.getTaxon(i).getDate();
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < options.alignment.getSequenceCount(); i++) {
                Date date = options.alignment.getTaxon(i).getDate();
                if (date != null) {
                    double height = timeScale.convertTime(date.getTimeValue(), date) - time0;
                    if (height > options.maximumTipHeight) options.maximumTipHeight = height;
                }
            }
        }
    }

	//Main method
	public static void main(String[] args) {

		new BeautiTester();
	}
}

