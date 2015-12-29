/*
 * BeagleSequenceSimulatorConsoleApp.java
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

package dr.app.bss;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.util.Arguments;
import dr.app.util.Arguments.ArgumentException;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.math.MathUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class BeagleSequenceSimulatorConsoleApp {

    private boolean VERBOSE = true;

    private Arguments arguments;
    private PartitionDataList dataList;

    private static final String SPLIT_PARTITION = ":";

    private static final String HELP = "help";

    private static final String TREE_FILE = "treeFile";
    private static final String TAXA_SET = "taxaSet";

    private static final String DEMOGRAPHIC_MODEL = "demographicModel";
    private static final String NO_DEMOGRAPHIC_MODEL = "noModel";
    private static final String CONSTANT_POPULATION = "constantPopulation";
    private static final String CONSTANT_POPULATION_PARAMETER_VALUES = "constantPopulationParameterValues";
    private static final String EXPONENTIAL_GROWTH_RATE = "exponentialGrowthRate";
    private static final String EXPONENTIAL_GROWTH_RATE_PARAMETER_VALUES = "exponentialGrowthRateParameterValues";
    private static final String EXPONENTIAL_DOUBLING_TIME = "exponentialDoublingTime";
    private static final String EXPONENTIAL_GROWTH_DOUBLING_TIME_PARAMETER_VALUES = "exponentialDoublingTimeParameterValues";

    private static final String BRANCH_SUBSTITUTION_MODEL = "branchSubstitutionModel";
    private static final String HKY = PartitionData.substitutionModels[0];
    private static final String HKY_SUBSTITUTION_PARAMETER_VALUES = "HKYsubstitutionParameterValues";
    private static final String GTR = PartitionData.substitutionModels[1];
    private static final String GTR_SUBSTITUTION_PARAMETER_VALUES = "GTRsubstitutionParameterValues";
    private static final String TN93 = PartitionData.substitutionModels[2];
    private static final String TN93_SUBSTITUTION_PARAMETER_VALUES = "TN93substitutionParameterValues";
    private static final String GY94_CODON_MODEL = PartitionData.substitutionModels[3];
    private static final String GY94_SUBSTITUTION_PARAMETER_VALUES = "GY94substitutionParameterValues";

    private static final String SITE_RATE_MODEL = "siteRateModel";
    private static final String NO_SITE_RATE_MODEL = "noModel";
    private static final String GAMMA_SITE_RATE_MODEL = "gammaSiteRateModel";
    private static final String GAMMA_SITE_RATE_MODEL_PARAMETER_VALUES = "gammaSiteRateModelParameterValues";

    private static final String CLOCK_RATE_MODEL = "clockRateModel";
    private static final String STRICT_CLOCK = "strictClock";
    private static final String STRICT_CLOCK_PARAMETER_VALUES = "strictClockParameterValues";
    private static final String LOGNORMAL_RELAXED_CLOCK = "lognormalRelaxedClock";
    private static final String LOGNORMAL_RELAXED_CLOCK_PARAMETER_VALUES = "lognormalRelaxedClockParameterValues";
    
    private static final String LRC_PARAMETERS_IN_REAL_SPACE = "lrcParametersInRealSpace";
    private static final String TRUE = "true";
    private static final String FALSE = "false";
    
    private static final String EXPONENTIAL_RELAXED_CLOCK = "exponentialRelaxedClock";
    private static final String EXPONENTIAL_RELAXED_CLOCK_PARAMETER_VALUES = "exponentialRelaxedClockParameterValues";
    private static final String INVERSE_GAUSSIAN_RELAXED_CLOCK = "inverseGaussianRelaxedClock";
    private static final String INVERSE_GAUSSIAN_RELAXED_CLOCK_PARAMETER_VALUES = "inverseGaussianRelaxedClockParameterValues";

    private static final String BASE_FREQUENCIES = "baseFrequencies";
    private static final String NUCLEOTIDE_FREQUENCIES = "nucleotideFrequencies";
    private static final String NUCLEOTIDE_FREQUENCY_PARAMETER_VALUES = "nucleotideFrequencyParameterValues";
    private static final String CODON_FREQUENCIES = "codonFrequencies";
    private static final String CODON_FREQUENCY_PARAMETER_VALUES = "codonFrequencyParameterValues";

    private static final String FROM = "from";
    private static final String TO = "to";
    private static final String EVERY = "every";

    private static final String ROOT_SEQUENCE = "rootSequence";

    public BeagleSequenceSimulatorConsoleApp() {

//        data = new PartitionData();
        dataList = new PartitionDataList();

        // //////////////////
        // ---DEFINITION---//
        // //////////////////

        arguments = new Arguments(
                new Arguments.Option[]{

                        new Arguments.Option(HELP,
                                "print this information and exit"),

                        new Arguments.StringOption(TREE_FILE, "tree file",
                                "specify tree topology"),

                        new Arguments.StringOption(TAXA_SET, "taxa set",
                                "specify taxa set"),

                        new Arguments.IntegerOption(FROM,
                                "specify 'from' attribute"),
                        new Arguments.IntegerOption(TO,
                                "specify 'to' attribute"),
                        new Arguments.IntegerOption(EVERY,
                                "specify 'every' attribute"),

                        new Arguments.StringOption(DEMOGRAPHIC_MODEL,
                                new String[]{NO_DEMOGRAPHIC_MODEL, //
                                        CONSTANT_POPULATION, //
                                        EXPONENTIAL_GROWTH_RATE, //
                                        EXPONENTIAL_DOUBLING_TIME //
                                }, false, "specify demographic model"),

                        new Arguments.RealArrayOption(CONSTANT_POPULATION_PARAMETER_VALUES, 1, "specify constant population model parameter values"),
                        new Arguments.RealArrayOption(EXPONENTIAL_GROWTH_RATE_PARAMETER_VALUES, 2, "specify exponential growth (growth rate) population model parameter values"),
                        new Arguments.RealArrayOption(EXPONENTIAL_GROWTH_DOUBLING_TIME_PARAMETER_VALUES, 2, "specify exponential growth (doubling time) population model parameter values"),

                        new Arguments.StringOption(BRANCH_SUBSTITUTION_MODEL,
                                new String[]{HKY, //
                                        GTR, //
                                        TN93, //
                                        GY94_CODON_MODEL //
                                }, false, "specify substitution model"),

                        new Arguments.RealArrayOption(HKY_SUBSTITUTION_PARAMETER_VALUES, 1, "specify HKY substitution model parameter values"),
                        new Arguments.RealArrayOption(GTR_SUBSTITUTION_PARAMETER_VALUES, 6, "specify GTR substitution model parameter values"),
                        new Arguments.RealArrayOption(TN93_SUBSTITUTION_PARAMETER_VALUES, 2, "specify TN93 substitution model parameter values"),
                        new Arguments.RealArrayOption(GY94_SUBSTITUTION_PARAMETER_VALUES, 2, "specify GY94 substitution model parameter values"),

                        new Arguments.StringOption(SITE_RATE_MODEL,
                                new String[]{NO_SITE_RATE_MODEL, //
                                        GAMMA_SITE_RATE_MODEL, //
                                }, false, "specify site rate model"),

                        new Arguments.RealArrayOption(
                                GAMMA_SITE_RATE_MODEL_PARAMETER_VALUES, 3, "specify Gamma Site Rate Model parameter values"),

                        new Arguments.StringOption(CLOCK_RATE_MODEL,
                                new String[]{STRICT_CLOCK, //
                                        LOGNORMAL_RELAXED_CLOCK, //
                                        EXPONENTIAL_RELAXED_CLOCK, //
                                        INVERSE_GAUSSIAN_RELAXED_CLOCK},
                                false, "specify clock rate model"),

                        new Arguments.RealArrayOption(STRICT_CLOCK_PARAMETER_VALUES, 1, "specify Strict Clock parameter values"),
                        new Arguments.RealArrayOption(LOGNORMAL_RELAXED_CLOCK_PARAMETER_VALUES, 3, "specify Lognormal Relaxed Clock parameter values"),
                        
//                        new Arguments.StringOption(LRC_PARAMETERS_IN_REAL_SPACE),
                        
                        new Arguments.StringOption(LRC_PARAMETERS_IN_REAL_SPACE,
                                new String[]{TRUE, //
                                        FALSE //
                               },
                                false, "Are the Lognormal Relaxed Clock parameter values in the real or log space"),
                        
                        new Arguments.RealArrayOption(EXPONENTIAL_RELAXED_CLOCK_PARAMETER_VALUES, 2, "specify Exponential Relaxed Clock parameter values"),
                        new Arguments.RealArrayOption(INVERSE_GAUSSIAN_RELAXED_CLOCK_PARAMETER_VALUES, 3, "specify Inverse Gaussian Relaxed Clock parameter values"),

                        new Arguments.StringOption(BASE_FREQUENCIES,
                                new String[]{NUCLEOTIDE_FREQUENCIES, //
                                        CODON_FREQUENCIES, //
                                }, false, "specify frequency model"),

                        new Arguments.RealArrayOption(NUCLEOTIDE_FREQUENCY_PARAMETER_VALUES, 4, "specify nucleotide frequency parameter values"),
                        new Arguments.RealArrayOption(CODON_FREQUENCY_PARAMETER_VALUES, 61, "specify codon frequency parameter values"),

                        new Arguments.StringOption(ROOT_SEQUENCE, "ancestral sequence",
                                "specify ancestral sequence"),

                });

    }// END: constructor

    public void simulate(String[] args) {

        try {

            // /////////////////////////////////
            // ---SPLIT PARTITION ARGUMENTS---//
            // /////////////////////////////////

            int from = 0;
            int to = 0;
            ArrayList<String[]> argsList = new ArrayList<String[]>();
            for (String arg : args) {

                if (arg.equalsIgnoreCase(SPLIT_PARTITION)) {
                    argsList.add(Arrays.copyOfRange(args, from, to));
                    from = to + 1;
                }// END: split check

                to++;
            }// END: args loop

            if (args[0].contains(HELP)) {

                gracefullyExit(null);

            } else if (argsList.size() == 0) {

                gracefullyExit("Empty or incorrect arguments list.");

            }// END: failed split check

            String[] leftoverArguments = Arrays.copyOfRange(args, from, args.length);
            if (leftoverArguments.length > 2) {
                gracefullyExit("Unrecognized option " + leftoverArguments[2]);
            }

            ArrayList<Partition> partitionsList = new ArrayList<Partition>();
            for (String partitionArgs[] : argsList) {

                // /////////////
                // ---PARSE---//
                // /////////////

                arguments.parseArguments(partitionArgs);

                // ///////////////////
                // ---INTERROGATE---//
                // ///////////////////

                PartitionData data = new PartitionData();
                String option = null;
                double[] values = null;

                if (partitionArgs.length == 0 || arguments.hasOption(HELP)) {

                    gracefullyExit(null);

                }// END: HELP option check

                // Tree / Taxa
                if (arguments.hasOption(TREE_FILE)) {

                    File treeFile = new File(arguments.getStringOption(TREE_FILE));
                    Tree tree = Utils.importTreeFromFile(treeFile);
                    data.record = new TreesTableRecord(treeFile.getName(), tree);

                } else if (arguments.hasOption(TAXA_SET)) {

                    File taxaFile = new File(arguments.getStringOption(TAXA_SET));
                    Taxa taxa = Utils.importTaxaFromFile(taxaFile);
                    data.record = new TreesTableRecord(taxaFile.getName(), taxa);

                } else {

                    throw new RuntimeException("Tree file / Taxa set not specified.");

                }// END: Tree / Taxa option check

                // Demographic Model
                if (arguments.hasOption(DEMOGRAPHIC_MODEL)) {

                    option = arguments.getStringOption(DEMOGRAPHIC_MODEL);

                    if (option.equalsIgnoreCase(NO_DEMOGRAPHIC_MODEL)) {

                        int index = 0;
                        data.demographicModelIndex = index;

                    } else if (option.equalsIgnoreCase(CONSTANT_POPULATION)) {

                        int index = 1;
                        data.demographicModelIndex = index;

                        if (arguments
                                .hasOption(CONSTANT_POPULATION_PARAMETER_VALUES)) {

                            values = arguments
                                    .getRealArrayOption(CONSTANT_POPULATION_PARAMETER_VALUES);
                            parseDemographicValues(index, values, data);

                        }

                    } else if (option.equalsIgnoreCase(EXPONENTIAL_GROWTH_RATE)) {

                        int index = 2;
                        data.demographicModelIndex = index;

                        if (arguments
                                .hasOption(EXPONENTIAL_GROWTH_RATE_PARAMETER_VALUES)) {

                            values = arguments
                                    .getRealArrayOption(EXPONENTIAL_GROWTH_RATE_PARAMETER_VALUES);
                            parseDemographicValues(index, values, data);

                        }

                    } else if (option
                            .equalsIgnoreCase(EXPONENTIAL_DOUBLING_TIME)) {

                        int index = 3;
                        data.demographicModelIndex = index;

                        if (arguments
                                .hasOption(EXPONENTIAL_GROWTH_DOUBLING_TIME_PARAMETER_VALUES)) {

                            values = arguments
                                    .getRealArrayOption(EXPONENTIAL_GROWTH_DOUBLING_TIME_PARAMETER_VALUES);
                            parseDemographicValues(index, values, data);

                        }

                    } else {
                        gracefullyExit("Unrecognized option.");
                    }

                }// END: DEMOGRAPHIC_MODEL option check

                // Branch Substitution Model
                if (arguments.hasOption(BRANCH_SUBSTITUTION_MODEL)) {

                    option = arguments
                            .getStringOption(BRANCH_SUBSTITUTION_MODEL);

                    if (option.equalsIgnoreCase(HKY)) {

                        int index = 0;
                        data.substitutionModelIndex = index;

                        if (arguments
                                .hasOption(HKY_SUBSTITUTION_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(HKY_SUBSTITUTION_PARAMETER_VALUES);
                            parseSubstitutionValues(index, values, data);
                        }

                    } else if (option.equalsIgnoreCase(GTR)) {

                        int index = 1;
                        data.substitutionModelIndex = index;

                        if (arguments
                                .hasOption(GTR_SUBSTITUTION_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(GTR_SUBSTITUTION_PARAMETER_VALUES);
                            parseSubstitutionValues(index, values, data);
                        }

                    } else if (option.equalsIgnoreCase(TN93)) {

                        int index = 2;
                        data.substitutionModelIndex = index;

                        if (arguments
                                .hasOption(TN93_SUBSTITUTION_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(TN93_SUBSTITUTION_PARAMETER_VALUES);
                            parseSubstitutionValues(index, values, data);
                        }

                    } else if (option.equalsIgnoreCase(GY94_CODON_MODEL)) {

                        int index = 3;
                        data.substitutionModelIndex = index;

                        if (arguments
                                .hasOption(GY94_SUBSTITUTION_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(GY94_SUBSTITUTION_PARAMETER_VALUES);
                            parseSubstitutionValues(index, values, data);
                        }

                    } else {
                        gracefullyExit("Unrecognized option.");
                    }

                }// END: BRANCH_SUBSTITUTION_MODEL option check

                // Site Rate Model
                if (arguments.hasOption(SITE_RATE_MODEL)) {

                    option = arguments.getStringOption(SITE_RATE_MODEL);

                    if (option.equalsIgnoreCase(NO_SITE_RATE_MODEL)) {

                        int index = 0;
                        data.siteRateModelIndex = index;

                    } else if (option.equalsIgnoreCase(GAMMA_SITE_RATE_MODEL)) {

                        int index = 1;
                        data.siteRateModelIndex = index;

                        if (arguments
                                .hasOption(GAMMA_SITE_RATE_MODEL_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(GAMMA_SITE_RATE_MODEL_PARAMETER_VALUES);
                            parseSiteRateValues(index, values, data);
                        }

                    } else {
                        gracefullyExit("Unrecognized option.");
                    }

                }// END: SITE_RATE_MODEL option check

                // Clock Rate Model
                if (arguments.hasOption(CLOCK_RATE_MODEL)) {

                    option = arguments.getStringOption(CLOCK_RATE_MODEL);

                    if (option.equalsIgnoreCase(STRICT_CLOCK)) {

                        int index = 0;
                        data.clockModelIndex = index;
                        if (arguments.hasOption(STRICT_CLOCK_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(STRICT_CLOCK_PARAMETER_VALUES);
                            parseClockValues(index, values, data);
                        }

                    } else if (option.equalsIgnoreCase(LOGNORMAL_RELAXED_CLOCK)) {

                        int index = 1;
                        data.clockModelIndex = index;
                        
                        if (arguments.hasOption(LOGNORMAL_RELAXED_CLOCK_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(LOGNORMAL_RELAXED_CLOCK_PARAMETER_VALUES);
                            parseClockValues(index, values, data);
                        
                        }

                        // set to true/false
						if (arguments.hasOption(LRC_PARAMETERS_IN_REAL_SPACE)) {

							boolean lrcParametersInRealSpace = (arguments
									.getStringOption(
											LRC_PARAMETERS_IN_REAL_SPACE)
									.equalsIgnoreCase(TRUE) ? true : false);
							data.lrcParametersInRealSpace = lrcParametersInRealSpace;

						}
                        
                        
                    } else if (option
                            .equalsIgnoreCase(EXPONENTIAL_RELAXED_CLOCK)) {

                        int index = 2;
                        data.clockModelIndex = index;
                        if (arguments
                                .hasOption(EXPONENTIAL_RELAXED_CLOCK_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(EXPONENTIAL_RELAXED_CLOCK_PARAMETER_VALUES);
                            parseClockValues(index, values, data);
                        }

                    } else if (option
                            .equalsIgnoreCase(INVERSE_GAUSSIAN_RELAXED_CLOCK)) {

                        int index = 3;
                        data.clockModelIndex = index;
                        if (arguments
                                .hasOption(INVERSE_GAUSSIAN_RELAXED_CLOCK_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(INVERSE_GAUSSIAN_RELAXED_CLOCK_PARAMETER_VALUES);
                            parseClockValues(index, values, data);
                        }

                    } else {
                        gracefullyExit("Unrecognized option.");
                    }

                }// END: CLOCK_RATE_MODEL option check

                // Frequency Model
                if (arguments.hasOption(BASE_FREQUENCIES)) {

                    option = arguments.getStringOption(BASE_FREQUENCIES);

                    if (option.equalsIgnoreCase(NUCLEOTIDE_FREQUENCIES)) {

                        int index = 0;
                        data.frequencyModelIndex = index;
                        if (arguments
                                .hasOption(NUCLEOTIDE_FREQUENCY_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(NUCLEOTIDE_FREQUENCY_PARAMETER_VALUES);
                            parseFrequencyValues(index, values, data);
                        }

                    } else if (option.equalsIgnoreCase(CODON_FREQUENCIES)) {

                        int index = 1;
                        data.frequencyModelIndex = index;
                        if (arguments
                                .hasOption(CODON_FREQUENCY_PARAMETER_VALUES)) {
                            values = arguments
                                    .getRealArrayOption(CODON_FREQUENCY_PARAMETER_VALUES);
                            parseFrequencyValues(index, values, data);
                        }

                    } else {
                        gracefullyExit("Unrecognized option.");
                    }

                }// END: FREQUENCY_MODEL option check

                if (arguments.hasOption(FROM)) {

                    data.from = arguments.getIntegerOption(FROM);

                }// END: FROM option check

                if (arguments.hasOption(TO)) {

                    data.to = arguments.getIntegerOption(TO);

                }// END: TO option check

                if (arguments.hasOption(EVERY)) {

                    data.every = arguments.getIntegerOption(EVERY);

                }// END: EVERY option check

                // create partition
                Partition partition = new Partition(
                        data.createTreeModel(), //
                        data.createBranchModel(), //
                        data.createSiteRateModel(), //
                        data.createClockRateModel(), //
                        data.createFrequencyModel(), //
                        data.from - 1, // from
                        data.to - 1, // to
                        data.every // every
                );

                if (arguments.hasOption(ROOT_SEQUENCE)) {

                    data.ancestralSequenceString = arguments.getStringOption(ROOT_SEQUENCE);
                    partition.setRootSequence(data.createAncestralSequence());

                }// END: ANCESTRAL_SEQUENCE option check

                partitionsList.add(partition);
                
//                System.err.println(data.from);
                
                dataList.add(data);

            }// END: partitionArgs loop

            // ////////////////
            // ---SIMULATE---//
            // ////////////////

            if (this.VERBOSE) {
            	
//        		System.out.println(dataList.get(0).from + " " + dataList.get(1).from);
            	
                Utils.printPartitionDataList(dataList);
                System.out.println();
            }

			SimpleAlignment alignment = new SimpleAlignment();
			String outputFile = null;
			if (leftoverArguments.length > 0) {

				outputFile = leftoverArguments[0];
				String[] file = outputFile.split("\\.", 2);
				if (file.length > 1) {

					String extension = outputFile.split("\\.", 2)[1];

					// TODO Delegate here to enum-class; switches are not generic
					if (extension.equalsIgnoreCase(SimpleAlignment.OutputType.FASTA
									.getText()) || extension.equalsIgnoreCase("fst")) {

						dataList.outputFormat = SimpleAlignment.OutputType.FASTA;

					} else if (extension.equalsIgnoreCase(SimpleAlignment.OutputType.NEXUS
									.getText()) || extension.equalsIgnoreCase("nxs")) {

						dataList.outputFormat = SimpleAlignment.OutputType.NEXUS;

					} else if (extension.equalsIgnoreCase(SimpleAlignment.OutputType.XML
									.getText())) {

						dataList.outputFormat = SimpleAlignment.OutputType.XML;

					} else {
						dataList.outputFormat = SimpleAlignment.OutputType.FASTA;
					}//END: extension check

				} else {
					outputFile = file[0] + "." + SimpleAlignment.OutputType.FASTA.toString().toLowerCase();
					dataList.outputFormat = SimpleAlignment.OutputType.FASTA;
				}//END: 

			}// END: leftoverArguments check

            if (leftoverArguments.length > 1) {
                dataList.startingSeed = Long.parseLong(leftoverArguments[1]);
                dataList.setSeed = true;
            }

            if (dataList.setSeed) {
                MathUtils.setSeed(dataList.startingSeed);
            }

            if (leftoverArguments.length > 2) {
                dataList.outputAncestralSequences = Boolean.parseBoolean(leftoverArguments[2]);
            }
            
            if (leftoverArguments.length > 3) {
                dataList.useParallel = Boolean.parseBoolean(leftoverArguments[3]);
            }

            BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
                    partitionsList);
            alignment = beagleSequenceSimulator.simulate(dataList.useParallel, dataList.outputAncestralSequences);
            alignment.setOutputType(dataList.outputFormat);

            PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
            writer.println(alignment.toString());
            writer.close();

		} catch (ArgumentException e) {

			System.out.println();
			printUsage(arguments);
			System.out.println();
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);

		} catch (IOException e) {
			
			System.out.println();
			printUsage(arguments);
			System.out.println();
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
			
		} catch (ImportException e) {
		
			System.out.println();
			printUsage(arguments);
			System.out.println();
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		
		}// END: try-catch block

    }// END: simulate

    private void parseDemographicValues(int demographicModelIndex,
                                        double[] values, 
                                        PartitionData data) {
        for (int i = 0; i < PartitionData.demographicParameterIndices[demographicModelIndex].length; i++) {

            int k = PartitionData.demographicParameterIndices[demographicModelIndex][i];
            data.demographicParameterValues[k] = values[i];

        }
    }// END: parseDemographicValues

    private void parseSubstitutionValues(int substitutionModelIndex,
                                         double[] values,
                                         PartitionData data) {
        for (int i = 0; i < PartitionData.substitutionParameterIndices[substitutionModelIndex].length; i++) {

            int k = PartitionData.substitutionParameterIndices[substitutionModelIndex][i];
            data.substitutionParameterValues[k] = values[i];

        }
    }// END: parseSubstitutionValues

    private void parseSiteRateValues(int siteRateModelIndex,
                                     double[] values,
                                     PartitionData data) {
        for (int i = 0; i < PartitionData.siteRateModelParameterIndices[siteRateModelIndex].length; i++) {

            int k = PartitionData.siteRateModelParameterIndices[siteRateModelIndex][i];
            data.siteRateModelParameterValues[k] = values[i];

        }
    }// END: parseSiteRateModelParameterValues

    private void parseClockValues(int clockModelIndex,
                                  double[] values,
                                  PartitionData data) {
        for (int i = 0; i < PartitionData.clockParameterIndices[clockModelIndex].length; i++) {

            int k = PartitionData.clockParameterIndices[clockModelIndex][i];
            data.clockParameterValues[k] = values[i];

        }
    }// END: parseClockValues

    private void parseFrequencyValues(int frequencyModelIndex,
                                      double[] values,
                                      PartitionData data) {
        for (int i = 0; i < data.frequencyParameterIndices[frequencyModelIndex].length; i++) {

            int k = data.frequencyParameterIndices[frequencyModelIndex][i];
            data.frequencyParameterValues[k] = values[i];

        }
    }// END: parseFrequencyValues

    private void gracefullyExit(String message) {
        printUsage(arguments);
        if (message != null) {
            System.out.println(message);
            System.out.println();
        }
        System.exit(0);
    }// END: gracefullyExit

    private void printUsage(Arguments arguments) {

        arguments.printUsage(
                "java -Djava.library.path=/usr/local/lib -jar buss.jar", " "
                + SPLIT_PARTITION + " " + "[<output-file-name>] [<seed>] [<true|false>] [<true|false>]");
        System.out.println();

        System.out
                .println("  Example: java -Djava.library.path=/usr/local/lib -jar pibuss.jar "
                        + "-treeFile SimTree.figtree -from 1 -to 500 -every 1 -branchSubstitutionModel HKY -HKYsubstitutionParameterValues 1.0"
                        + " "
                        + SPLIT_PARTITION
                        + " "
                        + "-treeFile SimTree.figtree -from 501 -to 1000 -every 1 -branchSubstitutionModel HKY -HKYsubstitutionParameterValues 10.0"
                        + " " + SPLIT_PARTITION + " " + "sequences.fasta");

        System.out.println();
    }// END: printUsage

}// END: class
