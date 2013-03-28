package dr.app.bss;

import java.io.File;
import java.util.ArrayList;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.util.Arguments;

public class BeagleSequenceSimulatorConsoleApp {

	private String[] args;
	private Arguments arguments;
	// private PartitionDataList dataList;
	private PartitionData data;

	private static final String HELP = "help";

	private static final String TREE_MODEL = "treeModel";

	private static final String BRANCH_SUBSTITUTION_MODEL = "branchSubstitutionModel";
	private static final String HKY = PartitionData.substitutionModels[0];
	private static final String GTR = PartitionData.substitutionModels[1];
	private static final String TN93 = PartitionData.substitutionModels[2];
	private static final String GY94_CODON_MODEL = PartitionData.substitutionModels[3];

	private static final String SITE_RATE_MODEL = "siteRateModel";
	private static final String NO_MODEL = "NoModel";
	private static final String GAMMA_SITE_RATE_MODEL = "gammaSiteRateModel";

	private static final String CLOCK_RATE_MODEL = "clockRateModel";
	private static final String STRICT_CLOCK = "strictClock";
	private static final String LOGNORMAL_RELAXED_CLOCK = "lognormalRelaxedClock";
	private static final String EXPONENTIAL_RELAXED_CLOCK = "exponentialRelaxedClock";

	private static final String FREQUENCY_MODEL = "frequencyModel";
	private static final String NUCLEOTIDE_FREQUENCIES = "nucleotideFrequencies";
	private static final String CODON_FREQUENCIES = "codonFrequencies";

	private static final String FROM = "from";
	private static final String TO = "to";
	private static final String EVERY = "every";

	public BeagleSequenceSimulatorConsoleApp(String[] args) {

		this.args = args;
		// dataList = new PartitionDataList();
		// dataList.add(new PartitionData());

		data = new PartitionData();

		// //////////////////
		// ---DEFINITION---//
		// //////////////////

		// TODO: verbose for printing state of data
		// TODO: set model parameter values
		arguments = new Arguments(
				new Arguments.Option[] {

						new Arguments.Option(HELP,
								"print this information and exit"),

						new Arguments.StringOption(TREE_MODEL, "tree model",
								"specify tree topology"),

						new Arguments.StringOption(BRANCH_SUBSTITUTION_MODEL,
								new String[] { HKY, //
										GTR, //
										TN93, //
										GY94_CODON_MODEL //
								}, false, "specify substitution model"),

						new Arguments.StringOption(SITE_RATE_MODEL,
								new String[] { NO_MODEL, //
										GAMMA_SITE_RATE_MODEL, //
								}, false, "specify site rate model"),

						new Arguments.StringOption(CLOCK_RATE_MODEL,
								new String[] { STRICT_CLOCK, //
										LOGNORMAL_RELAXED_CLOCK, //
										EXPONENTIAL_RELAXED_CLOCK //
								}, false, "specify clock rate model"),

						new Arguments.StringOption(FREQUENCY_MODEL,
								new String[] { NUCLEOTIDE_FREQUENCIES, //
										CODON_FREQUENCIES, //
								}, false, "specify frequency model"),

						new Arguments.IntegerOption(FROM,
								"specify 'from' attribute"),
						new Arguments.IntegerOption(TO,
								"specify 'to' attribute"),
						new Arguments.IntegerOption(EVERY,
								"specify 'every' attribute") });

	}// END: constructor

	// -tree /home/filip/SimTree.figtree -branchSubstitutionModel HKY
	// -siteRateModel GammaSiteRateModel -clockRateModel StrictClock
	// -frequencyModel NucleotideFrequencies
	// -from 1 -to 10 -every 1 sequences.fasta
	public void simulate() {

		try {

			// ///////////////
			// ---PARSING---//
			// ///////////////

			arguments.parseArguments(args);

			// ///////////////////
			// ---INTERROGATE---//
			// ///////////////////
			String option = null;

			if (args.length == 0 || arguments.hasOption(HELP)) {
				printUsage(arguments);
				System.exit(0);
			}// END: HELP option check

			// Tree Model
			if (arguments.hasOption(TREE_MODEL)) {

				data.treeFile = new File(arguments.getStringOption(TREE_MODEL));

			} else {

				throw new RuntimeException("TreeModel not specified.");

			}// END: TREE_MODEL option check

			// Branch Substitution Model
			if (arguments.hasOption(BRANCH_SUBSTITUTION_MODEL)) {

				option = arguments.getStringOption(BRANCH_SUBSTITUTION_MODEL);

				if (option.equalsIgnoreCase(HKY)) {
					data.substitutionModelIndex = 0;
				} else if (option.equalsIgnoreCase(GTR)) {
					data.substitutionModelIndex = 1;
				} else if (option.equalsIgnoreCase(TN93)) {
					data.substitutionModelIndex = 2;
				} else if (option.equalsIgnoreCase(GY94_CODON_MODEL)) {
					data.substitutionModelIndex = 3;
				} else {
					gracefullExit("Unrecognized option.");
				}

			}// END: BRANCH_SUBSTITUTION_MODEL option check

			// Site Rate Model
			if (arguments.hasOption(SITE_RATE_MODEL)) {

				option = arguments.getStringOption(SITE_RATE_MODEL);

				if (option.equalsIgnoreCase(NO_MODEL)) {
					data.siteRateModelIndex = 0;
				} else if (option.equalsIgnoreCase(GAMMA_SITE_RATE_MODEL)) {
					data.siteRateModelIndex = 1;
				} else {
					gracefullExit("Unrecognized option.");
				}

			}// END: SITE_RATE_MODEL option check

			// Clock Rate Model
			if (arguments.hasOption(CLOCK_RATE_MODEL)) {

				option = arguments.getStringOption(CLOCK_RATE_MODEL);

				if (option.equalsIgnoreCase(STRICT_CLOCK)) {
					data.clockModelIndex = 0;
				} else if (option.equalsIgnoreCase(LOGNORMAL_RELAXED_CLOCK)) {
					data.clockModelIndex = 1;
				} else if (option.equalsIgnoreCase(EXPONENTIAL_RELAXED_CLOCK)) {
					data.clockModelIndex = 2;
				} else {
					gracefullExit("Unrecognized option.");
				}

			}// END: CLOCK_RATE_MODEL option check

			// Frequency Model
			if (arguments.hasOption(FREQUENCY_MODEL)) {

				option = arguments.getStringOption(FREQUENCY_MODEL);

				if (option.equalsIgnoreCase(NUCLEOTIDE_FREQUENCIES)) {
					data.clockModelIndex = 0;
				} else if (option.equalsIgnoreCase(CODON_FREQUENCIES)) {
					data.clockModelIndex = 1;
				} else {
					gracefullExit("Unrecognized option.");
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

			
			String[] leftoverArguments = arguments.getLeftoverArguments();
			if(leftoverArguments.length > 1) {
				gracefullExit("Unrecognized option " + leftoverArguments[1]);
			}
			
			String outputFile = null;
			if(leftoverArguments.length > 0) {
				
				outputFile = leftoverArguments[0];
				
			} else {
				
				outputFile = "sequences.fasta";
				
			}
			
			
			
			
			// ////////////////
			// ---SIMULATE---//
			// ////////////////
			// TODO: loop over partitions here

			Utils.printPartitionData(data);

			ArrayList<Partition> partitionsList = new ArrayList<Partition>();

			// create partition
			Partition partition = new Partition(data.createTreeModel(), //
					data.createBranchModel(), //
					data.createSiteRateModel(), //
					data.createClockRateModel(), //
					data.createFrequencyModel(), //
					data.from - 1, // from
					data.to - 1, // to
					data.every // every
			);

			partitionsList.add(partition);

			BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
					partitionsList);

			beagleSequenceSimulator.simulate().toString();

			//TODO: write to this file
			System.out.println(outputFile);
			
		} catch (Arguments.ArgumentException ae) {
			System.out.println();
			System.out.println(ae.getMessage());
			System.out.println();
			printUsage(arguments);
			System.exit(1);
		}

	}

	private void gracefullExit(String message) {
		System.err.println(message);
		printUsage(arguments);
		System.exit(0);
	}//END: gracefullExit
	
	private void printUsage(Arguments arguments) {

		arguments.printUsage("java -jar bss.jar", "[<output-file-name>]");
		System.out.println();
		System.out.println("  Example: java -jar bss.jar -help");
		System.out.println();
	}

}// END: class
