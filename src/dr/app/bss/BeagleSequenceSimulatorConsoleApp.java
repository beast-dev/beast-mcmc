package dr.app.bss;

import java.io.File;
import java.util.ArrayList;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.util.Arguments;
import dr.evomodel.tree.TreeModel;

public class BeagleSequenceSimulatorConsoleApp {

	private static final String HELP = "help";
	private static final String TREE_MODEL = "treeModel";
	private static final String BRANCH_SUBSTITUTION_MODEL = "branchSubstitutionModel";
	private static final String SITE_RATE_MODEL = "siteRateModel";
	private static final String CLOCK_RATE_MODEL = "clockRateModel";
	private static final String FREQUENCY_MODEL = "frequencyModel";
	private static final String FROM = "from";
	private static final String TO = "to";
	private static final String EVERY = "every";
	
	
	
	private String[] args;
	private Arguments arguments;
//   private PartitionDataList dataList;
	private PartitionData data;
	
	public BeagleSequenceSimulatorConsoleApp(String[] args) {

		this.args = args;
//		dataList = new PartitionDataList();
//		dataList.add(new PartitionData());
		
		 data = new PartitionData();
		
		// //////////////////
		// ---DEFINITION---//
		// //////////////////

		// TODO: set model parameter values
		arguments = new Arguments(new Arguments.Option[] {

				new Arguments.Option(HELP, "print this information and exit"),

				new Arguments.StringOption(TREE_MODEL, "tree model",
						"specify tree topology"),

				new Arguments.StringOption(BRANCH_SUBSTITUTION_MODEL,
						new String[] { PartitionData.substitutionModels[0], // HKY
								PartitionData.substitutionModels[1], // GTR
								PartitionData.substitutionModels[2], // TN93
								PartitionData.substitutionModels[3] // GY94CodonModel
						}, false, "specify substitution model"),

				new Arguments.StringOption(SITE_RATE_MODEL, new String[] {
						"NoModel", // NoModel
						"GammaSiteRateModel", // GammaSiteRateModel
				}, false, "specify site rate model"),

				new Arguments.StringOption(CLOCK_RATE_MODEL, new String[] {
						"StrictClock", // StrictClock
						PartitionData.clockModels[1], // LognormalRelaxedClock
						PartitionData.clockModels[2] // ExponentialRelaxedClock
						}, false, "specify clock rate model"),
				
				new Arguments.StringOption(FREQUENCY_MODEL, new String[] {
						"NucleotideFrequencies", // NucleotideFrequencies
						PartitionData.frequencyModels[1], // CodonFrequencies
						}, false, "specify frequency model"),
				
				new Arguments.IntegerOption(FROM, "specify 'from' attribute"),
				new Arguments.IntegerOption(TO, "specify 'to' attribute"),
				new Arguments.IntegerOption(EVERY, "specify 'every' attribute")
		});

	}// END: constructor

// -tree /home/filip/SimTree.figtree -branchSubstitutionModel HKY -siteRateModel GammaSiteRateModel -clockRateModel StrictClock -frequencyModel NucleotideFrequencies
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

			if (args.length == 0 || arguments.hasOption(HELP)) {
				printUsage(arguments);
				System.exit(0);
			}


			// Tree Model
			if (arguments.hasOption(TREE_MODEL)) {

				data.treeFile = new File(arguments.getStringOption(TREE_MODEL));

			} else {
				
				throw new RuntimeException("TreeModel not specified.");
				
			}
			
			// Branch Substitution Model
	        if (arguments.hasOption(BRANCH_SUBSTITUTION_MODEL)) {
	        	
//	        	System.out.println(arguments.getStringOption(BRANCH_SUBSTITUTION_MODEL));
	        	
	        }
			
	        
	        if (arguments.hasOption(FROM)) {
	        	
	        	data.from = arguments.getIntegerOption(FROM);
	        	
	        }
	        
			
	        if (arguments.hasOption(TO)) {
	        	
	        	data.to = arguments.getIntegerOption(TO);
	        	
	        }
	        
	        if (arguments.hasOption(EVERY)) {
	        	
	        	data.every = arguments.getIntegerOption(EVERY);
	        	
	        }
	        
			// ////////////////
			// ---SIMULATE---//
			// ////////////////
	        //TODO: loop over partitions here

	        Utils.printPartitionData(data);
	        
	        ArrayList<Partition> partitionsList = new ArrayList<Partition>();
	        
			// create partition
			Partition partition = new Partition(
					data.createTreeModel(), //
					data.createBranchModel(), //
					data.createSiteRateModel(), //
					data.createBranchRateModel(), //
					data.createFrequencyModel(), //
					data.from - 1, // from
					data.to - 1, // to
					data.every // every
			);
	        
			partitionsList.add(partition);
			
			BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
					partitionsList
					);
			
			beagleSequenceSimulator.simulate().toString();
			
		} catch (Arguments.ArgumentException ae) {
			System.out.println();
			System.out.println(ae.getMessage());
			System.out.println();
			printUsage(arguments);
			System.exit(1);
		}

	}

	private void printUsage(Arguments arguments) {

		arguments.printUsage("java -jar bss.jar", "[<output-file-name>]");
		System.out.println();
		System.out.println("  Example: java -jar bss.jar -help");
		System.out.println();
	}

}// END: class
