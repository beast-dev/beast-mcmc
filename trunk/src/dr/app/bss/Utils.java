package dr.app.bss;

import java.awt.Desktop;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class Utils {

	// /////////////////
	// ---CONSTANTS---//
	// /////////////////

	public static final boolean VERBOSE = false;

	// public static final int TREE_MODEL_ELEMENT = 0;
	public static final int BRANCH_MODEL_ELEMENT = 1;
	public static final int SITE_RATE_MODEL_ELEMENT = 2;
	public static final int BRANCH_RATE_MODEL_ELEMENT = 3;
	public static final int FREQUENCY_MODEL_ELEMENT = 4;
	public static final int DEMOGRAPHIC_MODEL_ELEMENT = 5;
	
//	public static final String TAXA = "taxa";
	public static final String TOPOLOGY = "topology";
	public static final String ABSOLUTE_HEIGHT = "absoluteHeight";
	public static final String TREE_FILENAME = "treeFilename";
	public static final String SUBSTITUTION_MODEL = "substitutionModel";
	public static final String DEMOGRAPHIC_MODEL = "demographicModel";
	public static final String FREQUENCY_MODEL = "frequencyModel";
	public static final String CODON_UNIVERSAL = "codon-universal";
	public static final String CHOOSE_FILE = "Choose file...";
	public static final String EDIT_TAXA_SET = "Edit taxa set...";
	
	public static final String BSS_ICON = "icons/bss.png";
	public static final String CHECK_ICON = "icons/check.png";
	public static final String ERROR_ICON = "icons/error.png";
	public static final String HAMMER_ICON = "icons/hammer.png";
	public static final String CLOSE_ICON = "icons/close.png";
	public static final String BIOHAZARD_ICON = "icons/biohazard.png";
	public static final String BUBBLE_BLUE_ICON = "icons/bubble-blue.png";
	public static final String SAVE_ICON = "icons/save.png";
	public static final String TEXT_FILE_ICON = "icons/file.png";
	
	// ///////////////////////////////
	// ---GENERAL UTILITY METHODS---//
	// ///////////////////////////////

	public static Tree importTreeFromFile(File file) {

		Tree tree = null;
		
		try {

			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = reader.readLine();

			if (line.toUpperCase().startsWith("#NEXUS")) {

				NexusImporter importer = new NexusImporter(reader);
				tree = importer.importTree(null);

			} else {

				NewickImporter importer = new NewickImporter(reader);
				tree = importer.importTree(null);

			}

			reader.close();

		} catch (Exception e) {
			Utils.handleException(e);
		}// END: try-catch block

		return tree;
	}// END: importTreeFromFile
	
	public static void removeTaxaWithAttributeValue(PartitionDataList dataList,
			String attribute, String value) {

		synchronized (dataList.allTaxa) {
			for (int i = 0; i < dataList.allTaxa.getTaxonCount(); i++) {

				Taxon taxon = dataList.allTaxa.getTaxon(i);
				if (taxon.getAttribute(attribute).toString()
						.equalsIgnoreCase(value)) {
					dataList.allTaxa.removeTaxon(taxon);
					i--;
				}
			}
		}

	}// END: removeTaxaWithAttributeValue
	
	public static void centreLine(String line, int pageWidth) {
		int n = pageWidth - line.length();
		int n1 = n / 2;
		for (int i = 0; i < n1; i++) {
			System.out.print(" ");
		}
		System.out.println(line);
	}

	public static int getSiteCount(PartitionDataList dataList) {

		int siteCount = 0;
		int to = 0;
		for (PartitionData data : dataList) {
			// siteCount += data.createPartitionSiteCount();
			to = data.to;
			if (to > siteCount) {
				siteCount = to;
			}

		}

		return siteCount;// + 1;
	}// END: getSiteCount

	public static int arrayIndex(String[] array, String element) {

		List<String> vector = new ArrayList<String>();
		for (int i = 0; i < array.length; i++) {
			vector.add(array[i]);
		}

		return vector.indexOf(element);
	}// END: arrayIndex

	public static ArrayList<TreeModel> treesToList(PartitionDataList dataList) {

		ArrayList<TreeModel> treeModelsList = new ArrayList<TreeModel>();
		for (PartitionData data : dataList) {
			treeModelsList.add(data.createTreeModel());
		}

		return treeModelsList;
	}// END: treesToList

	public static boolean taxonExists(Taxon taxon, MutableTaxonList taxonList) {

		boolean exists = false;
		for (Taxon taxon2 : taxonList) {

			if (taxon.equals(taxon2) 
//					&& taxon.getAttribute(Utils.TREE_FILENAME).toString().equalsIgnoreCase(taxon2.getAttribute(Utils.TREE_FILENAME).toString())
					) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: taxonExists

//	 private boolean isFileInList(File file) {
//	 boolean exists = false;
//	
//	 for (File file2 : dataList.treesList) {
//	
//	 if (file.getName().equalsIgnoreCase(file2.getName())) {
//	 exists = true;
//	 break;
//	 }
//	
//	 }
//	
//	 return exists;
//	 }// END: isFileInList
	
	public static double getAbsoluteTaxonHeight(Taxon taxon, Tree tree) {

		double height = 0.0;
		for (int i = 0; i < tree.getExternalNodeCount(); i++) {

			NodeRef externalNode = tree.getExternalNode(i);
			Taxon externalNodeTaxon = tree.getNodeTaxon(externalNode);

			if (externalNodeTaxon.equals(taxon)) {
				height = tree.getNodeHeight(externalNode);
			}
		}// END: external node loop

		return height;
	}// END: getAbsoluteTaxonHeight

	public static boolean isRecordInList(TreesTableRecord record,
			ArrayList<TreesTableRecord> recordsList) {

		boolean exists = false;

		for (TreesTableRecord record2 : recordsList) {

			if (record.getName().equalsIgnoreCase(record2.getName())) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: isRecordInList
	
	//TODO: check if works
	public static boolean isTaxaInList(Taxa taxa,
			ArrayList<Taxa> taxaList) {

		boolean exists = false;

		for (Taxa taxa2 : taxaList) {

			if (taxaToString(taxa).equalsIgnoreCase(taxaToString(taxa2))) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: isTaxaInList
	
	public static int taxaIsIdenticalWith(Taxa taxa,
			ArrayList<Taxa> taxaList) {

		int index = -Integer.MAX_VALUE;

		for (Taxa taxa2 : taxaList) {

			if (taxaToString(taxa).equalsIgnoreCase(taxaToString(taxa2))) {
				index = taxaList.indexOf(taxa2);
				break;
			}

		}

		return index;
	}// END: treeModelIsIdenticalWith
	
	public static boolean isTreeModelInList(TreeModel treeModel,
			ArrayList<TreeModel> treeModelList) {

		boolean exists = false;

		for (TreeModel treeModel2 : treeModelList) {

			if (treeModel.getNewick().equalsIgnoreCase(treeModel2.getNewick())) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: isTreeModelInList

	public static int treeModelIsIdenticalWith(TreeModel treeModel,
			ArrayList<TreeModel> treeModelList) {

		int index = -Integer.MAX_VALUE;

		for (TreeModel treeModel2 : treeModelList) {

			if (treeModel.getNewick().equalsIgnoreCase(treeModel2.getNewick())) {
				index = treeModelList.indexOf(treeModel2);
				break;
			}

		}

		return index;
	}// END: treeModelIsIdenticalWith

	public static boolean isElementInList(PartitionData data,
			ArrayList<PartitionData> partitionList, int elementIndex) {

		boolean exists = false;

		switch (elementIndex) {

		case DEMOGRAPHIC_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (demographicModelToString(data).equalsIgnoreCase(
						demographicModelToString(data2))) {
					exists = true;
					break;
				}
			}

			break;
		
		case BRANCH_RATE_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (clockRateModelToString(data).equalsIgnoreCase(
						clockRateModelToString(data2))) {
					exists = true;
					break;
				}
			}

			break;

		case FREQUENCY_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (frequencyModelToString(data).equalsIgnoreCase(
						frequencyModelToString(data2))) {
					exists = true;
					break;
				}
			}

			break;

		case BRANCH_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (branchSubstitutionModelToString(data).equalsIgnoreCase(
						branchSubstitutionModelToString(data2))) {
					exists = true;
					break;
				}
			}

			break;

		case SITE_RATE_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (siteRateModelToString(data).equalsIgnoreCase(
						siteRateModelToString(data2))) {
					exists = true;
					break;
				}
			}

			break;

		default:

			throw new RuntimeException("Unknown element");

		}// END: switch

		return exists;
	}// END: isModelInList

	public static int isIdenticalWith(PartitionData data,
			ArrayList<PartitionData> partitionList, int elementIndex) {

		int index = -Integer.MAX_VALUE;

		switch (elementIndex) {

		case DEMOGRAPHIC_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (demographicModelToString(data).equalsIgnoreCase(
						demographicModelToString(data2))) {
					index = partitionList.indexOf(data2);
					break;
				}
			}

			break;
		
		case BRANCH_RATE_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (clockRateModelToString(data).equalsIgnoreCase(
						clockRateModelToString(data2))) {
					index = partitionList.indexOf(data2);
					break;
				}
			}

			break;

		case FREQUENCY_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (frequencyModelToString(data).equalsIgnoreCase(
						frequencyModelToString(data2))) {
					index = partitionList.indexOf(data2);
					break;
				}
			}

			break;

		case BRANCH_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (branchSubstitutionModelToString(data).equalsIgnoreCase(
						branchSubstitutionModelToString(data2))) {
					index = partitionList.indexOf(data2);
					break;
				}
			}

			break;

		case SITE_RATE_MODEL_ELEMENT:

			for (PartitionData data2 : partitionList) {
				if (siteRateModelToString(data).equalsIgnoreCase(
						siteRateModelToString(data2))) {
					index = partitionList.indexOf(data2);
					break;
				}
			}

			break;

		default:

			throw new RuntimeException("Unknown element");

		}// END: switch

		return index;
	}// END: isIdenticalWith


	// /////////////////
	// ---GUI UTILS---//
	// /////////////////

	public static ImageIcon createImageIcon(String path) {

		ImageIcon icon = null;

		try {

			URL imgURL = BeagleSequenceSimulatorApp.class.getResource(path);

			if (imgURL != null) {
				icon = new ImageIcon(imgURL);
			} else {
				System.err.println("Couldn't find file: " + path + "\n");
			}

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

		return icon;
	}// END: CreateImageIcon

	public static boolean isBrowsingSupported() {
		if (!Desktop.isDesktopSupported()) {
			return false;
		}
		boolean result = false;
		Desktop desktop = java.awt.Desktop.getDesktop();
		if (desktop.isSupported(Desktop.Action.BROWSE)) {
			result = true;
		}
		return result;

	}// END: isBrowsingSupported

	public static int getTabbedPaneComponentIndex(JTabbedPane tabbedPane,
			String title) {

		int index = -Integer.MAX_VALUE;

		int count = tabbedPane.getTabCount();
		for (int i = 0; i < count; i++) {
			if (tabbedPane.getTitleAt(i).toString().equalsIgnoreCase(title)) {
				index = i;
				break;
			}// END: title check

		}// END: i loop

		return index;
	}// END: getComponentIndex

	public static Frame getActiveFrame() {
		Frame result = null;
		Frame[] frames = Frame.getFrames();
		for (int i = 0; i < frames.length; i++) {
			Frame frame = frames[i];
			if (frame.isVisible()) {
				result = frame;
				break;
			}
		}
		return result;
	}

	public static String getMultipleWritePath(File outFile,
			String defaultExtension, int i) {

		String path = outFile.getParent();
		String[] nameArray = outFile.getName().split("\\.", 2);
		String name = ((i == 0) ? nameArray[0] : nameArray[0] + i);

		String extension = (nameArray.length == 1) ? (defaultExtension)
				: (nameArray[1]);
		String fullPath = path + System.getProperty("file.separator") + name
				+ "." + extension;

		return fullPath;
	}// END: getMultipleWritePath

	public static String getWritePath(File outFile, String defaultExtension) {

		String path = outFile.getParent();
		String[] nameArray = outFile.getName().split("\\.", 2);
		String name = nameArray[0];

		String extension = (nameArray.length == 1) ? (defaultExtension)
				: (nameArray[1]);
		String fullPath = path + System.getProperty("file.separator") + name
				+ "." + extension;

		return fullPath;
	}// END: getWritePath

	public static void showDialog(final String message) {

		if (SwingUtilities.isEventDispatchThread()) {

			JOptionPane.showMessageDialog(getActiveFrame(), message, "Message",
					JOptionPane.ERROR_MESSAGE,
					Utils.createImageIcon(Utils.BUBBLE_BLUE_ICON));

		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					JOptionPane.showMessageDialog(getActiveFrame(), message,
							"Message", JOptionPane.ERROR_MESSAGE,
							Utils.createImageIcon(Utils.BUBBLE_BLUE_ICON));

				}
			});
		}// END: edt check
	}// END: showDialog

	// ////////////////////////////////
	// ---EXCEPTION HANDLING UTILS---//
	// ////////////////////////////////

	public static void handleException(final Throwable e, final String message) {

		final Thread t = Thread.currentThread();

		if (SwingUtilities.isEventDispatchThread()) {
			showExceptionDialog(t, e, message);
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showExceptionDialog(t, e, message);
				}
			});
		}// END: edt check
	}// END: uncaughtException

	public static void handleException(final Throwable e) {

		final Thread t = Thread.currentThread();

		if (SwingUtilities.isEventDispatchThread()) {
			showExceptionDialog(t, e);
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showExceptionDialog(t, e);
				}
			});
		}// END: edt check
	}// END: handleException

	private static void showExceptionDialog(Thread t, Throwable e) {

		String msg = String.format("Unexpected problem on thread %s: %s",
				t.getName(), e.getMessage());

		logException(t, e);

		JOptionPane.showMessageDialog(Utils.getActiveFrame(), //
				msg, //
				"Error", //
				JOptionPane.ERROR_MESSAGE, //
				Utils.createImageIcon(Utils.ERROR_ICON));
	}// END: showExceptionDialog

	private static void showExceptionDialog(Thread t, Throwable e,
			String message) {

		String msg = String.format("Unexpected problem on thread %s: %s" + "\n"
				+ message, t.getName(), e.getMessage());

		logException(t, e);

		JOptionPane.showMessageDialog(Utils.getActiveFrame(), //
				msg, //
				"Error", //
				JOptionPane.ERROR_MESSAGE, //
				Utils.createImageIcon(Utils.ERROR_ICON));
	}// END: showExceptionDialog

	private static void logException(Thread t, Throwable e) {
		e.printStackTrace();
	}// END: logException

	// ///////////////////
	// ---PRINT UTILS---//
	// ///////////////////

	public static void printMap(Map<?, ?> mp) {
		Iterator<?> it = mp.entrySet().iterator();
		while (it.hasNext()) {
			Entry<?, ?> pairs = (Entry<?, ?>) it.next();
			Object obj = pairs.getValue();
			if (obj instanceof int[]) {
				int[] seq = (int[]) obj;
				System.out.print(pairs.getKey() + " =");
				for (int i = 0; i < seq.length; ++i) {
					System.out.print(" " + seq[i]);
				}
				System.out.println();
			} else {
				System.out.println(pairs.getKey() + " = " + pairs.getValue());
			}
		}
	}// END: printMap
	
	public static void printHashMap(ConcurrentHashMap<?, ?> hashMap) {

		Iterator<?> iterator = hashMap.entrySet().iterator();
		while (iterator.hasNext()) {

			Entry<?, ?> pairs = (Entry<?, ?>) iterator.next();

			Taxon taxon = (Taxon) pairs.getKey();
			int[] sequence = (int[]) pairs.getValue();

			System.out.println(taxon.toString());
			Utils.printArray(sequence);

		}// END: while has next

	}// END: printHashMap

	public static void printArray(Object[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.println(x[i]);
		}
	}// END: printArray

	public static void printArray(int[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.print(x[i] + " ");
		}
		System.out.println();
	}// END: printArray

	public static void printArray(double[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.println(x[i]);
		}
	}// END: printArray
	
	public static void printArray(String[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.println(x[i]);
		}
	}// END: printArray

	public static void print2DArray(double[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(array[row][col] + " ");
			}
			System.out.print("\n");
		}
	}// END: print2DArray

	public static void printBranchSubstitutionModel(PartitionData data) {
		System.out.print("\tBranch Substitution model: ");
		System.out.print(branchSubstitutionModelToString(data));
		System.out.print("\n");
	}// END: printBranchSubstitutionModel

	public static void printClockRateModel(PartitionData data) {
		System.out.print("\tClock rate model: ");
		System.out.print(clockRateModelToString(data));
		System.out.print("\n");
	}// END: printClockRateModel

	public static void printFrequencyModel(PartitionData data) {
		System.out.print("\tFrequency model: ");
		System.out.print(frequencyModelToString(data));
		System.out.print("\n");
	}// END: printFrequencyModel

	public static void printSiteRateModel(PartitionData data) {
		System.out.print("\tSite rate model: ");
		System.out.print(siteRateModelToString(data));
		System.out.print("\n");
	}// END: printFrequencyModel

	public static void printPartitionData(PartitionData data) {

		// System.out.println("\tData type: "+
		// PartitionData.dataTypes[data.dataTypeIndex]);
//		System.out.println("\tTree model: " + data.treeFile);
		System.out.println("\tTree model: " + data.createTreeModel().toString());
		System.out.println("\tFrom: " + data.from);
		System.out.println("\tTo: " + data.to);
		System.out.println("\tEvery: " + data.every);
		printBranchSubstitutionModel(data);
		printSiteRateModel(data);
		printClockRateModel(data);
		printFrequencyModel(data);

	}// END: printPartitionData

	public static void printPartitionDataList(PartitionDataList dataList) {

		// printTaxonList(dataList);
		System.out.println("\tSite count: " + getSiteCount(dataList));
		if (dataList.setSeed) {
			System.out.println("\tStarting seed: " + dataList.startingSeed);
		}

		int row = 1;
		for (PartitionData data : dataList) {

			System.out.println("Partition: " + row);
			printPartitionData(data);
			row++;

		}// END: data list loop

	}// END: printDataList

//	public static void printTreeFileList(PartitionDataList dataList) {
//
//		int i = 0;
//		for (File treeFile : dataList.treesList) {
//			System.out.println(treeFile + " (" + dataList.taxaCounts.get(i)
//					+ ")");
//			i++;
//		}
//
//	}// END: printForestList

	public static void printTaxonList(PartitionDataList dataList) {
		for (int i = 0; i < dataList.allTaxa.getTaxonCount(); i++) {

			System.out.println(dataList.allTaxa.getTaxon(i).getId()
					+ ": "
					+ dataList.allTaxa.getTaxon(i).getAttribute(
							Utils.ABSOLUTE_HEIGHT));

		}// END: taxon loop
	}// END: printTaxonList

	// //////////////////////
	// ---TOSTRING UTILS---//
	// //////////////////////

	public static String taxonToString(Taxon taxon) {
		String string = taxon.getId() + " ("
				+ taxon.getAttribute(Utils.ABSOLUTE_HEIGHT) + ","
				+ taxon.getAttribute(Utils.TREE_FILENAME) + ")";
		return string;
	}// END: taxonToString
	
	public static String taxaToString(Taxa taxa) {

		String string = "";

		for (int i = 0; i < taxa.getTaxonCount(); i++) {
			
			Taxon taxon = taxa.getTaxon(i);
			string += taxonToString(taxon) + ("\n");
			
		}
		
		return string;
	}// END: taxaToString
	
	public static String partitionDataToString(PartitionData data) {

		String string = "";

		//TODO: for coalescent this simulates again, so that's bad
		string += ("Tree model: " + data.createTreeModel().toString())+ ("\n");
		string += ("From: " + data.from)+ ("\n");
		string += ("To: " + data.to)+ ("\n");
		string += ("Every: " + data.every)+ ("\n");

		string += ("Branch Substitution model: ") + branchSubstitutionModelToString(data) + ("\n");
		string += ("Site Rate model: ") + siteRateModelToString(data) + ("\n");
		string += ("Clock Rate model: ") + clockRateModelToString(data) + ("\n");
		string += ("Frequency model: ") + frequencyModelToString(data) + ("\n");

		return string;
	}// END: partitionDataToString
	
	public static String partitionDataListToString(PartitionDataList dataList) {

		String string = "";

		string += ("Site count: " + getSiteCount(dataList)) + ("\n");
		if (dataList.setSeed) {
			string += ("Starting seed: " + dataList.startingSeed) + ("\n");
		}

		int row = 1;
		for (PartitionData data : dataList) {

			string += ("Partition: " + row)+ ("\n");
			string += partitionDataToString(data);
			string += ("\n");
			row++;

		}// END: data list loop

		return string;
	}// END: partitionDataListToString
	
	public static String demographicModelToString(PartitionData data) {

		String string = PartitionData.demographicModels[data.demographicModelIndex];

		string += (" ( ");
		for (int i = 0; i < PartitionData.demographicParameterIndices[data.demographicModelIndex].length; i++) {
			string += data.demographicParameterValues[PartitionData.demographicParameterIndices[data.demographicModelIndex][i]];
			string += " ";
		}// END: indices loop
		string += ")";

		return string;
	}
	
	public static String clockRateModelToString(PartitionData data) {

		String string = PartitionData.clockModels[data.clockModelIndex];

		string += (" ( ");
		for (int i = 0; i < PartitionData.clockParameterIndices[data.clockModelIndex].length; i++) {
			string += data.clockParameterValues[PartitionData.clockParameterIndices[data.clockModelIndex][i]];
			string += " ";
		}// END: indices loop
		string += ")";

		return string;
	}

	public static String frequencyModelToString(PartitionData data) {

		String string = PartitionData.frequencyModels[data.frequencyModelIndex];

		string += (" ( ");
		for (int i = 0; i < data.frequencyParameterIndices[data.frequencyModelIndex].length; i++) {
			string += data.frequencyParameterValues[data.frequencyParameterIndices[data.frequencyModelIndex][i]];
			string += " ";
		}// END: indices loop
		string += ")";

		return string;
	}

	public static String branchSubstitutionModelToString(PartitionData data) {

		String string = PartitionData.substitutionModels[data.substitutionModelIndex];

		string += (" ( ");
		for (int i = 0; i < PartitionData.substitutionParameterIndices[data.substitutionModelIndex].length; i++) {
			string += data.substitutionParameterValues[PartitionData.substitutionParameterIndices[data.substitutionModelIndex][i]];
			string += " ";
		}// END: indices loop
		string += ")";

		return string;
	}

	public static String siteRateModelToString(PartitionData data) {

		String string = PartitionData.siteRateModels[data.siteRateModelIndex];

		string += (" ( ");
		for (int i = 0; i < PartitionData.siteRateModelParameterIndices[data.siteRateModelIndex].length; i++) {
			string += data.siteRateModelParameterValues[PartitionData.siteRateModelParameterIndices[data.siteRateModelIndex][i]];
			string += " ";
		}// END: indices loop
		string += ")";

		return string;
	}
	
}// END: class
