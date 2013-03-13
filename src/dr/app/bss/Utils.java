package dr.app.bss;

import java.awt.Frame;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;

public class Utils {

	// /////////////////
	// ---CONSTANTS---//
	// /////////////////

	// public static final int TREE_MODEL_ELEMENT = 0;
	public static final int BRANCH_MODEL_ELEMENT = 1;
	public static final int SITE_RATE_MODEL_ELEMENT = 2;
	public static final int BRANCH_RATE_MODEL_ELEMENT = 3;
	public static final int FREQUENCY_MODEL_ELEMENT = 4;
	public static final String ABSOLUTE_HEIGHT = "absoluteHeight";

	// ///////////////////////////////
	// ---GENERAL UTILITY METHODS---//
	// ///////////////////////////////

	public static void centreLine(String line, int pageWidth) {
		int n = pageWidth - line.length();
		int n1 = n / 2;
		for (int i = 0; i < n1; i++) {
			System.out.print(" ");
		}
		System.out.println(line);
	}

	public static void printMap(Map<?, ?> mp) {
		Iterator<?> it = mp.entrySet().iterator();
		while (it.hasNext()) {
			Entry<?, ?> pairs = (Entry<?, ?>) it.next();
			System.out.println(pairs.getKey() + " = " + pairs.getValue());
			// it.remove(); // avoids a ConcurrentModificationException
		}
	}// END: printMap

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
			treeModelsList.add(data.treeModel);
		}

		return treeModelsList;
	}// END: treesToList

	public static boolean taxonExists(Taxon taxon, MutableTaxonList taxonList) {

		boolean exists = false;
		for (Taxon taxon2 : taxonList) {

			if (taxon.equals(taxon2)) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: taxonExists

	public static double getAbsoluteTaxonHeight(Taxon taxon, TreeModel tree) {

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

	public static boolean isTreeModelInList(TreeModel treeModel,
			ArrayList<TreeModel> treeModelList) {

		boolean exists = false;

		for (TreeModel treeModel2 : treeModelList) {

			if (treeModel.equals(treeModel2)) {
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

			if (treeModel.equals(treeModel2)) {
				index = treeModelList.indexOf(treeModel2);
				break;
			}

		}

		return index;
	}// END: treeModelIsIdenticalWith

	// TODO: horrible amount of code duplication
	public static boolean isElementInList(PartitionData data,
			ArrayList<PartitionData> partitionList, int elementIndex) {

		boolean exists = false;

		switch (elementIndex) {

		case BRANCH_RATE_MODEL_ELEMENT:

			int clockModelIndex = data.clockModelIndex;
			for (PartitionData data2 : partitionList) {

				if (clockModelIndex == data2.clockModelIndex) {

					for (int i = 0; i < data2.clockParameterValues.length; i++) {

						if (data.clockParameterValues[i] == data2.clockParameterValues[i]) {

							exists = true;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		case FREQUENCY_MODEL_ELEMENT:

			int frequencyModelIndex = data.frequencyModelIndex;
			for (PartitionData data2 : partitionList) {

				if (frequencyModelIndex == data2.frequencyModelIndex) {

					for (int i = 0; i < data2.frequencyParameterValues.length; i++) {

						if (data.frequencyParameterValues[i] == data2.frequencyParameterValues[i]) {

							exists = true;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		case BRANCH_MODEL_ELEMENT:

			int substitutionModelIndex = data.substitutionModelIndex;
			for (PartitionData data2 : partitionList) {

				if (substitutionModelIndex == data2.substitutionModelIndex) {

					for (int i = 0; i < data2.substitutionParameterValues.length; i++) {

						if (data.substitutionParameterValues[i] == data2.substitutionParameterValues[i]) {

							exists = true;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		case SITE_RATE_MODEL_ELEMENT:

			int siteRateModelIndex = data.siteRateModelIndex;
			for (PartitionData data2 : partitionList) {

				if (siteRateModelIndex == data2.siteRateModelIndex) {

					for (int i = 0; i < data2.siteRateModelParameterValues.length; i++) {

						if (data.siteRateModelParameterValues[i] == data2.siteRateModelParameterValues[i]) {

							exists = true;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		default:

			throw new RuntimeException("Unknown element");

		}// END: switch

		return exists;
	}// END: isModelInList

	// TODO: horrible amount of code duplication
	public static int isIdenticalWith(PartitionData data,
			ArrayList<PartitionData> partitionList, int elementIndex) {

		int index = -Integer.MAX_VALUE;

		switch (elementIndex) {

		case BRANCH_RATE_MODEL_ELEMENT:

			int clockModelIndex = data.clockModelIndex;
			for (PartitionData data2 : partitionList) {

				if (clockModelIndex == data2.clockModelIndex) {

					for (int i = 0; i < data2.clockParameterValues.length; i++) {

						if (data.clockParameterValues[i] == data2.clockParameterValues[i]) {

							index = partitionList.indexOf(data2);
							break;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		case FREQUENCY_MODEL_ELEMENT:

			int frequencyModelIndex = data.frequencyModelIndex;
			for (PartitionData data2 : partitionList) {

				if (frequencyModelIndex == data2.frequencyModelIndex) {

					for (int i = 0; i < data2.frequencyParameterValues.length; i++) {

						if (data.frequencyParameterValues[i] == data2.frequencyParameterValues[i]) {

							index = partitionList.indexOf(data2);
							break;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		case BRANCH_MODEL_ELEMENT:

			int substitutionModelIndex = data.substitutionModelIndex;
			for (PartitionData data2 : partitionList) {

				if (substitutionModelIndex == data2.substitutionModelIndex) {

					for (int i = 0; i < data2.substitutionParameterValues.length; i++) {

						if (data.substitutionParameterValues[i] == data2.substitutionParameterValues[i]) {

							index = partitionList.indexOf(data2);
							break;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		case SITE_RATE_MODEL_ELEMENT:

			int siteRateModelIndex = data.siteRateModelIndex;
			for (PartitionData data2 : partitionList) {

				if (siteRateModelIndex == data2.siteRateModelIndex) {

					for (int i = 0; i < data2.siteRateModelParameterValues.length; i++) {

						if (data.siteRateModelParameterValues[i] == data2.siteRateModelParameterValues[i]) {

							index = partitionList.indexOf(data2);
							break;

						}// END: parameters check
					}// END: parameters loop
				}// END: model index check
			}// END: list loop

			break;

		default:

			throw new RuntimeException("Unknown element");

		}// END: switch

		return index;
	}// END: isIdenticalWith

	// /////////////////
	// ---GUI UTILS---//
	// /////////////////

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
	}// END: uncaughtException

	private static void showExceptionDialog(Thread t, Throwable e) {

		String msg = String.format("Unexpected problem on thread %s: %s",
				t.getName(), e.getMessage());

		logException(t, e);

		JOptionPane.showMessageDialog(Utils.getActiveFrame(), //
				msg, //
				"Error", //
				JOptionPane.ERROR_MESSAGE, //
				BeagleSequenceSimulatorApp.errorIcon);
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
				BeagleSequenceSimulatorApp.errorIcon);
	}// END: showExceptionDialog

	private static void logException(Thread t, Throwable e) {
		e.printStackTrace();
	}// END: logException

	// ///////////////////////
	// ---DEBUGGING UTILS---//
	// ///////////////////////

	public static void printArray(int[] x) {
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

	public static void printPartitionData(PartitionData data) {

		System.out.println("\tTree model: " + data.treeFile);
		System.out.println("\tData type: "
				+ PartitionData.dataTypes[data.dataTypeIndex]);
		System.out.println("\tFrom: " + data.from);
		System.out.println("\tTo: " + data.to);
		System.out.println("\tEvery: " + data.every);
		System.out
				.println("\tBranch Substitution model: "
						+ PartitionData.substitutionModels[data.substitutionModelIndex]);
		System.out.println("\tSite rate model: "
				+ PartitionData.siteRateModels[data.siteRateModelIndex]);
		System.out.println("\tClock rate model: "
				+ PartitionData.clockModels[data.clockModelIndex]);
		System.out.println("\tFrequency model: "
				+ PartitionData.frequencyModels[data.frequencyModelIndex]);

	}// END: printPartitionData

	public static void printPartitionDataList(PartitionDataList dataList) {

		System.out.println("\tReplications: " + dataList.siteCount);

		int row = 1;
		for (PartitionData data : dataList) {

			System.out.println("Partition: " + row);

			printPartitionData(data);

			row++;
		}// END: data list loop

	}// END: printDataList

}// END: class
