package dr.app.bss;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;

public class Utils {

	// ///////////////////////////////
	// ---GENERAL UTILITY METHODS---//
	// ///////////////////////////////

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
	
	public static boolean taxonExists(Taxon taxon1, MutableTaxonList taxonList) {
		
		boolean exists = false;
		for(Taxon taxon2 : taxonList) {

			if(taxon1.equals(taxon2)) {
				exists = true;
				break;
			}
			
		}
		
		return exists;
	}// END: taxonExists
	
	public static boolean isModelInList(PartitionData data,
			ArrayList<PartitionData> partitionList) {

		boolean exists = false;
		
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

		return exists;
	}// END: isModelInList
	
	public static int isIdenticalWith(PartitionData data,
			ArrayList<PartitionData> partitionList) {

		int index = -Integer.MAX_VALUE;

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

	// ////////////////////////////////
	// ---EXCEPTION HANDLING UTILS---//
	// ////////////////////////////////

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

	public static void print2DArray(double[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(array[row][col] + " ");
			}
			System.out.print("\n");
		}
	}// END: print2DArray
	
	public static void printDataList(PartitionDataList dataList) {

		int row = 1;
		for (PartitionData data : dataList) {

			System.out.println("Partition: " + row);
			System.out.println("\tReplications: " + dataList.siteCount);
			System.out.println("\tFrom: " + data.from);
			System.out.println("\tTo: " + data.to);
			System.out.println("\tEvery: " + data.every);
			System.out.println("\tTree model: " + data.treeFile);
			System.out.println("\tData type: " + PartitionData.dataTypes[data.dataTypeIndex]);
			System.out.println("\tSubstitution model: "
					+ PartitionData.substitutionModels[data.substitutionModelIndex]);
			System.out.println("\tSite rate model: "
					+ PartitionData.siteRateModels[data.siteRateModelIndex]);
			System.out.println("\tClock rate model: "
					+ PartitionData.clockModels[data.clockModelIndex]);
			System.out.println("\tFrequency model: "
					+ PartitionData.frequencyModels[data.frequencyModelIndex]);

			// System.out.println("Possible trees: ");
			// for (int i = 0; i < dataList.treeFilesList.size(); i++) {
			// System.out.println(dataList.treeFilesList.get(i).getName());
			// }

			row++;
		}// END: data list loop

	}// END: printDataList

}// END: class
