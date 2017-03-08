/*
 * Utils.java
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

import java.awt.Desktop;
import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.charset.Charset;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import dr.evolution.tree.TreeUtils;
import org.apache.commons.math.random.MersenneTwister;

import dr.app.bss.test.AncestralSequenceTrait;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class Utils {

	// ////////////////////
	// ---THREAD UTILS---//
	// ////////////////////

	public static void sleep(int seconds) {
		try {

			Thread.sleep(seconds * 1000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}// END: sleep
	
	// ////////////////////////////////
	// ---RANDOM NUMB3R GENERATION---//
	// ////////////////////////////////
	
	private static MersenneTwister random = new MersenneTwister(
			MathUtils.nextLong());

    public static double rLogNormal(double stdev, double mean) {
    	
    	double rNorm = random.nextGaussian() * stdev + mean;
    	double rLognormal = Math.exp(rNorm);
    	
    	return rLognormal;
	}// END: drawRandom
	
	public static int rMultinom(double[] probabilities) {

		int range = probabilities.length + 1;
		double[] distribution = new double[range];
		double sumProb = 0;

		for (double value : probabilities) {
			sumProb += value;
		}// END: probabilities loop

		distribution[0] = 0;
		for (int i = 1; i < range; ++i) {

			distribution[i] = distribution[i - 1]
					+ (probabilities[i - 1] / sumProb);

		}// END: i loop

		distribution[range - 1] = 1.0;

		double key = random.nextDouble();

		int mindex = 1;
		int maxdex = range - 1;
		int midpoint = mindex + (maxdex - mindex) / 2;
		while (mindex <= maxdex) {

			if (key < distribution[midpoint - 1]) {
				
				maxdex = midpoint - 1;
				
			} else if (key > distribution[midpoint]) {
				
				mindex = midpoint + 1;
				
			} else {
				
				return midpoint - 1;

			}
			
			midpoint = mindex + (int) Math.ceil((maxdex - mindex) / 2);

		}//END: mindex loop
		
		System.out.println("Error in rMultinom!");
		
		return range - 1;
	}//END: rMultinom
    
	// ////////////
	// ---MATH---//
	// ////////////
    
	public static int sample(double[] probabilities) {

		int samplePos = -Integer.MAX_VALUE;
		double cumProb = 0.0;
		double u = random.nextDouble();

		for (int i = 0; i < probabilities.length; i++) {

			cumProb += probabilities[i];

			if (u <= cumProb) {
				samplePos = i;
				break;
			}
		}

		return samplePos;
	}// END: randomChoicePDF
	
    public static void rescale(double[] logX) {
        double max = max(logX);
        for (int i = 0; i < logX.length; i++) {
            logX[i] -= max;
		}
	}// END: rescale

	public double getParameterVariance(Parameter param) {
		
		int n = param.getSize();
		double mean = getParameterMean(param);
		double var = 0;
		
		for (int i = 0; i < n; i++) {

			var+= Math.pow( (param.getValue(i) - mean), 2);

		}
		
		var/= (n-1);
		
		return var;
	}// END: getParameterVariance

	public double getParameterMean(Parameter param) {
		double mean = 0;
		int n = param.getSize();
		for (int i = 0; i < n; i++) {
			mean += param.getValue(i);
		}
		mean /= n;

		return mean;
	}// END: getParameterMean
	
	public static double getNorm(double[] vector) {
		double norm = 0;
		for (int i = 0; i < vector.length; i++) {
			norm += Math.pow(vector[i], 2);
		}

		return Math.sqrt(norm);
	}// END: getNorm

	public static void normalize(double[] vector) {

		double norm = getNorm(vector);
		for (int i = 0; i < vector.length; i++) {
			vector[i] = vector[i] / norm;
		}

	}// END: normalize
	
	// ////////////////////
	// ---ARRAYS UTILS---//
	// ////////////////////

	public static void exponentiate(double[] array) {
		for (int i = 0; i < array.length; i++) {
			array[i] = Math.exp(array[i]);
		}
	}// END: exponentiate
    
	public static int max(int[] array) {

		int max = -Integer.MAX_VALUE;

		for (int i = 0; i < array.length; i++) {
		
			if (array[i] > max) {

				max = (int)array[i];

			}// END: if check

		}// END: i loop

		return max;
	}// END: findMaximum
	
	public static int max(double[] array) {

		int max = -Integer.MAX_VALUE;

		for (int i=0; i< array.length;i++) {
		
			if (array[i] > max) {

				max = (int)array[i];

			}// END: if check

		}// END: i loop

		return max;
	}// END: findMaximum
	
	public static int max(ArrayList<Integer> array) {

		int max = -Integer.MAX_VALUE;

		for (Integer element : array) {

			if (element > max) {

				max = element;

			}// END: if check

		}// END: i loop

		return max;
	}// END: findMaximum

	public static double sumArray(int[] array) {

		double sum = 0.0;
		for (int i = 0; i < array.length; i++) {

			sum += array[i];

		}

		return sum;
	}// END: sumArray
    
	public static double sumArray(double[] array) {

		double sum = 0.0;
		for (int i = 0; i < array.length; i++) {

			sum += array[i];

		}

		return sum;
	}// END: sumArray
	
	// /////////////////
	// ---CONSTANTS---//
	// /////////////////

	// public static final int TREE_MODEL_ELEMENT = 0;
	public static final int BRANCH_MODEL_ELEMENT = 1;
	public static final int SITE_RATE_MODEL_ELEMENT = 2;
	public static final int BRANCH_RATE_MODEL_ELEMENT = 3;
	public static final int FREQUENCY_MODEL_ELEMENT = 4;
	public static final int DEMOGRAPHIC_MODEL_ELEMENT = 5;

	public static final String TOPOLOGY = "topology";
	public static final String ABSOLUTE_HEIGHT = "absoluteHeight";
	public static final String TREE_FILENAME = "treeFilename";
	public static final String SUBSTITUTION_MODEL = "substitutionModel";
	public static final String DEMOGRAPHIC_MODEL = "demographicModel";
	public static final String FREQUENCY_MODEL = "frequencyModel";
	public static final String CODON_UNIVERSAL = "codon-universal";
	public static final String CHOOSE_FILE = "Choose file...";
	public static final String EDIT_TAXA_SET = "Edit taxa set...";
	public static final String ANCESTRAL_SEQUENCE = "ancestralSequence";

	public static final String BSS_ICON = "icons/bss.png";
	public static final String CHECK_ICON = "icons/check.png";
	public static final String ERROR_ICON = "icons/error.png";
	public static final String HAMMER_ICON = "icons/hammer.png";
	public static final String CLOSE_ICON = "icons/close.png";
	public static final String BIOHAZARD_ICON = "icons/biohazard.png";
	public static final String BUBBLE_BLUE_ICON = "icons/bubble-blue.png";
	public static final String SAVE_ICON = "icons/save.png";
	public static final String TEXT_FILE_ICON = "icons/file.png";

	public static final double[] UNIFORM_CODON_FREQUENCIES = new double[] {
			0.0163936, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344,
			0.01639344 };

	public static final String STOP_CODONS[] = new String[] { "TAA", "TAG", "TGA" };
	
	// ///////////////////////////////
	// ---GENERAL UTILITY METHODS---//
	// ///////////////////////////////

	public static double logfactor(int n) {

		double logfactor = 0.0;
		for (int i = 1; i <= n; i++) {
			logfactor += Math.log(i);
		}

		return logfactor;
	}

	public static double map(double value, double low1, double high1,
			double low2, double high2) {
		/**
		 * maps a single value from its range into another interval
		 * 
		 * @param low1
		 *            , high1 - range of value; low2, high2 - interval
		 * @return the mapped value
		 */
		return (value - low1) / (high1 - low1) * (high2 - low2) + low2;
	}// END: map

	public static String[] loadStrings(String filename) throws IOException {

		int linesCount = countLines(filename);
		String[] lines = new String[linesCount];

		FileInputStream inputStream;
		BufferedReader reader;
		String line;

		inputStream = new FileInputStream(filename);
		reader = new BufferedReader(new InputStreamReader(inputStream,
				Charset.forName("UTF-8")));
		int i = 0;
		while ((line = reader.readLine()) != null) {
			lines[i] = line;
			i++;
		}

		// Clean up
		reader.close();
		reader = null;
		inputStream = null;

		return lines;
	}// END: loadStrings

	public static int countLines(String filename) throws IOException {

		InputStream is = new BufferedInputStream(new FileInputStream(filename));

		byte[] c = new byte[1024];
		int count = 0;
		int readChars = 0;
		boolean empty = true;
		while ((readChars = is.read(c)) != -1) {
			empty = false;
			for (int i = 0; i < readChars; ++i) {
				if (c[i] == '\n') {
					++count;
				}
			}
		}

		is.close();

		return (count == 0 && !empty) ? 1 : count;

	}// END: countLines

	public static Taxa importTaxaFromFile(File file) throws IOException {

		Taxa taxa = new Taxa();
		Taxon taxon;
		String[] lines = Utils.loadStrings(file.getAbsolutePath());

		for (int i = 0; i < lines.length; i++) {

			String[] line = lines[i].split("\\s+");

			taxon = new Taxon(line[TaxaEditorTableModel.NAME_INDEX]);
			taxon.setAttribute(Utils.ABSOLUTE_HEIGHT,
					Double.valueOf(line[TaxaEditorTableModel.HEIGHT_INDEX]));

			taxa.addTaxon(taxon);

		}// END: i loop

		return taxa;
	}// END: importTaxaFromFile

	public static Tree importTreeFromFile(File file) throws IOException,
			ImportException {

		Tree tree = null;

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

		return tree;
	}// END: importTreeFromFile

	public static void removeTaxaWithAttributeValue(PartitionDataList dataList,
			String attribute, String value) {

		for (int i = 0; i < dataList.allTaxa.getTaxonCount(); i++) {

			Taxon taxon = dataList.allTaxa.getTaxon(i);
			if (taxon.getAttribute(attribute).toString()
					.equalsIgnoreCase(value)) {

				dataList.allTaxa.removeTaxon(taxon);
				i--;

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
			// &&
			// taxon.getAttribute(Utils.TREE_FILENAME).toString().equalsIgnoreCase(taxon2.getAttribute(Utils.TREE_FILENAME).toString())
			) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: taxonExists

	// private boolean isFileInList(File file) {
	// boolean exists = false;
	//
	// for (File file2 : dataList.treesList) {
	//
	// if (file.getName().equalsIgnoreCase(file2.getName())) {
	// exists = true;
	// break;
	// }
	//
	// }
	//
	// return exists;
	// }// END: isFileInList

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

	public static boolean isTaxaInList(Taxa taxa, ArrayList<Taxa> taxaList) {

		boolean exists = false;

		for (Taxa taxa2 : taxaList) {

			if (taxaToString(taxa, true).equalsIgnoreCase(
					taxaToString(taxa2, true))) {
				exists = true;
				break;
			}

		}

		return exists;
	}// END: isTaxaInList

	public static int taxaIsIdenticalWith(Taxa taxa, ArrayList<Taxa> taxaList) {

		int index = -Integer.MAX_VALUE;

		for (Taxa taxa2 : taxaList) {

			if (taxaToString(taxa, true).equalsIgnoreCase(
					taxaToString(taxa2, true))) {
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

		URL imgURL = BeagleSequenceSimulatorApp.class.getResource(path);

		if (imgURL != null) {
			icon = new ImageIcon(imgURL);
		} else {
			System.err.println("Couldn't find file: " + path + "\n");
		}

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

			} else if (obj instanceof double[]) {

				double[] seq = (double[]) obj;
				System.out.print(pairs.getKey() + " =");

				for (int i = 0; i < seq.length; ++i) {
					System.out.print(" " + seq[i]);
				}
				System.out.println();

			} else {
				System.out.println(pairs.getKey() + " = " + pairs.getValue());
			}// END: obj class check

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
			System.out.print(x[i] + " ");
		}
		System.out.println();
	}// END: printArray

	public static void printArray(int[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.print(x[i] + " ");
		}
		System.out.println();
	}// END: printArray

	public static void printArray(double[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.print(x[i] + " ");
		}
		System.out.println();
	}// END: printArray

	public static void printArray(boolean[] x) {
		for (int i = 0; i < x.length; i++) {
			System.out.print(x[i] + " ");
		}
		System.out.println();
		
	}
	
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

	public static void print2DArray(int[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {
				System.out.print(array[row][col] + " ");
			}
			System.out.print("\n");
		}
	}// END: print2DArray
	
	public static void print2Arrays(int[] array1, double[] array2, int nrow) {
		for (int row = 0; row < nrow; row++) {
			System.out.print(array1[row] + " " + array2[row] + " ");
			System.out.print("\n");
		}
	}// END: print2DArray

	public static void print2DArray(double[][] array, int formatEvery) {

		int i = 0;
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[row].length; col++) {

				if (i == formatEvery) {
					System.out.print("\n");
					i = 0;
				}

				System.out.print(array[row][col] + " ");
				i++;
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

	public static void printDemographicModel(PartitionData data) {
		System.out.print("\tDemographic model: ");
		System.out.print(demographicModelToString(data));
		System.out.print("\n");
	}// END: printFrequencyModel

	private static void printDataType(PartitionData data) {
		System.out.print("\tData type: ");
		System.out.print(dataTypeToString(data));
		System.out.print("\n");
	}// END: printDataType

	public static void printTaxaSet(Taxa taxa) {
		for (int i = 0; i < taxa.getTaxonCount(); i++) {
			Taxon taxon = taxa.getTaxon(i);
			System.out.print("\t\t " + taxonToString(taxon, false) + ("\n"));
		}
	}// END: printTaxaSet

	public static void printTree(TreesTableRecord record) {
		System.out.print(record.getTree().toString());
		System.out.print("\n");
	}// END: printTree

	public static void printRecord(TreesTableRecord record) {

		if (record == null) {

			System.out.println("\tRecord: NOT SET");

		} else if (record.isTreeSet()) {

			System.out.print("\t" + record.getName() + ": ");
			printTree(record);

		} else if (record.isTaxaSet()) {

			System.out.println("\t" + record.getName() + ":");
			printTaxaSet(record.getTaxa());

		} else {
			//
		}
	}// END: printRecord

	public static void printRecords(PartitionDataList dataList) {
		for (TreesTableRecord record : dataList.recordsList) {
			printRecord(record);
		}// END: record loop
	}// END: printRecords

	public static void printPartitionData(PartitionData data) {
		printRecord(data.record);
		printDataType(data);
		printDemographicModel(data);
		System.out.println("\tFrom: " + data.from);
		System.out.println("\tTo: " + data.to);
		System.out.println("\tEvery: " + data.every);
		printBranchSubstitutionModel(data);
		printSiteRateModel(data);
		printClockRateModel(data);
		printFrequencyModel(data);
	}// END: printPartitionData

	public static void printPartitionDataList(PartitionDataList dataList) {

		// System.out.println(dataList.get(0).from + " " +
		// dataList.get(1).from);

		if (BeagleSequenceSimulatorApp.DEBUG) {
			System.out.println("Possible records: ");
			printRecords(dataList);
		}

		System.out.println("\tSite count: " + getSiteCount(dataList));
		System.out.println("\tOutput type: " + dataList.outputFormat);

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

	public static void printTaxonList(PartitionDataList dataList) {
		System.out.println(taxaToString(dataList.allTaxa, true));
	}// END: printTaxonList

	public static Sequence intArray2Sequence(Taxon taxon, int[] seq,
			int gapFlag, DataType dataType) {

		StringBuilder sSeq = new StringBuilder();
		int partitionSiteCount = seq.length;

		if (dataType instanceof Codons) {

			for (int i = 0; i < partitionSiteCount; i++) {

				int state = seq[i];

				if (state == gapFlag) {
					sSeq.append(dataType.getTriplet(dataType.getGapState()));
				} else {
					sSeq.append(dataType.getTriplet(seq[i]));
				}// END: gap check

			}// END: replications loop

		} else {

			for (int i = 0; i < partitionSiteCount; i++) {

				int state = seq[i];

				if (state == gapFlag) {
					sSeq.append(dataType.getCode(dataType.getGapState()));
				} else {
					sSeq.append(dataType.getCode(seq[i]));
				}// END: gap check

			}// END: replications loop

		}// END: dataType check

		return new Sequence(taxon, sSeq.toString());
	}// END: intArray2Sequence

	// //////////////////////
	// ---TOSTRING UTILS---//
	// //////////////////////

	public static String taxonToString(Taxon taxon, boolean printNames) {

		String string = null;

		if (printNames) {
			string = taxon.getId() + " ("
					+ taxon.getAttribute(Utils.ABSOLUTE_HEIGHT) + ","
					+ taxon.getAttribute(Utils.TREE_FILENAME) + ")";

		} else {
			string = taxon.getId() + " ("
					+ taxon.getAttribute(Utils.ABSOLUTE_HEIGHT) + ")";
		}

		return string;
	}// END: taxonToString

	public static String taxaToString(Taxa taxa, boolean printNames) {

		String string = "";

		for (int i = 0; i < taxa.getTaxonCount(); i++) {

			Taxon taxon = taxa.getTaxon(i);
			string += taxonToString(taxon, printNames) + ("\n");

		}

		return string;
	}// END: taxaToString

	public static String partitionDataToString(PartitionData data,
			TreeModel simulatedTreeModel
	// , LinkedHashMap<NodeRef, int[]> sequencesMap
	) {

		String string = "";

		// if (data.record.isTreeSet()) {
		//
		// string += ("Tree: " + data.record.getTree().toString())+ ("\n");
		//
		// } else if (data.record.isTaxaSet()) {
		//
		// string += ("Taxa Set: \n" + taxaToString(data.record.getTaxa(),
		// false));//+ ("\n");
		//
		// } else {
		// //
		// }

		string += ("Tree model: " + simulatedTreeModel.toString()) + ("\n");
		// string += ("Tree model: "
		// +annotatedTreeModelToString(simulatedTreeModel, sequencesMap,
		// data.createDataType()) ) + ("\n");
		string += ("From: " + data.from) + ("\n");
		string += ("To: " + data.to) + ("\n");
		string += ("Every: " + data.every) + ("\n");
		string += ("Data type: ") + dataTypeToString(data) + ("\n");
		string += ("Demographic model: ") + demographicModelToString(data)
				+ ("\n");
		string += ("Branch Substitution model: ")
				+ branchSubstitutionModelToString(data) + ("\n");
		string += ("Frequency model: ") + frequencyModelToString(data) + ("\n");
		string += ("Site Rate model: ") + siteRateModelToString(data) + ("\n");
		string += ("Clock Rate model: ") + clockRateModelToString(data)
				+ ("\n");

		return string;
	}// END: partitionDataToString

	public static String partitionDataListToString(PartitionDataList dataList, //
			ArrayList<TreeModel> simulatedTreeModelList
	// ,LinkedHashMap<Integer,LinkedHashMap<NodeRef, int[]>>
	// partitionSequencesMap
	) {

		String string = "";
		TreeModel simulatedTreeModel;
		// LinkedHashMap<NodeRef, int[]> sequencesMap;

		string += ("Site count: " + getSiteCount(dataList)) + ("\n");
		if (dataList.setSeed) {
			string += ("Starting seed: " + dataList.startingSeed) + ("\n");
		}

		int row = 0;
		for (PartitionData data : dataList) {

			simulatedTreeModel = simulatedTreeModelList.get(row);
			// sequencesMap = partitionSequencesMap.get(row);

			string += ("Partition: " + (row + 1)) + ("\n");
			string += partitionDataToString(data, simulatedTreeModel
			// , sequencesMap
			);
			string += ("\n");
			row++;

		}// END: data list loop

		return string;
	}// END: partitionDataListToString

	// TODO: doesn't work
	public static String annotatedTreeModelToString(TreeModel treeModel,
			LinkedHashMap<NodeRef, int[]> sequencesMap, DataType dataType) {

		StringBuffer buffer = new StringBuffer();
		NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
		boolean useTipLabels = true;

		AncestralSequenceTrait ancestralSequence = new AncestralSequenceTrait(
				sequencesMap, dataType);
		TreeTraitProvider[] treeTraitProviders = new TreeTraitProvider[] { ancestralSequence };

		TreeUtils.newick(treeModel, //
				treeModel.getRoot(), //
				useTipLabels, //
				TreeUtils.BranchLengthType.LENGTHS_AS_TIME, //
				format, //
				null, //
				treeTraitProviders, //
				null, buffer);

		return buffer.toString();
	}

	private static String dataTypeToString(PartitionData data) {
		String string = PartitionData.dataTypes[data.dataTypeIndex];
		return string;
	}

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

		if(data.clockModelIndex == data.LRC_INDEX) {
			
			String space = (data.lrcParametersInRealSpace == true ? "real" : "log");
			string += " ( " +  "Parameters in " + space + " space )";
			
		}
		
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

	public static String demographyModelToString(PartitionData data) {

		String string = PartitionData.demographicModels[data.demographicModelIndex];

		string += (" ( ");
		for (int i = 0; i < PartitionData.demographicParameterIndices[data.demographicModelIndex].length; i++) {
			string += data.demographicParameterValues[PartitionData.demographicParameterIndices[data.demographicModelIndex][i]];
			string += " ";
		}// END: indices loop
		string += ")";

		return string;
	}


}// END: class
