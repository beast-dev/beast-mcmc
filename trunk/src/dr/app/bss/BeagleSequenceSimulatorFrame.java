package dr.app.bss;

import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.plaf.BorderUIResource;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.coalescent.CoalescentSimulatorParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evoxml.NewickParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

@SuppressWarnings("serial")
public class BeagleSequenceSimulatorFrame extends DocumentFrame {

	private JTabbedPane tabbedPane = new JTabbedPane();
	private TaxaPanel taxaPanel;
	private TreePanel treePanel;
	private PartitionsPanel partitionsPanel;
	private SimulationPanel simulationPanel;

	private PartitionDataList dataList;

	private JLabel statusLabel;
	private JProgressBar progressBar;
	private File workingDirectory = null;

	public BeagleSequenceSimulatorFrame(String title) {

		super();

		setTitle(title);
		dataList = new PartitionDataList();
		dataList.add(new PartitionData());

	}// END: Constructor

	@Override
	protected void initializeComponents() {

		try {

			setSize(new Dimension(900, 600));
			setMinimumSize(new Dimension(260, 100));

			taxaPanel = new TaxaPanel(this, dataList);
			treePanel = new TreePanel(this, dataList);
			partitionsPanel = new PartitionsPanel(this, dataList);
			simulationPanel = new SimulationPanel(this, dataList);

			tabbedPane.addTab("Taxa", null, taxaPanel);
			tabbedPane.addTab("Tree", null, treePanel);
			tabbedPane.addTab("Partitions", null, partitionsPanel);
			tabbedPane.addTab("Simulation", null, simulationPanel);

			statusLabel = new JLabel("No taxa loaded");

			JPanel progressPanel = new JPanel(new BorderLayout(0, 0));
			progressBar = new JProgressBar();
			progressPanel.add(progressBar, BorderLayout.CENTER);

			JPanel statusPanel = new JPanel(new BorderLayout(0, 0));
			statusPanel.add(statusLabel, BorderLayout.CENTER);
			statusPanel.add(progressPanel, BorderLayout.EAST);
			statusPanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
					new Insets(0, 6, 0, 6)));

			JPanel tabbedPanePanel = new JPanel(new BorderLayout(0, 0));
			tabbedPanePanel.add(tabbedPane, BorderLayout.CENTER);
			tabbedPanePanel.add(statusPanel, BorderLayout.SOUTH);
			tabbedPanePanel
					.setBorder(new BorderUIResource.EmptyBorderUIResource(
							new Insets(12, 12, 12, 12)));

			getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
			getContentPane().add(tabbedPanePanel, BorderLayout.CENTER);

			tabbedPane.setSelectedComponent(treePanel);

		} catch (NumberFormatException e) {
			e.printStackTrace();
		}

	}// END: initializeComponents

	public void fireTaxaChanged() {

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {

				taxaPanel.updateUI();
				setStatus(Integer.toString(dataList.taxonList.getTaxonCount())
						+ " taxa loaded.");

			}
		});

	}// END: fireTaxaChanged

	// ////////////////
	// ---SIMULATE---//
	// ////////////////

	public Action getExportAction() {
		return simulateAction;
	}// END: getExportAction

	private AbstractAction simulateAction = new AbstractAction("Simulate...") {
		public void actionPerformed(ActionEvent ae) {

			doExport();

		}// END: actionPerformed
	};

	public final void doExport() {

		try {

			if (dataList.forestMap.size() == 0) {

				tabbedPane.setSelectedComponent(treePanel);
				treePanel.doImport();

			} else {

				JFileChooser chooser = new JFileChooser();
				chooser.setDialogTitle("Simulate...");
				chooser.setMultiSelectionEnabled(false);
				chooser.setCurrentDirectory(workingDirectory);

				chooser.showSaveDialog(Utils.getActiveFrame());
				File file = chooser.getSelectedFile();

				if (file != null) {

					collectAllSettings();
					generateFile(file);

					File tmpDir = chooser.getCurrentDirectory();
					if (tmpDir != null) {
						workingDirectory = tmpDir;
					}

				}// END: file selected check

			}// END: tree loaded check

		} catch (Exception e) {
			Utils.handleException(e);
		} // END: try catch block

	}// END: doExport

	private void generateFile(final File outFile) throws IOException,
			ImportException {

		setBusy();

		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					if (BeagleSequenceSimulatorApp.DEBUG) {
						Utils.printDataList(dataList);
					}

					PrintWriter writer;

					writer = new PrintWriter(new FileWriter(outFile));

					ArrayList<Partition> partitionsList = new ArrayList<Partition>();

					for (PartitionData data : dataList) {

						// create partition
						Partition partition = new Partition(data.treeModel, //
								data.createBranchModel(), //
								data.createSiteRateModel(), //
								data.createBranchRateModel(), //
								data.createFrequencyModel(), //
								0, // from
								data.to - 1, // to
								1 // every
						);

						partitionsList.add(partition);

					}// END: data list loop

					BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
							partitionsList, dataList.sequenceLength);

					writer.println(beagleSequenceSimulator.simulate()
							.toString());
					writer.close();

				} catch (Exception e) {
					Utils.handleException(e);
				}

				return null;
			}// END: doInBackground

			// Executed in event dispatch thread
			public void done() {

				setStatus("Generated " + dataList.sequenceLength
						+ " replicates.");
				setIdle();

			}// END: done
		};

		worker.execute();

	}// END: generateFile

	// ////////////////////
	// ---GENERATE XML---//
	// ////////////////////

	public final void doGenerateXML() {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Generate XML...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(workingDirectory);

			chooser.showSaveDialog(Utils.getActiveFrame());
			File file = chooser.getSelectedFile();

			if (file != null) {

				collectAllSettings();
				generateXML(file);

				File tmpDir = chooser.getCurrentDirectory();
				if (tmpDir != null) {
					workingDirectory = tmpDir;
				}

			}// END: file selected check

		} catch (IOException e) {
			e.printStackTrace();
		}// END: try-catch block

	}// END: doGenerateXML

	public void generateXML(File file) throws IOException {

		// System.out.println("TODO: generateXML");

		XMLWriter writer = new XMLWriter(new BufferedWriter(
				new FileWriter(file)));

		// //////////////
		// ---header---//
		// //////////////

		writer.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>");
		writer.writeComment("Generated by "
				+ BeagleSequenceSimulatorApp.BEAGLE_SEQUENCE_SIMULATOR + " "
				+ BeagleSequenceSimulatorApp.VERSION);

		writer.writeOpenTag("beast");
		writer.writeBlankLine();

		// ////////////////////
		// ---taxa element---//
		// ////////////////////

		try {

			writeTaxa(dataList.taxonList, writer);
			writer.writeBlankLine();

		} catch (Exception e) {

			System.err.println(e);
			throw new RuntimeException("Taxon list generation has failed:\n"
					+ e.getMessage());

		}// END: try-catch block

		// /////////////////////////////
		// ---starting tree element---//
		// /////////////////////////////

		try {

			for (PartitionData data : dataList) {

				TreeModel tree = data.treeModel;
				writeStartingTree(tree, writer);
				writer.writeBlankLine();

			}// END: partitions loop

		} catch (Exception e) {

			System.err.println(e);
			throw new RuntimeException("Starting tree generation has failed:\n"
					+ e.getMessage());

		}// END: try-catch block

		// //////////////////////////
		// ---tree model element---//
		// //////////////////////////

		try {

			for (PartitionData data : dataList) {

				TreeModel tree = data.treeModel;
				writeTreeModel(tree, writer);
				writer.writeBlankLine();

			}// END: partitions loop

		} catch (Exception e) {

			System.err.println(e);
			throw new RuntimeException("Tree model generation has failed:\n"
					+ e.getMessage());

		}// END: try-catch block

		// //////////////////////////////////
		// ---branch rates model element---//
		// //////////////////////////////////

		try {

			for (PartitionData data : dataList) {

				writeBranchRatesModel(data, writer);
				writer.writeBlankLine();

			}// END: partitions loop

		} catch (Exception e) {

			System.err.println(e);
			throw new RuntimeException("Clock model generation has failed:\n"
					+ e.getMessage());

		}// END: try-catch block

		// ///////////////////////////////
		// ---frequency model element---//
		// ///////////////////////////////

		try {

			for (PartitionData data : dataList) {

				writeFrequencyModel(data, writer);
				writer.writeBlankLine();

			}// END: partitions loop

		} catch (Exception e) {

			System.err.println(e);
			throw new RuntimeException(
					"Frequency model generation has failed:\n" + e.getMessage());

		}// END: try-catch block

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		writer.writeCloseTag("beast");
		writer.flush();
		writer.close();
	}// END: generateXML

	private void writeFrequencyModel(PartitionData data, XMLWriter writer) {

//		 DataType dataType = data.
		
		
		
		
		
		
		
		
		
		
		
		
		
	}// END: writeFrequencyModel

	private void writeBranchRatesModel(PartitionData data, XMLWriter writer) {

		int clockModel = data.clockModel;
		switch (clockModel) {

		case 0: // StrictClock

			writer.writeOpenTag(
					StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
					new Attribute[] { new Attribute.Default<String>(
							XMLParser.ID, BranchRateModel.BRANCH_RATES) });

			writer.writeOpenTag("rate");

			ArrayList<Attribute.Default> attributes = new ArrayList<Attribute.Default>();

			attributes.add(new Attribute.Default<String>(XMLParser.ID,
					"clock.rate"));
			attributes.add(new Attribute.Default<String>(ParameterParser.VALUE,
					String.valueOf(data.clockParameterValues[0])));

			Attribute[] attrArray = new Attribute[attributes.size()];
			for (int i = 0; i < attrArray.length; i++) {
				attrArray[i] = attributes.get(i);
			}

			writer.writeTag(ParameterParser.PARAMETER, attrArray, true);

			writer.writeCloseTag("rate");

			writer.writeCloseTag(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES);

			break;

		}// END: switch

	}// END: writeBranchRatesModel

	private void writeTreeModel(TreeModel tree, XMLWriter writer) {

		final String treeModelName = TreeModel.TREE_MODEL;

		writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(
				XMLParser.ID, treeModelName), false);

		writer.writeIDref("tree", BeagleSequenceSimulatorApp.STARTING_TREE);

		writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
		writer.writeTag(ParameterParser.PARAMETER,
				new Attribute.Default<String>(XMLParser.ID, treeModelName + "."
						+ CoalescentSimulatorParser.ROOT_HEIGHT), true);
		writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);

		writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
				new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES,
						"true"));
		writer.writeTag(ParameterParser.PARAMETER,
				new Attribute.Default<String>(XMLParser.ID, treeModelName + "."
						+ "internalNodeHeights"), true);
		writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

		writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
				new Attribute[] {
						new Attribute.Default<String>(
								TreeModelParser.INTERNAL_NODES, "true"),
						new Attribute.Default<String>(
								TreeModelParser.ROOT_NODE, "true") });
		writer.writeTag(ParameterParser.PARAMETER,
				new Attribute.Default<String>(XMLParser.ID, treeModelName + "."
						+ "allInternalNodeHeights"), true);
		writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

		writer.writeCloseTag(TreeModel.TREE_MODEL);

	}// END: writeTreeModel

	private void writeTaxa(TaxonList taxonList, XMLWriter writer) {

		writer.writeOpenTag(TaxaParser.TAXA, // tagname
				new Attribute[] { // attributes[]
				new Attribute.Default<String>(XMLParser.ID, TaxaParser.TAXA) });

		for (int i = 0; i < taxonList.getTaxonCount(); i++) {

			Taxon taxon = taxonList.getTaxon(i);

			writer.writeTag(
					TaxonParser.TAXON, // tagname
					new Attribute[] { // attributes[]
					new Attribute.Default<String>(XMLParser.ID, taxon.getId()) },
					true // close
			);

		}// END: i loop

		writer.writeCloseTag(TaxaParser.TAXA);
	}// END: writeTaxa

	private void writeStartingTree(TreeModel tree, XMLWriter writer) {

		writer.writeOpenTag(NewickParser.NEWICK,
				new Attribute[] { new Attribute.Default<String>(XMLParser.ID,
						BeagleSequenceSimulatorApp.STARTING_TREE),
				// new Attribute.Default<String>(DateParser.UNITS,
				// options.datesUnits.getAttribute()),
				// new Attribute.Default<Boolean>(SimpleTreeParser.USING_DATES,
				// options.clockModelOptions.isTipCalibrated())
				});
		writer.writeText(Tree.Utils.newick(tree));
		writer.writeCloseTag(NewickParser.NEWICK);

	}// END: writeStartingTree

	@Override
	public JComponent getExportableComponent() {
		JComponent exportable = null;
		Component component = tabbedPane.getSelectedComponent();

		if (component instanceof Exportable) {
			exportable = ((Exportable) component).getExportableComponent();
		} else if (component instanceof JComponent) {
			exportable = (JComponent) component;
		}

		return exportable;
	}// END: getExportableComponent

	@Override
	protected boolean readFromFile(File arg0) throws IOException {
		return false;
	}

	@Override
	protected boolean writeToFile(File arg0) throws IOException {
		return false;
	}

	private void collectAllSettings() {

		// frequencyPanel.collectSettings();
		// substModelPanel.collectSettings();
		// clockPanel.collectSettings();
		// sitePanel.collectSettings();
		simulationPanel.collectSettings();

	}// END: collectAllSettings

	public void dataSelectionChanged(boolean isSelected) {
		if (isSelected) {
			getDeleteAction().setEnabled(true);
		} else {
			getDeleteAction().setEnabled(false);
		}
	}// END: dataSelectionChanged

	public File getWorkingDirectory() {
		return workingDirectory;
	}// END: getWorkingDirectory

	public void setWorkingDirectory(File workingDirectory) {
		this.workingDirectory = workingDirectory;
	}// END: setWorkingDirectory

	public void fireModelChanged() {
		collectAllSettings();
	}// END: fireModelChanged

	public void setBusy() {
		progressBar.setIndeterminate(true);
	}

	public void setIdle() {
		progressBar.setIndeterminate(false);
	}

	public void setStatus(String status) {
		statusLabel.setText(status);
	}

}// END: class
