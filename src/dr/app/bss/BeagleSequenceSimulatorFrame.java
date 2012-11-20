package dr.app.bss;

import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.event.ActionEvent;
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
import dr.evolution.io.Importer.ImportException;

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
			simulationPanel = new SimulationPanel(
					this, 
					dataList);

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
			tabbedPanePanel.setBorder(new BorderUIResource.EmptyBorderUIResource(
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
				setStatus(Integer.toString(dataList.taxonList.getTaxonCount()) + " taxa loaded.");

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

					writer.println(beagleSequenceSimulator.simulate().toString());
					writer.close();

				} catch (Exception e) {
					Utils.handleException(e);
				}

				return null;
			}// END: doInBackground

			// Executed in event dispatch thread
			public void done() {

				setStatus("Generated " + dataList.sequenceLength + " replicates.");
				setIdle();

			}// END: done
		};

		worker.execute();

	}// END: generateFile

	
	
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
