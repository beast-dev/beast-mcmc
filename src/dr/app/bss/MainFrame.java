package dr.app.bss;

import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

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
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evomodel.tree.TreeModel;

@SuppressWarnings("serial")
public class MainFrame extends DocumentFrame {

	private JTabbedPane tabbedPane = new JTabbedPane();
	private TaxaPanel taxaPanel;
	private TreePanel treePanel;
	private PartitionsPanel partitionsPanel;
	private SimulationPanel simulationPanel;

	private PartitionDataList dataList;

	private JLabel statusLabel;
	private JProgressBar progressBar;
	private File workingDirectory = null;

	public MainFrame(String title) {

		super();

		setTitle(title);
		dataList = new PartitionDataList();
		dataList.add(new PartitionData());

	}// END: Constructor

	@Override
	protected void initializeComponents() {

		try {

			setSize(new Dimension(1100, 600));
			setMinimumSize(new Dimension(260, 100));

			taxaPanel = new TaxaPanel(dataList);
			treePanel = new TreePanel(this, dataList);
			partitionsPanel = new PartitionsPanel(this, dataList);
			simulationPanel = new SimulationPanel(this, dataList);

			tabbedPane.addTab("Taxa", null, taxaPanel);
			tabbedPane.addTab("Trees", null, treePanel);
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

	// ////////////////
	// ---SIMULATE---//
	// ////////////////

	// file chooser
	public void doExport() {

		if (dataList.forestMap.size() == 0
				&& simulationPanel.simulationType == SimulationPanel.FIRST_SIMULATION_TYPE) {

			tabbedPane.setSelectedComponent(treePanel);
			treePanel.doImportTree();

		} else if (dataList.treesFilename == null
				&& simulationPanel.simulationType == SimulationPanel.SECOND_SIMULATION_TYPE) {

			tabbedPane.setSelectedComponent(treePanel);
			treePanel.doSelectTreesFilename();

		} else {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Simulate...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(workingDirectory);

			int returnVal = chooser.showSaveDialog(Utils.getActiveFrame());
			if (returnVal == JFileChooser.APPROVE_OPTION) {

				File file = chooser.getSelectedFile();

				collectAllSettings();

				switch (simulationPanel.simulationType) {

				case SimulationPanel.FIRST_SIMULATION_TYPE:

					generateNumberOfSimulations(file);
					break;

				case SimulationPanel.SECOND_SIMULATION_TYPE:
					generateForEachTree(file);
					break;

				default:
					throw new RuntimeException("Unknown analysis type!");
				}// END: switch

				File tmpDir = chooser.getCurrentDirectory();
				if (tmpDir != null) {
					workingDirectory = tmpDir;
				}

			}// END: approve check

		}// END: tree loaded check

	}// END: doExport

	// threading, UI, exceptions handling
	private void generateForEachTree(final File outFile) {

		setBusy();

		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					Partition partition;
					ArrayList<Partition> partitionsList;
					PrintWriter writer;
					TreeModel treeModel;
					BeagleSequenceSimulator beagleSequenceSimulator;
					TreeImporter importer = null;

					BufferedReader reader = new BufferedReader(new FileReader(
							dataList.treesFilename));

					String line = reader.readLine();
					if (line.toUpperCase().startsWith("#NEXUS")) {

						importer = new NexusImporter(reader);

					} else {

						importer = new NewickImporter(reader);

					}// END: nexus check

					int treesRead = 0;
					while (importer.hasTree()) {

						setStatus("Generating for tree # " + treesRead + 1);

						treeModel = new TreeModel(importer.importNextTree());

						String path = ((treesRead == 0) ? outFile.toString()
								: outFile.toString() + treesRead);

						writer = new PrintWriter(new FileWriter(path));

						partitionsList = new ArrayList<Partition>();

						for (PartitionData data : dataList) {

							// create partition
							partition = new Partition(treeModel, //
									data.createBranchModel(), //
									data.createSiteRateModel(), //
									data.createBranchRateModel(), //
									data.createFrequencyModel(), //
									0, // from
									data.to - 1, // to
									1 // every
							);

							partitionsList.add(partition);

						}// END: partition loop

						beagleSequenceSimulator = new BeagleSequenceSimulator(
								partitionsList, dataList.siteCount);

						writer.println(beagleSequenceSimulator.simulate()
								.toString());
						writer.close();

						treesRead++;
					}// END: trees loop

				} catch (Exception e) {
					Utils.handleException(e);
				}// END: try-catch block

				return null;
			}// END: doInBackground

			// Executed in event dispatch thread
			public void done() {

				setStatus("Finished.");
				setIdle();

			}// END: done
		};

		worker.execute();

	}// END: generateForEachTree

	// threading, UI, exceptions handling
	private void generateNumberOfSimulations(final File outFile) {

		setBusy();

		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					for (int i = 0; i < dataList.simulationsCount; i++) {

						if (BeagleSequenceSimulatorApp.DEBUG) {
							Utils.printDataList(dataList);
						}

						String path = ((i == 0) ? outFile.toString() : outFile
								.toString() + i);

						PrintWriter writer = new PrintWriter(new FileWriter(
								path));

						ArrayList<Partition> partitionsList = new ArrayList<Partition>();

						for (PartitionData data : dataList) {

							if (data.treeModel == null) {

								throw new RuntimeException(
										"Set Tree Model in Partitions tab for "
												+ (partitionsList.size() + 1)
												+ " partition.");

							} else {

								// create partition
								Partition partition = new Partition(
										data.treeModel, //
										data.createBranchModel(), //
										data.createSiteRateModel(), //
										data.createBranchRateModel(), //
										data.createFrequencyModel(), //
										0, // from
										data.to - 1, // to
										1 // every
								);

								partitionsList.add(partition);

							}

						}// END: data list loop

						BeagleSequenceSimulator beagleSequenceSimulator = new BeagleSequenceSimulator(
								partitionsList, dataList.siteCount);

						writer.println(beagleSequenceSimulator.simulate()
								.toString());
						writer.close();

					}// END: simulationsCount loop

				} catch (Exception e) {
					Utils.handleException(e);
				}

				return null;
			}// END: doInBackground

			// Executed in event dispatch thread
			public void done() {

				setStatus("Generated " + dataList.siteCount + " replicates.");
				setIdle();

			}// END: done
		};

		worker.execute();

	}// END: generateNumberOfSimulations

	// ////////////////////
	// ---GENERATE XML---//
	// ////////////////////

	// file chooser
	public final void doGenerateXML() {

		JFileChooser chooser = new JFileChooser();
		chooser.setDialogTitle("Generate XML...");
		chooser.setMultiSelectionEnabled(false);
		chooser.setCurrentDirectory(workingDirectory);

		int returnVal = chooser.showSaveDialog(Utils.getActiveFrame());
		if (returnVal == JFileChooser.APPROVE_OPTION) {

			File file = chooser.getSelectedFile();

			generateXML(file);

			File tmpDir = chooser.getCurrentDirectory();
			if (tmpDir != null) {
				workingDirectory = tmpDir;
			}

		}// END: approve check

	}// END: doGenerateXML

	private void generateXML(final File outFile) {

		setBusy();

		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					collectAllSettings();
					XMLGenerator xmlGenerator = new XMLGenerator(dataList);
					xmlGenerator.generateXML(outFile);

				} catch (Exception e) {
					Utils.handleException(e);
				}

				return null;
			}// END: doInBackground

			// Executed in event dispatch thread
			public void done() {

				setStatus("Generated " + outFile);
				setIdle();

			}// END: done
		};

		worker.execute();

	}// END: generateNumberOfSimulations

	// /////////////////
	// ---MAIN MENU---//
	// /////////////////

	// public Action getExportAction() {
	// return simulateAction;
	// }// END: getExportAction
	//
	// private AbstractAction simulateAction = new AbstractAction("Simulate...")
	// {
	// public void actionPerformed(ActionEvent ae) {
	//
	// doExport();
	//
	// }// END: actionPerformed
	// };

	@Override
	protected boolean readFromFile(File arg0) throws IOException {
		return false;
	}

	@Override
	protected boolean writeToFile(File arg0) throws IOException {
		return false;
	}

	// //////////////////////
	// ---SHARED METHODS---//
	// //////////////////////

	public File getWorkingDirectory() {
		return workingDirectory;
	}// END: getWorkingDirectory

	public void setWorkingDirectory(File workingDirectory) {
		this.workingDirectory = workingDirectory;
	}// END: setWorkingDirectory

	public void collectAllSettings() {

		// frequencyPanel.collectSettings();
		// substModelPanel.collectSettings();
		// clockPanel.collectSettings();
		// sitePanel.collectSettings();
		simulationPanel.collectSettings();

	}// END: collectAllSettings

	public void fireTaxaChanged() {

		if (SwingUtilities.isEventDispatchThread()) {

			taxaPanel.fireTableDataChanged();
			setStatus(Integer.toString(dataList.taxonList.getTaxonCount())
					+ " taxa loaded.");

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					taxaPanel.fireTableDataChanged();
					setStatus(Integer.toString(dataList.taxonList
							.getTaxonCount()) + " taxa loaded.");

				}
			});
		}// END: edt check

	}// END: fireTaxaChanged

	public void setBusy() {

		if (SwingUtilities.isEventDispatchThread()) {

			simulationPanel.setBusy();
			progressBar.setIndeterminate(true);

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					simulationPanel.setBusy();
					progressBar.setIndeterminate(true);

				}
			});
		}// END: edt check

	}// END: setBusy

	public void setIdle() {

		if (SwingUtilities.isEventDispatchThread()) {

			simulationPanel.setIdle();
			progressBar.setIndeterminate(false);

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					simulationPanel.setIdle();
					progressBar.setIndeterminate(false);

				}
			});
		}// END: edt check

	}// END: setIdle

	public void setStatus(final String status) {

		if (SwingUtilities.isEventDispatchThread()) {

			statusLabel.setText(status);

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					statusLabel.setText(status);

				}
			});
		}// END: edt check

	}// END: setStatus

	public void enableTreeFileButton() {

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {

				treePanel.enableTreeFileButton();

			}
		});
	}// END: enableTreeFileButton

	public void disableTreeFileButton() {

		if (SwingUtilities.isEventDispatchThread()) {

			treePanel.disableTreeFileButton();

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					treePanel.disableTreeFileButton();

				}
			});
		}// END: edt check

	}// END: disableTreeFileButton

	public void enableTreesFileButton() {

		if (SwingUtilities.isEventDispatchThread()) {

			treePanel.enableTreesFileButton();

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					treePanel.enableTreesFileButton();

				}
			});
		}// END: edt check

	}// END: enableTreesFileButton

	public void disableTreesFileButton() {

		if (SwingUtilities.isEventDispatchThread()) {

			treePanel.disableTreesFileButton();

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					treePanel.disableTreesFileButton();

				}
			});
		}// END: edt check

	}// END: disableTreesFileButton

	public void hideTreeColumn() {

		if (SwingUtilities.isEventDispatchThread()) {

			partitionsPanel.hideTreeColumn();

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					partitionsPanel.hideTreeColumn();

				}
			});
		}// END: edt check

	}// END: hideTreeColumn

	public void showTreeColumn() {

		if (SwingUtilities.isEventDispatchThread()) {

			partitionsPanel.showTreeColumn();

		} else {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {

					partitionsPanel.showTreeColumn();

				}
			});
		}// END: edt check

	}// END: showTreeColumn

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

}// END: class
