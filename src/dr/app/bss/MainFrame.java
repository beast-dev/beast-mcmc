/*
 * MainFrame.java
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

import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
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
import dr.evolution.alignment.SimpleAlignment;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class MainFrame extends DocumentFrame implements FileMenuHandler {

    private PartitionDataList dataList;
private BeagleSequenceSimulator beagleSequenceSimulator;
    
    private final String TAXA_TAB_NAME = "Taxa";
    private final String TREES_TAB_NAME = "Data";
    private final String PARTITIONS_TAB_NAME = "Partitions";
    private final String SIMULATION_TAB_NAME = "Simulation";
    private final String TERMINAL_TAB_NAME = "Terminal";

    private JTabbedPane tabbedPane = new JTabbedPane();
    private TaxaPanel taxaPanel;
    private TreesPanel treesPanel;
    private PartitionsPanel partitionsPanel;
    private SimulationPanel simulationPanel;
    private TerminalPanel terminalPanel;

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

        setSize(new Dimension(1300, 600));
        setMinimumSize(new Dimension(260, 100));

        taxaPanel = new TaxaPanel(dataList);
        treesPanel = new TreesPanel(this, dataList);
        partitionsPanel = new PartitionsPanel(dataList);
        simulationPanel = new SimulationPanel(this, dataList);
        terminalPanel = new TerminalPanel();

        tabbedPane.addTab(TAXA_TAB_NAME, null, taxaPanel);
        tabbedPane.addTab(TREES_TAB_NAME, null, treesPanel);
        tabbedPane.addTab(PARTITIONS_TAB_NAME, null, partitionsPanel);
        tabbedPane.addTab(SIMULATION_TAB_NAME, null, simulationPanel);
        tabbedPane.addTab(TERMINAL_TAB_NAME, null, terminalPanel);

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

        tabbedPane.setSelectedComponent(treesPanel);

    }// END: initializeComponents

    // ////////////////
    // ---SIMULATE---//
    // ////////////////

    // file chooser
    public void doExport() {

        if (getTreesCount() == 0) {

            tabbedPane.setSelectedComponent(treesPanel);
            Utils.showDialog("Please load at least one tree file before generating alignment.");

        } else {

            JFileChooser chooser = new JFileChooser();
            chooser.setDialogTitle("Simulate...");
            chooser.setMultiSelectionEnabled(false);
            chooser.setCurrentDirectory(workingDirectory);

            int returnVal = chooser.showSaveDialog(Utils.getActiveFrame());
            if (returnVal == JFileChooser.APPROVE_OPTION) {

                File file = chooser.getSelectedFile();

                collectAllSettings();
                generateNumberOfSimulations(file);

                File tmpDir = chooser.getCurrentDirectory();
                if (tmpDir != null) {
                    workingDirectory = tmpDir;
                }

            }// END: approve check

        }// END: tree loaded check

    }// END: doExport

    // threading, UI, exceptions handling
    private void generateNumberOfSimulations(final File outFile) {

        setBusy();

        SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

            ArrayList<TreeModel> simulatedTreeModelList = new ArrayList<TreeModel>();

            // Executed in background thread
            public Void doInBackground() {

                try {

                    if (BeagleSequenceSimulatorApp.VERBOSE) {
                        Utils.printPartitionDataList(dataList);
                        System.out.println();
                    }

                    long startingSeed = dataList.startingSeed;
                    for (int i = 0; i < dataList.simulationsCount; i++) {

                        String fullPath = Utils.getMultipleWritePath(outFile,
                               dataList.outputFormat.toString().toLowerCase(), i);
                        PrintWriter writer = new PrintWriter(new FileWriter(
                                fullPath));

                        ArrayList<Partition> partitionsList = new ArrayList<Partition>();

                        for (PartitionData data : dataList) {

                            if (data.record == null) {

                            	writer.close();
                                throw new RuntimeException(
                                        "Set data in Partitions tab for "
                                                + (partitionsList.size() + 1)
                                                + " partition.");
                                
                            } else {

                                TreeModel treeModel = data.createTreeModel();
                                simulatedTreeModelList.add(treeModel);

                                // create partition
                                Partition partition = new Partition(
                                        treeModel, //
                                        data.createBranchModel(), //
                                        data.createSiteRateModel(), //
                                        data.createClockRateModel(), //
                                        data.createFrequencyModel(), //
                                        data.from - 1, // from
                                        data.to - 1, // to
                                        data.every // every
                                );

                                if (data.ancestralSequenceString != null) {
                                    partition.setRootSequence(data.createAncestralSequence());
                                }

                                partitionsList.add(partition);

                            }

                        }// END: data list loop

                        if (dataList.setSeed) {
                            MathUtils.setSeed(startingSeed);
                            startingSeed += 1;
                        }

                         beagleSequenceSimulator = new BeagleSequenceSimulator(
                                partitionsList
                        );

                        SimpleAlignment alignment = beagleSequenceSimulator.simulate(dataList.useParallel, dataList.outputAncestralSequences);
                        alignment.setOutputType(dataList.outputFormat);
                        
//                        if (dataList.outputFormat == SimpleAlignment.OutputType.NEXUS) {
//                            alignment.setOutputType(dataList.outputFormat);
//                        } else if(dataList.outputFormat == SimpleAlignment.OutputType.XML) {
//                        	alignment.setOutputType(dataList.outputFormat);
//                        }else {
//                        	//
//                        }

                        writer.println(alignment.toString());
                        writer.close();

                    }// END: simulationsCount loop

                } catch (Exception e) {
                    Utils.handleException(e);
                    setStatus("Exception occured.");
                    setIdle();
                }

                return null;
            }// END: doInBackground

            // Executed in event dispatch thread
            public void done() {

//            	LinkedHashMap<Integer, LinkedHashMap<NodeRef, int[]>> partitionSequencesMap = beagleSequenceSimulator.getPartitionSequencesMap();
            	
                terminalPanel.setText(Utils.partitionDataListToString(dataList, simulatedTreeModelList
//                		, beagleSequenceSimulator.getPartitionSequencesMap()
                		));
                setStatus("Generated " + Utils.getSiteCount(dataList) + " sites.");
                setIdle();

            }// END: done
        };

        worker.execute();

    }// END: generateNumberOfSimulations

    // ////////////////////
    // ---GENERATE XML---//
    // ////////////////////

    public final void doGenerateXML() {

        if (getTreesCount() == 0) {

            tabbedPane.setSelectedComponent(treesPanel);
            Utils.showDialog("Please load at least one tree file before generating XML.");


        } else {

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

        }// END: tree loaded check

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
                    setStatus("Exception occured.");
                    setIdle();
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

    //	@Override
    public Action getSaveSettingsAction() {
        return new AbstractAction("Save settings...") {
            public void actionPerformed(ActionEvent ae) {
                doSaveSettings();
            }
        };
    }// END: generateXMLAction

    private void doSaveSettings() {

        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle("Save as...");
        chooser.setMultiSelectionEnabled(false);
        chooser.setCurrentDirectory(workingDirectory);

        int returnVal = chooser.showSaveDialog(Utils.getActiveFrame());
        if (returnVal == JFileChooser.APPROVE_OPTION) {

            File file = chooser.getSelectedFile();

            saveSettings(file);

            File tmpDir = chooser.getCurrentDirectory();
            if (tmpDir != null) {
                workingDirectory = tmpDir;
            }

            setStatus("Saved as " + file.getAbsolutePath());

        }// END: approve check

    }// END: saveSettings

    private void saveSettings(File file) {

        try {

            String fullPath = Utils.getWritePath(file, "bss");
            OutputStream fileOut = new FileOutputStream(new File(fullPath));
            ObjectOutputStream out = new ObjectOutputStream(fileOut);

            out.writeObject(dataList);
            out.close();
            fileOut.close();

        } catch (FileNotFoundException e) {
            Utils.handleException(e);
        } catch (IOException e) {
            Utils.handleException(e);
        }

    }// END: saveSettings

    //	@Override
    public Action getLoadSettingsAction() {
        return new AbstractAction("Load settings...") {
            public void actionPerformed(ActionEvent ae) {
                doLoadSettings();
            }
        };
    }// END: generateXMLAction

    private void doLoadSettings() {

        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle("Load...");
        chooser.setMultiSelectionEnabled(false);
        chooser.setCurrentDirectory(workingDirectory);

        int returnVal = chooser.showOpenDialog(Utils.getActiveFrame());
        if (returnVal == JFileChooser.APPROVE_OPTION) {

            File file = chooser.getSelectedFile();

            loadSettings(file);

            File tmpDir = chooser.getCurrentDirectory();
            if (tmpDir != null) {
                workingDirectory = tmpDir;
            }

            setStatus("Loaded " + file.getAbsolutePath());

        }// END: approve check

    }// END: doLoadSettings

    private void loadSettings(File file) {

        try {

            FileInputStream fileIn = new FileInputStream(file);
            ObjectInputStream in = new ObjectInputStream(fileIn);

            dataList = (PartitionDataList) in.readObject();

            in.close();
            fileIn.close();

            if (BeagleSequenceSimulatorApp.VERBOSE) {
                Utils.printPartitionDataList(dataList);
                System.out.println();
            }

			partitionsPanel.updatePartitionTable(dataList);
			taxaPanel.updateTaxaTable(dataList);
			treesPanel.updateTreesTable(dataList);
			simulationPanel.updateSimulationPanel(dataList);
            
        } catch (IOException ioe) {

            Utils.handleException(
                    ioe,
                    "Unable to read BSS file. "
                            + "BSS can only read files created by 'Saving' within BSS. "
                            + "It cannot read XML files.");

        } catch (ClassNotFoundException cnfe) {

            Utils.handleException(cnfe);

        }// END: try-catch block

    }// END: loadSettings

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

            doTaxaChanged();

        } else {

            SwingUtilities.invokeLater(new Runnable() {
                public void run() {

                    doTaxaChanged();

                }
            });
        }// END: edt check

    }// END: fireTaxaChanged

    private void doTaxaChanged() {

        treesPanel.fireTableDataChanged();
        taxaPanel.fireTaxaChanged();
        setStatus(Integer.toString(dataList.allTaxa
                .getTaxonCount()) + " taxa loaded.");

    }//END: doTaxaChanged

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

    public void disableTaxaPanel() {

        final int index = Utils.getTabbedPaneComponentIndex(tabbedPane, TAXA_TAB_NAME);

        if (SwingUtilities.isEventDispatchThread()) {

            tabbedPane.setEnabledAt(index, false);

        } else {

            SwingUtilities.invokeLater(new Runnable() {
                public void run() {

                    tabbedPane.setEnabledAt(index, false);

                }
            });
        }// END: edt check

    }// END: disableTaxaPanel

    public void enableTaxaPanel() {

        final int index = Utils.getTabbedPaneComponentIndex(tabbedPane, TAXA_TAB_NAME);

        if (SwingUtilities.isEventDispatchThread()) {

            tabbedPane.setEnabledAt(index, true);

        } else {

            SwingUtilities.invokeLater(new Runnable() {
                public void run() {

                    tabbedPane.setEnabledAt(index, true);

                }
            });
        }// END: edt check

    }// END: enableTaxaPanel

    private int getTreesCount() {
        return dataList.recordsList.size();
    }// END: getTreesCount

    //	@Override  // Please use Java 1.5
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
