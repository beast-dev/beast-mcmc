/*
 * CoalGenFrame.java
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

package dr.app.coalgen;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import jam.framework.DocumentFrame;
import jam.framework.Exportable;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.io.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class CoalGenFrame extends DocumentFrame {

    private CoalGenData data = null;

    private final JTabbedPane tabbedPane = new JTabbedPane();
    private InputsPanel inputsPanel;
    private DataPanel dataPanel;
    private ModelPanel modelPanel;
    private SimulationsPanel simulationsPanel;

    private JLabel statusLabel;

    public CoalGenFrame(String title) {
        super();

        setTitle(title);

        data = new CoalGenData();
    }

    public void initializeComponents() {

        setSize(new java.awt.Dimension(800, 600));

        inputsPanel = new InputsPanel(this, data);
        dataPanel = new DataPanel(this, data);
        modelPanel = new ModelPanel(this, data);
        simulationsPanel = new SimulationsPanel(this, data);

        tabbedPane.addTab("Inputs", null, inputsPanel);
        tabbedPane.addTab("Taxa", null, dataPanel);
        tabbedPane.addTab("Model", null, modelPanel);
        tabbedPane.addTab("Simulations", null, simulationsPanel);

        statusLabel = new JLabel("No taxa loaded");

        JPanel progressPanel = new JPanel(new BorderLayout(0, 0));
        JLabel progressLabel = new JLabel("");
        JProgressBar progressBar = new JProgressBar();
        progressPanel.add(progressLabel, BorderLayout.NORTH);
        progressPanel.add(progressBar, BorderLayout.CENTER);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(statusLabel, BorderLayout.CENTER);
        panel1.add(progressPanel, BorderLayout.EAST);
        panel1.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 6, 0, 6)));

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.add(tabbedPane, BorderLayout.CENTER);
        panel.add(panel1, BorderLayout.SOUTH);
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);
    }

    public void fireTaxaChanged() {
    }

    public void fireModelChanged() {
        modelPanel.collectSettings();
    }


    public void fireTracesChanged() {
        inputsPanel.tracesChanged();
        modelPanel.tracesChanged();
        simulationsPanel.tracesChanged();
    }

    public final void dataSelectionChanged(boolean isSelected) {
        if (isSelected) {
            getDeleteAction().setEnabled(true);
        } else {
            getDeleteAction().setEnabled(false);
        }
    }

    public void doDelete() {
//        if (tabbedPane.getSelectedComponent() == dataPanel) {
//            dataPanel.deleteSelection();
//        } else {
//            throw new RuntimeException("Delete should only be accessable from the Data panel");
//        }
    }

    protected boolean writeToFile(File file) {
        inputsPanel.collectSettings();
        modelPanel.collectSettings();
        simulationsPanel.collectSettings();
        return false;
    }

    protected boolean readFromFile(final File file) throws IOException {

        try {
            final String fileName = file.getName();
            final ProgressMonitorInputStream in = new ProgressMonitorInputStream(
                    this,
                    "Reading " + fileName,
                    new FileInputStream(file));

//            final Reader reader = new InputStreamReader(in);
            final JFrame frame = this;

            // the monitored activity must be in a new thread.
            Thread readThread = new Thread() {
                public void run() {
                    try {
                        final File file1 = new File(fileName);
                        final LogFileTraces traces = new LogFileTraces(fileName, file1);
                        traces.loadTraces(in);

                        EventQueue.invokeLater(
                                new Runnable() {
                                    public void run() {
                                        data.logFile = file;
                                        data.traces = traces;
                                        fireTracesChanged();
                                    }
                                });

                    } catch (final TraceException tex) {
                        EventQueue.invokeLater(
                                new Runnable() {
                                    public void run() {
                                        JOptionPane.showMessageDialog(frame, "Error reading trace file: " + tex,
                                                "Error reading trace file",
                                                JOptionPane.ERROR_MESSAGE);
                                    }
                                });
                    } catch (final InterruptedIOException iioex) {
                        // The cancel dialog button was pressed - do nothing
                    } catch (final IOException ioex) {
                        EventQueue.invokeLater(
                                new Runnable() {
                                    public void run() {
                                        JOptionPane.showMessageDialog(frame, "File I/O Error: " + ioex,
                                                "File I/O Error",
                                                JOptionPane.ERROR_MESSAGE);
                                    }
                                });
                    } catch (final Exception ex) {
                        EventQueue.invokeLater(
                                new Runnable() {
                                    public void run() {
                                        JOptionPane.showMessageDialog(frame, "Fatal exception: " + ex,
                                                "Error reading file",
                                                JOptionPane.ERROR_MESSAGE);
                                    }
                                });
                    }

                }
            };
            readThread.start();

        } catch (IOException ioex) {
            JOptionPane.showMessageDialog(this, "File I/O Error: " + ioex,
                    "File I/O Error",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Fatal exception: " + ex,
                    "Error reading file",
                    JOptionPane.ERROR_MESSAGE);
            return false;
        }
        return true;
    }

    public final void doImport() {

        FileDialog dialog = new FileDialog(this,
                "Import Tree or Alignment...",
                FileDialog.LOAD);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                importFromFile(file);


            } catch (Importer.ImportException ie) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ie,
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                        "Unable to open file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe,
                        "Unable to read file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    protected void importFromFile(File file) throws IOException, Importer.ImportException {

        BufferedReader reader = new BufferedReader(new FileReader(file));

        String line = reader.readLine();
        Tree tree;

        if (line.toUpperCase().startsWith("#NEXUS")) {
            NexusImporter importer = new NexusImporter(reader);
            tree = importer.importTree(null);
        } else {
            NewickImporter importer = new NewickImporter(reader);
            tree = importer.importTree(null);
        }

        data.taxonList = tree;
        statusLabel.setText(Integer.toString(data.taxonList.getTaxonCount()) + " taxa loaded.");
        reader.close();

        fireTaxaChanged();
    }

    public final void doExport() {

        FileDialog outDialog = new FileDialog(this,
                "Save Log File As...",
                FileDialog.SAVE);

        outDialog.setVisible(true);
        if (outDialog.getFile() != null) {

            File outFile = new File(outDialog.getDirectory(), outDialog.getFile());

            try {
                generateFile(outFile);

            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to write file: " + ioe,
                        "Unable to read/write file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    protected void generateFile(File outFile) throws IOException {

        PrintWriter writer = new PrintWriter(new FileWriter(outFile));

        dr.evolution.coalescent.CoalescentSimulator simulator = new dr.evolution.coalescent.CoalescentSimulator();

        int count = 0;
        while (data.hasNext()) {
            DemographicFunction demo = data.nextDemographic();

            Tree tree = simulator.simulateTree(data.taxonList, demo);

            writer.println(count + "\t" + TreeUtils.newick(tree));
            count += 1;
        }

        writer.close();
    }

    public void doCopy() {
    }

    public JComponent getExportableComponent() {

        JComponent exportable = null;
        Component comp = tabbedPane.getSelectedComponent();

        if (comp instanceof Exportable) {
            exportable = ((Exportable) comp).getExportableComponent();
        } else if (comp instanceof JComponent) {
            exportable = (JComponent) comp;
        }

        return exportable;
    }

    public boolean useImportAction() {
        return true;
    }

    public Action getImportAction() {
        return importTaxaAction;
    }

    protected AbstractAction importTaxaAction = new AbstractAction("Import Taxa...") {
        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doImport();
        }
    };

    public boolean useExportAction() {
        return true;
    }

    public Action getExportAction() {
        return generateSimulationsAction;
    }

    protected AbstractAction generateSimulationsAction = new AbstractAction("Generate Simulations...") {
        public void actionPerformed(java.awt.event.ActionEvent ae) {
            doExport();
        }
    };
}
