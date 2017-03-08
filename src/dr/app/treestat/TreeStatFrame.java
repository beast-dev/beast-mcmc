/*
 * TreeStatFrame.java
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

package dr.app.treestat;

import dr.evolution.tree.TreeUtils;
import jam.framework.Application;
import jam.framework.DocumentFrame;
import jam.util.IconUtils;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import java.awt.*;
import java.io.*;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.app.treestat.statistics.TreeSummaryStatistic;

public class TreeStatFrame extends DocumentFrame {

    /**
	 *
	 */
	private static final long serialVersionUID = -1775448072034877658L;

	private TreeStatData treeStatData = null;

    private JTabbedPane tabbedPane = new JTabbedPane();

    private TaxonSetsPanel taxonSetsPanel;
    private StatisticsPanel statisticsPanel;

    private JLabel statusLabel;
    private JLabel progressLabel;
    private JProgressBar progressBar;

    final Icon gearIcon = IconUtils.getIcon(this.getClass(), "images/gear.png");

    public TreeStatFrame(Application application, String title) {
        super();

        setTitle(title);

        treeStatData = new TreeStatData();

        setImportAction(importTaxaAction);
        setExportAction(processTreeFileAction);

        getOpenAction().setEnabled(false);
        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);
    }

    public void initializeComponents() {

        setSize(new java.awt.Dimension(800, 600));

        taxonSetsPanel = new TaxonSetsPanel(this, treeStatData);
        statisticsPanel = new StatisticsPanel(this, treeStatData);

        tabbedPane.addTab("Statistics", null, statisticsPanel);
        tabbedPane.addTab("Taxon Sets", null, taxonSetsPanel);

        statusLabel = new JLabel("No statistics selected");
        processTreeFileAction.setEnabled(false);

        JPanel progressPanel = new JPanel(new BorderLayout(0,0));
        progressLabel = new JLabel("");
        progressBar = new JProgressBar();
        progressPanel.add(progressLabel, BorderLayout.NORTH);
        progressPanel.add(progressBar, BorderLayout.CENTER);

        JPanel panel2 = new JPanel(new FlowLayout());
        JButton goButton = new JButton(processTreeFileAction);
        goButton.setFocusable(false);
        goButton.putClientProperty("JButton.buttonType", "textured");
        goButton.setMargin(new Insets(4,4,4,4));
        panel2.add(goButton);
        panel2.add(progressPanel);

        JPanel panel1 = new JPanel(new BorderLayout(0,0));
        panel1.add(statusLabel, BorderLayout.WEST);
        panel1.add(panel2, BorderLayout.EAST);
        panel1.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(0, 6, 0, 6)));

        JPanel panel = new JPanel(new BorderLayout(0,0));
        panel.add(tabbedPane, BorderLayout.CENTER);
        panel.add(panel1, BorderLayout.SOUTH);
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        getContentPane().setLayout(new java.awt.BorderLayout(0, 0));
        getContentPane().add(panel, BorderLayout.CENTER);
    }

    public void fireDataChanged() {
        if (treeStatData.statistics.size() > 0) {
            statusLabel.setText("" + treeStatData.statistics.size() + " statistics selected");
            processTreeFileAction.setEnabled(true);
        } else {
            statusLabel = new JLabel("No statistics selected");
            processTreeFileAction.setEnabled(false);
        }
        taxonSetsPanel.dataChanged();
        statisticsPanel.dataChanged();

    }

    protected boolean readFromFile(File file) throws IOException {
        return false;
    }

    protected boolean writeToFile(File file) {
        return false;
    }

    public final void doImport() {

        FileDialog dialog = new FileDialog(this,
                                       "Import Tree File...",
                                       FileDialog.LOAD);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                importFromFile(file);

            } catch (Importer.ImportException ie) {
                JOptionPane.showMessageDialog(this, "Unable to read tree file: " + ie,
                                                            "Unable to read tree file",
                                                            JOptionPane.ERROR_MESSAGE);
            } catch (FileNotFoundException fnfe) {
                JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                                                            "Unable to open file",
                                                            JOptionPane.ERROR_MESSAGE);
            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to read file: " + ioe,
                                                            "Unable to read file",
                                                            JOptionPane.ERROR_MESSAGE);
            } catch (Exception e) {
                JOptionPane.showMessageDialog(this, "Error: " + e,
                                                            "Error",
                                                            JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    protected void importFromFile(File file) throws IOException, Importer.ImportException {

        BufferedReader reader = new BufferedReader(new FileReader(file));

        String line = reader.readLine();
        Tree tree = null;

        if (line.toUpperCase().startsWith("#NEXUS")) {
            NexusImporter importer = new NexusImporter(reader);
            tree = importer.importTree(null);
        } else {
            reader.close();
            reader = new BufferedReader(new FileReader(file));
            NewickImporter importer = new NewickImporter(reader);
            tree = importer.importTree(null);
        }

        treeStatData.allTaxa = TreeUtils.getLeafSet(tree);
        statusLabel.setText(Integer.toString(treeStatData.allTaxa.size()) + " taxa loaded.");
        reader.close();

        fireDataChanged();
    }

    public final void doExport() {

        FileDialog inDialog = new FileDialog(this,
                                       "Import Tree File...",
                                       FileDialog.LOAD);

        inDialog.setVisible(true);
        if (inDialog.getFile() != null) {
            File inFile = new File(inDialog.getDirectory(), inDialog.getFile());


            FileDialog outDialog = new FileDialog(this,
                                           "Save Log File As...",
                                           FileDialog.SAVE);

            outDialog.setVisible(true);
            if (outDialog.getFile() != null) {

                File outFile = new File(outDialog.getDirectory(), outDialog.getFile());

                try {
                    processTreeFile(inFile, outFile);

                } catch (FileNotFoundException fnfe) {
                    JOptionPane.showMessageDialog(this, "Unable to open file: File not found",
                                                                "Unable to open file",
                                                                JOptionPane.ERROR_MESSAGE);
                } catch (IOException ioe) {
                    JOptionPane.showMessageDialog(this, "Unable to read/write file: " + ioe,
                                                                "Unable to read/write file",
                                                                JOptionPane.ERROR_MESSAGE);
                } catch (Importer.ImportException ie) {
                    JOptionPane.showMessageDialog(this, "Unable to import file: " + ie,
                                                                "Unable to import tree file",
                                                                JOptionPane.ERROR_MESSAGE);
                } catch (Exception e) {
                    JOptionPane.showMessageDialog(this, "Error: " + e,
                                                                "Error",
                                                                JOptionPane.ERROR_MESSAGE);
                }
            }
        }
    }

    protected void processTreeFile(File inFile, File outFile) throws IOException, Importer.ImportException {
        processTreeFileAction.setEnabled(false);

        BufferedReader r = new BufferedReader(new FileReader(inFile));
        String line = r.readLine();
        r.close();

        final ProgressMonitorInputStream in = new ProgressMonitorInputStream(
                this,
                "Reading " + inFile.getName(),
                new FileInputStream(inFile));
        in.getProgressMonitor().setMillisToDecideToPopup(0);
        in.getProgressMonitor().setMillisToPopup(0);

        final Reader reader = new InputStreamReader(new BufferedInputStream(in));

//        final Reader reader = new FileReader(inFile);
        final TreeImporter importer;

        if (line.toUpperCase().startsWith("#NEXUS")) {
            importer = new NexusImporter(reader);
        } else {
            reader.close();
            importer = new NewickImporter(reader);
        }

        final Tree firstTree = importer.importNextTree();
        boolean isUltrametric = TreeUtils.isUltrametric(firstTree);
        boolean isBinary = TreeUtils.isBinary(firstTree);
        boolean stop = false;

        // check that the trees conform with the requirements of the selected statistics
        for (int i = 0; i < treeStatData.statistics.size(); i++) {
            TreeSummaryStatistic tss = (TreeSummaryStatistic)treeStatData.statistics.get(i);
            String label = tss.getSummaryStatisticName();

            if (!isUltrametric && !tss.allowsNonultrametricTrees()) {
                if (JOptionPane.showConfirmDialog(
                        this, "Warning: These trees may not be ultrametric and this is\na requirement of the " +
                        label + " statistic. Do you wish to continue?", "Warning", JOptionPane.YES_NO_OPTION) != JOptionPane.YES_OPTION) {
                    stop = true;
                    break;
                }
                // don't ask the question again...
                isUltrametric = true;
            }

            if (!isBinary && !tss.allowsPolytomies()) {
                if (JOptionPane.showConfirmDialog(
                        this, "Warning: These trees may not be strictly bifurcating and this is\na requirement of the " +
                        label + " statistic. Do you wish to continue?", "Warning", JOptionPane.YES_NO_OPTION) != JOptionPane.YES_OPTION) {
                    stop = true;
                    break;
                }
                // don't ask the question again...
                isBinary = true;
            }
        }

        if (stop) {
            processTreeFileAction.setEnabled(true);
            return;
        }

        final PrintWriter writer = new PrintWriter(new FileWriter(outFile));

//        Thread readThread = new Thread() {
//            public void run() {
                Tree tree = firstTree;

                writer.print("state");
                for (int i = 0; i < treeStatData.statistics.size(); i++) {
                    TreeSummaryStatistic tss = (TreeSummaryStatistic)treeStatData.statistics.get(i);

                    int dim = tss.getStatisticDimensions(tree);
                    for (int j = 0; j < dim; j++) {
                        writer.print("\t" + tss.getStatisticLabel(tree, j));
                    }

                }
                writer.println();

                state = 0;
                do {
                    writer.print(state);

                    for (int i = 0; i < treeStatData.statistics.size(); i++) {
                        TreeSummaryStatistic tss = (TreeSummaryStatistic)treeStatData.statistics.get(i);
                        double[] stats = tss.getSummaryStatistic(tree);
                        for (int j = 0; j < stats.length; j++) {
                            writer.print("\t" + stats[j]);
                        }
                    }
                    writer.println();
                    state += 1;

                    final int currentState = state;

                    in.getProgressMonitor().setNote("Processing Tree " + currentState + "...");
//                    EventQueue.invokeLater(
//                            new Runnable() {
//                                public void run() {
//                                    progressLabel.setText("Processing Tree " + currentState + "...");
//                                }
//                            });


//                    try {
                        tree = importer.importNextTree();
//                    } catch (final IOException e) {
//                        EventQueue.invokeLater(
//                                new Runnable() {
//                                    public void run() {
//                                        JOptionPane.showMessageDialog(TreeStatFrame.this, "File I/O Error: " + e.getMessage(),
//                                                "File I/O Error",
//                                                JOptionPane.ERROR_MESSAGE);
//                                    }
//                                });
//                    } catch (final Importer.ImportException e) {
//                        EventQueue.invokeLater(
//                                new Runnable() {
//                                    public void run() {
//                                        JOptionPane.showMessageDialog(TreeStatFrame.this, "Error importing tree: " + e.getMessage(),
//                                                "Tree Import Error",
//                                                JOptionPane.ERROR_MESSAGE);
//                                    }
//                                });
//                    }
                } while (tree != null);
//            }
//        };
//
//        readThread.start();
//        while (readThread.isAlive()) {
//            Thread.yield();
//        }

        reader.close();
        writer.close();

        progressLabel.setText("" + state + " trees processed.");
        processTreeFileAction.setEnabled(true);
    }

    private int state = 0;

    public void doCopy() {
//		statisticsPanel.doCopy();
    }

    public JComponent getExportableComponent() {

        return statisticsPanel.getExportableComponent();
    }

      protected AbstractAction importTaxaAction = new AbstractAction("Import Taxa...") {
          /**
		 *
		 */
		private static final long serialVersionUID = -3185667996732228702L;

		public void actionPerformed(java.awt.event.ActionEvent ae) {
              doImport();
          }
      };

      protected AbstractAction processTreeFileAction = new AbstractAction("Process Tree File...", gearIcon) {
          /**
		 *
		 */
		private static final long serialVersionUID = -8285433136692586532L;

		public void actionPerformed(java.awt.event.ActionEvent ae) {
              doExport();
          }
      };

    TreeImporter treeImporter;
}
