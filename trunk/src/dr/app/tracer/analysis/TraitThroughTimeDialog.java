/*
 * TraitThroughTimeDialog.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tracer.analysis;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import dr.app.gui.util.LongTask;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import jam.framework.DocumentFrame;
import jam.panels.OptionsPanel;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;

public class TraitThroughTimeDialog {

    private JFrame frame;

    private final JButton button = new JButton("Choose File...");
    private ActionListener buttonListener;

    private final JTextField fileNameText = new JTextField("not selected", 16);
    private File treeFile = null;
    private String traitName = null;
    private WholeNumberField binCountField;

    private JTextField traitNameField;
    private JComboBox maxHeightCombo = new JComboBox(new String[]{
            "Lower 95% HPD", "Median", "Mean", "Upper 95% HPD"});
    private JComboBox rootHeightCombo;
    private JCheckBox manualRangeCheckBox;
    private RealNumberField minTimeField;
    private RealNumberField maxTimeField;
    private String rootHeightTrace = "None selected";

    private RealNumberField ageOfYoungestField = new RealNumberField();

    private OptionsPanel optionPanel;

    public TraitThroughTimeDialog(JFrame frame) {
        this.frame = frame;

        traitNameField = new JTextField();
        traitNameField.setColumns(24);

        rootHeightCombo = new JComboBox();

        binCountField = new WholeNumberField(2, 2000);
        binCountField.setValue(100);
        binCountField.setColumns(4);

        manualRangeCheckBox = new JCheckBox("Use manual range for bins:");

        maxTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
        maxTimeField.setColumns(12);

        minTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
        minTimeField.setColumns(12);

        ageOfYoungestField.setValue(0.0);
        ageOfYoungestField.setColumns(12);

        optionPanel = new OptionsPanel(12, 12);
    }

    private int findArgument(JComboBox comboBox, String argument) {
        for (int i = 0; i < comboBox.getItemCount(); i++) {
            String item = ((String) comboBox.getItemAt(i)).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    public int showDialog(TraceList traceList, TemporalAnalysisFrame temporalAnalysisFrame) {

        setArguments(temporalAnalysisFrame);

        for (int j = 0; j < traceList.getTraceCount(); j++) {
            String statistic = traceList.getTraceName(j);
            rootHeightCombo.addItem(statistic);
        }
        int index = findArgument(rootHeightCombo, rootHeightTrace);
        if (index == -1) index = findArgument(rootHeightCombo, "root");
        if (index == -1) index = findArgument(rootHeightCombo, "height");
        if (index == -1) index = 0;
        rootHeightCombo.setSelectedIndex(index);

        final JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        button.removeActionListener(buttonListener);
        buttonListener = new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                FileDialog dialog = new FileDialog(frame,
                        "Open Trees Log File...",
                        FileDialog.LOAD);

                dialog.setVisible(true);
                if (dialog.getFile() == null) {
                    // the dialog was cancelled...
                    return;
                }

                treeFile = new File(dialog.getDirectory(), dialog.getFile());
                fileNameText.setText(treeFile.getName());
            }
        };

        button.addActionListener(buttonListener);

        final JDialog dialog = optionPane.createDialog(frame, "Trait Through Time Analysis");
        dialog.pack();

        int result = JOptionPane.CANCEL_OPTION;

        boolean done;
        do {
            done = true;
            dialog.setVisible(true);

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }

            if (result == JOptionPane.OK_OPTION) {
                if (treeFile == null) {
                    JOptionPane.showMessageDialog(frame, "A tree file was not selected",
                            "Error parsing file",
                            JOptionPane.ERROR_MESSAGE);
                    done = false;
                } else {
                    rootHeightTrace = (String) rootHeightCombo.getSelectedItem();
                }

                traitName = traitNameField.getText();
                if (traitName == null || traitName.length() == 0) {
                    JOptionPane.showMessageDialog(frame, "A trait name was not specified",
                            "Error parsing file",
                            JOptionPane.ERROR_MESSAGE);
                    done = false;
                }
            }
        } while (!done);

        return result;
    }

    private void setArguments(TemporalAnalysisFrame temporalAnalysisFrame) {
        optionPanel.removeAll();

        if (treeFile != null) {
            fileNameText.setText(treeFile.getName());
        }
        fileNameText.setEditable(false);

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.add(fileNameText, BorderLayout.CENTER);
        panel.add(button, BorderLayout.EAST);
        optionPanel.addComponentWithLabel("Trees Log File: ", panel);

        optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Name of Trait:", traitNameField);

        optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Maximum time is the root height's:", maxHeightCombo);

        optionPanel.addComponentWithLabel("Select the trace of the root height:", rootHeightCombo);

        if (temporalAnalysisFrame == null) {
            optionPanel.addSeparator();
            optionPanel.addComponentWithLabel("Number of bins:", binCountField);

            optionPanel.addSpanningComponent(manualRangeCheckBox);
            final JLabel label1 = optionPanel.addComponentWithLabel("Minimum time:", minTimeField);
            final JLabel label2 = optionPanel.addComponentWithLabel("Maximum time:", maxTimeField);

            if (manualRangeCheckBox.isSelected()) {
                label1.setEnabled(true);
                minTimeField.setEnabled(true);
                label2.setEnabled(true);
                maxTimeField.setEnabled(true);
            } else {
                label1.setEnabled(false);
                minTimeField.setEnabled(false);
                label2.setEnabled(false);
                maxTimeField.setEnabled(false);
            }

            manualRangeCheckBox.addChangeListener(new ChangeListener() {
                public void stateChanged(ChangeEvent changeEvent) {
                    if (manualRangeCheckBox.isSelected()) {
                        label1.setEnabled(true);
                        minTimeField.setEnabled(true);
                        label2.setEnabled(true);
                        maxTimeField.setEnabled(true);
                    } else {
                        label1.setEnabled(false);
                        minTimeField.setEnabled(false);
                        label2.setEnabled(false);
                        maxTimeField.setEnabled(false);
                    }
                }
            });
        }

        optionPanel.addSeparator();
        optionPanel.addComponentWithLabel("Age of youngest tip:", ageOfYoungestField);
        JLabel label3 = new JLabel(
                "<html>You can set the age of sampling of the most recent tip in<br>" +
                        "the tree. If this is set to zero then the plot is shown going<br>" +
                        "backwards in time, otherwise forwards in time.</html>");
        label3.setFont(label3.getFont().deriveFont(((float) label3.getFont().getSize() - 2)));
        optionPanel.addSpanningComponent(label3);
    }

    Timer timer = null;

    public void createTraitThroughTimeFrame(TraceList traceList, DocumentFrame parent) {

        TemporalAnalysisFrame frame;

        int binCount = binCountField.getValue();
        double minTime;
        double maxTime;
        boolean manualRange = manualRangeCheckBox.isSelected();
        if (manualRange) {
            minTime = minTimeField.getValue();
            maxTime = maxTimeField.getValue();

            if (minTime >= maxTime) {
                JOptionPane.showMessageDialog(parent,
                        "The minimum time value should be less than the maximum.",
                        "Error creating Trait-Through-Time plot",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            frame = new TemporalAnalysisFrame(parent, "", binCount, minTime, maxTime);
        } else {
            frame = new TemporalAnalysisFrame(parent, "", binCount);
        }
        frame.initialize();
        addToTemporalAnalysis(traceList, frame);
    }

    public void addToTemporalAnalysis(TraceList traceList, TemporalAnalysisFrame frame) {

        final AnalyseTTTTask analyseTask = new AnalyseTTTTask(traceList,
                treeFile,
                frame);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Analysing Trait-Through-Time",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);
        progressMonitor.setProgress(0);

        timer = new Timer(1000, new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                progressMonitor.setProgress(analyseTask.getCurrent());
                if (progressMonitor.isCanceled() || analyseTask.done()) {
                    progressMonitor.close();
                    analyseTask.stop();
                    timer.stop();
                }
            }
        });

        analyseTask.go();
        timer.start();
    }

    class AnalyseTTTTask extends LongTask {

        TraceList traceList;
        TemporalAnalysisFrame frame;
        File treeFile;
        int binCount;
        boolean rangeSet;
        double minTime;
        double maxTime;
        double ageOfYoungest;

        int stateCount;

        private int lengthOfTask = 0;
        private int current = 1;

        public AnalyseTTTTask(TraceList traceList, File treeFile, TemporalAnalysisFrame frame) {
            this.traceList = traceList;
            this.frame = frame;
            this.treeFile = treeFile;

            this.binCount = frame.getBinCount();
            this.rangeSet = frame.isRangeSet();

            ageOfYoungest = ageOfYoungestField.getValue();

            lengthOfTask = traceList.getStateCount();

            stateCount = traceList.getStateCount();

        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Calculating Traits...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {

            List heights = traceList.getValues(traceList.getTraceIndex(rootHeightTrace));

            TraceDistribution distribution = new TraceDistribution(heights,
                    traceList.getTrace(traceList.getTraceIndex(rootHeightTrace)).getTraceType(), traceList.getStepSize());

            double timeMean = distribution.getMean();
            double timeMedian = distribution.getMedian();
            double timeUpper = distribution.getUpperHPD();
            double timeLower = distribution.getLowerHPD();

            double maxHeight = 0.0;
            switch (maxHeightCombo.getSelectedIndex()) {
                // setting a timeXXXX to -1 means that it won't be displayed...
                case 0:
                    maxHeight = timeLower;
                    break;
                case 1:
                    maxHeight = timeMedian;
                    break;
                case 2:
                    maxHeight = timeMean;
                    break;
                case 3:
                    maxHeight = timeUpper;
                    break;
            }

            if (rangeSet) {
                minTime = frame.getMinTime();
                maxTime = frame.getMaxTime();
            } else {
                if (ageOfYoungest > 0.0) {
                    minTime = ageOfYoungest - maxHeight;
                    maxTime = ageOfYoungest;
                } else {
                    minTime = 0.0;
                    maxTime = maxHeight - ageOfYoungest;
                }

                frame.setRange(minTime, maxTime);
            }

            if (ageOfYoungest > 0.0) {
                // reverse them if ageOfYoungest is set positive
                timeMean = ageOfYoungest - timeMean;
                timeMedian = ageOfYoungest - timeMedian;
                timeUpper = ageOfYoungest - timeUpper;
                timeLower = ageOfYoungest - timeLower;

                // setting a timeXXXX to -1 means that it won't be displayed...
                if (minTime >= timeLower) timeLower = -1;
                if (minTime >= timeMean) timeMean = -1;
                if (minTime >= timeMedian) timeMedian = -1;
                if (minTime >= timeUpper) timeUpper = -1;
            } else {
                // otherwise use use ageOfYoungest as an offset
                timeMean = timeMean - ageOfYoungest;
                timeMedian = timeMedian - ageOfYoungest;
                timeUpper = timeUpper - ageOfYoungest;
                timeLower = timeLower - ageOfYoungest;

                // setting a timeXXXX to -1 means that it won't be displayed...
                if (maxTime <= timeLower) timeLower = -1;
                if (maxTime <= timeMean) timeMean = -1;
                if (maxTime <= timeMedian) timeMedian = -1;
                if (maxTime <= timeUpper) timeUpper = -1;
            }

            double delta = (maxTime - minTime) / (binCount - 1);

            try {
                BufferedReader reader = new BufferedReader(new FileReader(treeFile));

                String line = reader.readLine();

                TreeImporter importer;
                if (line.toUpperCase().startsWith("#NEXUS")) {
                    importer = new NexusImporter(reader);
                } else {
                    importer = new NewickImporter(reader, false);
                }

                int burnin = traceList.getBurnIn();
                int skip = burnin / traceList.getStepSize();
                int state = 0;

                while (importer.hasTree() && state < skip) {
                    importer.importNextTree();
                    state += 1;
                }

                current = 0;

                Variate.D[] bins = new Variate.D[binCount];
                for (int k = 0; k < binCount; k++) {
                    bins[k] = new Variate.D();
                }

                double startHeight;
                if (ageOfYoungest > 0.0) {
                    startHeight = ageOfYoungest - maxTime;
                } else {
                    startHeight = ageOfYoungest;
                }

                //int tips = 0;
                state = 0;

                try {
                    while (importer.hasTree()) {
                        boolean traitFound = false;
                        RootedTree tree = (RootedTree) importer.importNextTree();

                        for (Node node : tree.getNodes()) {
                            if (!tree.isRoot(node)) {
                                double h1 = tree.getHeight(node);
                                double h0 = tree.getHeight(tree.getParent(node));

                                Double value = (Double) node.getAttribute(traitName);
                                if (value != null) {
                                    traitFound = true;
                                    double height = startHeight;

                                    for (int k = 0; k < binCount; k++) {
                                        if (height >= h1) {
                                            if (height >= h0) {
                                                break;
                                            }
                                            bins[k].add(value);
                                        }
                                        height += delta;
                                    }
                                }


                            }
                        }

                        if (!traitFound) {
                            throw new ImportException("Tree missing values for trait, " + traitName);
                        }
                        state += 1;
                        current += 1;
                    }

                } catch (ImportException ie) {
                    JOptionPane.showMessageDialog(frame, "Error parsing file: " + ie.getMessage(),
                            "Error parsing file",
                            JOptionPane.ERROR_MESSAGE);
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(frame, "Fatal exception (email the authors):" + ex.getMessage(),
                            "Fatal exception",
                            JOptionPane.ERROR_MESSAGE);
                    ex.printStackTrace(System.out);
                }

                Variate.D xData = new Variate.D();
                Variate.D yDataMean = new Variate.D();
                Variate.D yDataMedian = new Variate.D();
                Variate.D yDataUpper = new Variate.D();
                Variate.D yDataLower = new Variate.D();

                double t;
                if (ageOfYoungest > 0.0) {
                    t = maxTime;
                } else {
                    t = minTime;
                }
                for (Variate.D bin : bins) {
                    xData.add(t);
                    if (bin.getCount() > 0) {
                        yDataMean.add(bin.getMean());
                        yDataMedian.add(bin.getQuantile(0.5));
                        yDataLower.add(bin.getQuantile(0.025));
                        yDataUpper.add(bin.getQuantile(0.975));
                    } else {
                        yDataMean.add(Double.NaN);
                        yDataMedian.add(Double.NaN);
                        yDataLower.add(Double.NaN);
                        yDataUpper.add(Double.NaN);
                    }
                    if (ageOfYoungest > 0.0) {
                        t -= delta;
                    } else {
                        t += delta;
                    }
                }

                frame.addDemographic("Trait Through Time: " + traceList.getName(), xData,
                        yDataMean, yDataMedian,
                        yDataUpper, yDataLower,
                        timeMean, timeMedian,
                        timeUpper, timeLower);

            } catch (java.io.IOException ioe) {
                JOptionPane.showMessageDialog(frame, "Error reading file: " + ioe.getMessage(),
                        "Error reading file",
                        JOptionPane.ERROR_MESSAGE);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(frame, "Fatal exception (email the authors):" + ex.getMessage(),
                        "Fatal exception",
                        JOptionPane.ERROR_MESSAGE);
                ex.printStackTrace(System.out);
            }

            return null;
        }

    }
}