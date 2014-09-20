/*
 * SkyGridDialog.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.tracer.analysis;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import dr.app.gui.util.LongTask;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceFactory;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import jam.framework.DocumentFrame;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SkyGridDialog {

    private JFrame frame;

    private String[][] argumentGuesses = {
            {"popsize"},
            {"cutoff"},
    };

    private String[] argumentNames = new String[]{
            "Population Size",
            "Grid Height",
    };

    //    private final JButton button = new JButton("Choose File...");
    private ActionListener buttonListener;

    private final JTextField fileNameText = new JTextField("not selected", 16);
    private File treeFile = null;
    private WholeNumberField binCountField;

    private String[] argumentTraces = new String[argumentNames.length];
    private JComboBox[] argumentCombos = new JComboBox[argumentNames.length];
    private JComboBox maxHeightCombo = new JComboBox(new String[]{
            "Lower 95% HPD", "Median", "Mean", "Upper 95% HPD"});
    private JComboBox rootHeightCombo;
    private JCheckBox manualRangeCheckBox;
    private RealNumberField minTimeField;
    private RealNumberField maxTimeField;
    private String rootHeightTrace = "None selected";

    private RealNumberField ageOfYoungestField = new RealNumberField();

    private OptionsPanel optionPanel;

    public SkyGridDialog(JFrame frame) {
        this.frame = frame;

        for (int i = 0; i < argumentNames.length; i++) {
            argumentCombos[i] = new JComboBox();
            argumentTraces[i] = "None selected";
        }

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

        maxHeightCombo.setSelectedIndex(3);
    }

    private int findArgument(JComboBox comboBox, String argument) {
        for (int i = 0; i < comboBox.getItemCount(); i++) {
            String item = ((String) comboBox.getItemAt(i)).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    private String getNumericalSuffix(String argument) {
        int i = argument.length() - 1;

        if (i < 0) return "";

        char ch = argument.charAt(i);

        if (!Character.isDigit(ch)) return "";

        while (i > 0 && Character.isDigit(ch)) {
            i -= 1;
            ch = argument.charAt(i);
        }

        return argument.substring(i + 1, argument.length());
    }

    private int getTraceRange(TraceList traceList, int first) {
        int i = 1;
        int k = first;

        String name = traceList.getTraceName(first);
        String root = name.substring(0, name.length() - 1);
        while (k < traceList.getTraceCount() && traceList.getTraceName(k).equals(root + i)) {
            i++;
            k++;
        }

        return i - 1;
    }

    public int showDialog(TraceList traceList, TemporalAnalysisFrame temporalAnalysisFrame) {

        Set<String> roots = new HashSet<String>();
        for (int j = 0; j < traceList.getTraceCount(); j++) {
            String statistic = traceList.getTraceName(j);
            String suffix = getNumericalSuffix(statistic);
            if (suffix.equals("1")) {
                roots.add(statistic.substring(0, statistic.length() - 1));
            } else if (statistic.equals("skygrid.cutOff")) {
                roots.add(statistic);
            }
        }

        if (roots.size() == 0) {
            JOptionPane.showMessageDialog(frame, "No traces found with a range of numerical suffixes (1-n).",
                    "Probably not a Sky Grid analysis",
                    JOptionPane.ERROR_MESSAGE);
            return JOptionPane.CANCEL_OPTION;
        }

        for (int i = 0; i < argumentCombos.length; i++) {
            argumentCombos[i].removeAllItems();

            for (String root : roots) {
                argumentCombos[i].addItem(root);
            }

            int index = findArgument(argumentCombos[i], argumentTraces[i]);

            for (int j = 0; j < argumentGuesses[i].length; j++) {
                if (index != -1) break;

                index = findArgument(argumentCombos[i], argumentGuesses[i][j]);
            }
            if (index == -1) index = 0;

            argumentCombos[i].setSelectedIndex(index);
        }

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

//        button.removeActionListener(buttonListener);
//        buttonListener = new ActionListener() {
//            public void actionPerformed(ActionEvent ae) {
//                FileDialog dialog = new FileDialog(frame,
//                        "Open Trees Log File...",
//                        FileDialog.LOAD);
//
//                dialog.setVisible(true);
//                if (dialog.getFile() == null) {
//                    // the dialog was cancelled...
//                    return;
//                }
//
//                treeFile = new File(dialog.getDirectory(), dialog.getFile());
//                fileNameText.setText(treeFile.getName());
//            }
//        };
//
//        button.addActionListener(buttonListener);

        final JDialog dialog = optionPane.createDialog(frame, "Sky Grid Analysis");
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
                for (int i = 0; i < argumentCombos.length; i++) {
                    argumentTraces[i] = argumentCombos[i].getSelectedItem() + "1";
                }
                rootHeightTrace = (String) rootHeightCombo.getSelectedItem();
            }
        } while (!done);

        return result;
    }

    private void setArguments(TemporalAnalysisFrame temporalAnalysisFrame) {
        optionPanel.removeAll();

        JLabel label = new JLabel(
                "<html>Warning! This analysis should only be run on traces where<br>" +
                        "the Sky Grid model was specified as the demographic in BEAST.<br>" +
                        "<em>Any other model will produce meaningless results.</em></html>");
        label.setFont(label.getFont().deriveFont(((float) label.getFont().getSize() - 2)));
        optionPanel.addSpanningComponent(label);
        optionPanel.addSeparator();

        if (treeFile != null) {
            fileNameText.setText(treeFile.getName());
        }
        fileNameText.setEditable(false);

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.add(fileNameText, BorderLayout.CENTER);
//        panel.add(button, BorderLayout.EAST);
//        optionPanel.addComponentWithLabel("Trees Log File: ", panel);
//
//        optionPanel.addSeparator();

        optionPanel.addLabel("Select the traces to use for the arguments:");

        for (int i = 0; i < argumentNames.length; i++) {
            optionPanel.addComponentWithLabel(argumentNames[i] + ":",
                    argumentCombos[i]);
        }

        optionPanel.addSeparator();

        optionPanel.addComponentWithLabel("Maximum time is the root height's:", maxHeightCombo);

        optionPanel.addComponentWithLabel("Select the trace of the root height:", rootHeightCombo);

        if (temporalAnalysisFrame == null) {
//            optionPanel.addSeparator();
//            optionPanel.addComponentWithLabel("Number of bins:", binCountField);
//
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

    public void createSkyGridFrame(TraceList traceList, DocumentFrame parent) {

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
                        "Error creating Bayesian SkyGrid",
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

        int firstPopSize = traceList.getTraceIndex(argumentTraces[0]);
        int popSizeCount = getTraceRange(traceList, firstPopSize);

        double gridHeight = (Double) traceList.getTrace(traceList.getTraceIndex(
                (String) argumentCombos[1].getSelectedItem())).getValue(0);

        final AnalyseSkyGridTask analyseTask = new AnalyseSkyGridTask(traceList,
                treeFile,
                firstPopSize,
                popSizeCount,
                gridHeight,
                frame);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Analysing Sky Grid",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

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

    class AnalyseSkyGridTask extends LongTask {

        TraceList traceList;
        TemporalAnalysisFrame frame;
        File treeFile;
        int firstPopSize;
        int popSizeCount;
        int binCount;
        boolean rangeSet;
        double minTime;
        double maxTime;
        double ageOfYoungest;
        double gridHeight;

        int stateCount;

        ArrayList<ArrayList> popSizes;

        private int lengthOfTask = 0;
        private int current = 0;

        public AnalyseSkyGridTask(TraceList traceList, File treeFile, int firstPopSize, int popSizeCount,
                                  double gridHeight, TemporalAnalysisFrame frame) {
            this.traceList = traceList;
            this.frame = frame;
            this.treeFile = treeFile;
            this.firstPopSize = firstPopSize;
            this.popSizeCount = popSizeCount;
            this.gridHeight = gridHeight;

            this.binCount = frame.getBinCount();
            this.rangeSet = frame.isRangeSet();

            ageOfYoungest = ageOfYoungestField.getValue();

            lengthOfTask = traceList.getStateCount() /* + binCount*/;

            stateCount = traceList.getStateCount();
        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Calculating SkyGrid...";
        }

        public String getMessage() {
            return null;
        }


        private double transform(double x) {
            return Math.exp(x);
        }

        public Object doWork() {

            popSizes = new ArrayList<ArrayList>();
            for (int i = 0; i < popSizeCount; i++) {
                ArrayList column = new ArrayList(traceList.getValues(firstPopSize + i));
                popSizes.add(column);
            }

            List heights = traceList.getValues(traceList.getTraceIndex(rootHeightTrace));

            TraceDistribution distribution = new TraceDistribution(heights,
                    traceList.getTrace(traceList.getTraceIndex(rootHeightTrace)).getTraceType(), traceList.getStepSize());

            double timeMean = distribution.getMean();
            double timeMedian = distribution.getMedian();
            double timeUpper = distribution.getUpperHPD();
            double timeLower = distribution.getLowerHPD();
            double timeGrid = gridHeight;

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

            // Check for large difference between maxHeight and cutOff
            double factor = 10.0;
            boolean condition = (gridHeight > maxHeight * factor || gridHeight < maxHeight / factor);
            if (condition) {

                String pt1 = "<html><body width='";
                String pt2 = "'><h1>Bayesian SkyGrid: Caution</h1>" +
                                "<p>Inferred root height is considerably smaller or larger than the SkyGrid cut-off value (" +
                                gridHeight +
                                ").  " +
                                "For improved interpretability, it is advisable to re-run the SkyGrid posterior inference using a cut-off commensurate with the root height.";

                int width = 400;
                String s = pt1 + width + pt2;

                JOptionPane.showMessageDialog(frame,
                        s, "Reconstruction Warning",
                        JOptionPane.OK_OPTION);
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
                timeGrid = ageOfYoungest - timeGrid;


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
                timeGrid = timeGrid - ageOfYoungest;

                // setting a timeXXXX to -1 means that it won't be displayed...
                if (maxTime <= timeLower) timeLower = -1;
                if (maxTime <= timeMean) timeMean = -1;
                if (maxTime <= timeMedian) timeMedian = -1;
                if (maxTime <= timeUpper) timeUpper = -1;
            }

            int numGridPoints = popSizeCount - 1;
            double delta = gridHeight / (double) numGridPoints;

            try {
                double height;
                if (ageOfYoungest > 0.0) {
                    height = ageOfYoungest;
                    delta = -delta;
                } else {
                    height = 0.0;
                }

                Variate.D xData = new Variate.D();
                Variate.D yDataMean = new Variate.D();
                Variate.D yDataMedian = new Variate.D();
                Variate.D yDataUpper = new Variate.D();
                Variate.D yDataLower = new Variate.D();

                double plotMin = minTime;
                double plotMax = maxTime;

                for (int i = 0; i < popSizeCount; ++i) {

                    if (height >= plotMin && height <= plotMax) {

                        xData.add(height);
                        TraceDistribution dist = new TraceDistribution(popSizes.get(i), TraceFactory.TraceType.DOUBLE);
                        yDataMean.add(transform(dist.getMean()));
                        yDataMedian.add(transform(dist.getMedian()));
                        yDataUpper.add(transform(dist.getUpperHPD()));
                        yDataLower.add(transform(dist.getLowerHPD()));

                        if (i == popSizeCount - 1) {
//                        double fillTime = (popSizeCount + 1.5) * gridSpacing;
                            while (height >= plotMin && height <= plotMax) {
                                xData.add(height);
                                yDataMean.add(transform(dist.getMean()));
                                yDataMedian.add(transform(dist.getMedian()));
                                yDataUpper.add(transform(dist.getUpperHPD()));
                                yDataLower.add(transform(dist.getLowerHPD()));
                                height += delta;
                            }
                        }
                    }
                    height += delta;
                }

                frame.addDemographic("Bayesian SkyGrid: " + traceList.getName(), xData,
                        yDataMean, yDataMedian,
                        yDataUpper, yDataLower,
                        timeMean, timeMedian,
                        timeUpper, timeLower);

            } catch (IllegalArgumentException ile) {
                JOptionPane.showMessageDialog(frame, ile.getMessage(),
                        "Invalid log file",
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

    private boolean inPlot(double height, double plotMin, double plotMax, double ageOfYoungest) {
        if (ageOfYoungest > 0.0) {
            return (height <= plotMin) && (height >= plotMax);
        } else {
            return (height >= plotMin) && (height <= plotMax);
        }
    }
}