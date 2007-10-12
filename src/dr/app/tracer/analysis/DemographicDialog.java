/*
 * DemographicDialog.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;
import dr.util.Variate;
import jebl.evolution.coalescent.*;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.util.LongTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

public class DemographicDialog {

    private JFrame frame;

    private JComboBox demographicCombo;
    private WholeNumberField binCountField;

    public static String[] demographicModels = {"Constant Population",
            "Exponential Growth (Growth Rate)",
            "Exponential Growth (Doubling Time)",
            "Logistic Growth (Growth Rate)",
            "Logistic Growth (Doubling Time)",
            "Expansion (Growth Rate)",
            "Expansion (Doubling Time)",
            "Boom-Bust"};

    private String[][] argumentGuesses = {
            {"populationsize", "population", "popsize", "n0", "size", "pop"},
            {"ancestralsize", "ancestralproportion", "proportion", "ancestral", "n1"},
            {"growthrate", "growth", "rate", "r"},
            {"doublingtime", "doubling", "time", "t"},
            {"logisticshape", "halflife", "t50", "time50", "shape"},
            {"spikefactor", "spike", "factor", "f"},
            {"cataclysmtime", "cataclysm", "time", "t"}};

    private String[] argumentNames = new String[]{
            "Population Size", "Ancestral Proportion", "Growth Rate", "Doubling Time", "Logistic Shape", "Spike Factor", "Spike Time"
    };

    private int[][] argumentIndices = {{0}, {0, 2}, {0, 3}, {0, 2, 4}, {0, 3, 4}, {0, 1, 2}, {0, 1, 3}, {0, 2, 5, 6}};

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

    public DemographicDialog(JFrame frame) {
        this.frame = frame;

        demographicCombo = new JComboBox(demographicModels);

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
    }

    public int getDemographicModel() {
        return demographicCombo.getSelectedIndex();
    }

    private int findArgument(JComboBox comboBox, String argument) {
        for (int i = 0; i < comboBox.getItemCount(); i++) {
            String item = ((String) comboBox.getItemAt(i)).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    public int showDialog(TraceList traceList, final TemporalAnalysisFrame temporalAnalysisFrame) {

        for (int i = 0; i < argumentCombos.length; i++) {
            argumentCombos[i].removeAllItems();
            for (int j = 0; j < traceList.getTraceCount(); j++) {
                String statistic = traceList.getTraceName(j);
                argumentCombos[i].addItem(statistic);
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

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Demographic Analysis");
        dialog.pack();

        demographicCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                setArguments(temporalAnalysisFrame);
                dialog.pack();
            }
        });

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        if (result == JOptionPane.OK_OPTION) {
            for (int i = 0; i < argumentCombos.length; i++) {
                argumentTraces[i] = (String) argumentCombos[i].getSelectedItem();
            }
            rootHeightTrace = (String) rootHeightCombo.getSelectedItem();
        }

        return result;
    }

    private void setArguments(TemporalAnalysisFrame temporalAnalysisFrame) {
        optionPanel.removeAll();

        optionPanel.addComponents(new JLabel("Demographic Model:"), demographicCombo);

        JLabel label = new JLabel("<html>Warning! Do not select a model other than that which was<br>" +
                "specified in BEAST to generate the trace being analysed.<br>" +
                "<em>Any other model will produce meaningless results.</em></html>");
        label.setFont(label.getFont().deriveFont(((float) label.getFont().getSize() - 2)));
        optionPanel.addSpanningComponent(label);
        optionPanel.addSeparator();

        optionPanel.addLabel("Select the traces to use for the arguments:");

        int demo = demographicCombo.getSelectedIndex();

        for (int i = 0; i < argumentIndices[demo].length; i++) {
            int k = argumentIndices[demo][i];
            optionPanel.addComponentWithLabel(argumentNames[k] + ":",
                    argumentCombos[k]);
        }

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

    javax.swing.Timer timer = null;

    public void createDemographicFrame(TraceList traceList, DocumentFrame parent) {

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
                        "Error creating Bayesian skyline",
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

        final AnalyseDemographicTask analyseTask = new AnalyseDemographicTask(traceList, frame);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Analysing Demographic Model",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

        timer = new javax.swing.Timer(1000, new ActionListener() {
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

    class AnalyseDemographicTask extends LongTask {

        TraceList traceList;
        TemporalAnalysisFrame frame;
        int binCount;
        boolean rangeSet;
        double minTime;
        double maxTime;
        double ageOfYoungest;

        private int lengthOfTask = 0;
        private int current = 0;

        public AnalyseDemographicTask(TraceList traceList, TemporalAnalysisFrame frame) {
            this.traceList = traceList;
            this.frame = frame;
            this.binCount = frame.getBinCount();
            this.rangeSet = frame.isRangeSet();

            ageOfYoungest = ageOfYoungestField.getValue();

        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Calculating demographic reconstruction...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {


            int n = traceList.getStateCount();

            current = 0;

            int[] argIndices = argumentIndices[demographicCombo.getSelectedIndex()];
            double[][] values = new double[argIndices.length][n];

            Variate[] bins = new Variate[binCount];
            for (int k = 0; k < binCount; k++) {
                bins[k] = new Variate.Double();
            }

            int index = traceList.getTraceIndex(rootHeightTrace);
            double[] heights = new double[n];
            traceList.getValues(index, heights);

            TraceDistribution distribution = new TraceDistribution(heights, traceList.getStepSize());

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
            String title = "";

            for (int j = 0; j < argIndices.length; j++) {
                index = traceList.getTraceIndex(argumentTraces[argIndices[j]]);
                traceList.getValues(index, values[j]);
            }

            if (demographicCombo.getSelectedIndex() == 0) { // Constant Size
                title = "Constant Population Size";
                ConstantPopulation demo = new ConstantPopulation();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else
            if (demographicCombo.getSelectedIndex() == 1) { // Exponential Growth (Growth Rate)
                title = "Exponential Growth";
                ExponentialGrowth demo = new ExponentialGrowth();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setGrowthRate(values[1][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);

                    current++;
                }


            } else
            if (demographicCombo.getSelectedIndex() == 2) { // Exponential Growth (Doubling Time)
                title = "Exponential Growth";
                ExponentialGrowth demo = new ExponentialGrowth();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setDoublingTime(values[1][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else
            if (demographicCombo.getSelectedIndex() == 3) { // Logistic Growth (Growth Rate)
                title = "Logistic Growth";
                LogisticGrowth demo = new LogisticGrowth();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setGrowthRate(values[1][i]);
                    demo.setTime50(values[2][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else
            if (demographicCombo.getSelectedIndex() == 4) { // Logistic Growth (Doubling Time)
                title = "Logistic Growth";
                LogisticGrowth demo = new LogisticGrowth();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setDoublingTime(values[1][i]);
                    demo.setTime50(values[2][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else
            if (demographicCombo.getSelectedIndex() == 5) { // Expansion (Growth Rate)
                title = "Expansion";
                Expansion demo = new Expansion();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setProportion(values[1][i]);
                    demo.setGrowthRate(values[2][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else
            if (demographicCombo.getSelectedIndex() == 6) { // Expansion (Doubling Time)
                title = "Expansion";
                Expansion demo = new Expansion();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setProportion(values[1][i]);
                    demo.setDoublingTime(values[2][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 7) { // Cataclysm
                title = "Boom-Bust";
                CataclysmicDemographic demo = new CataclysmicDemographic();
                for (int i = 0; i < n; i++) {
                    demo.setN0(values[0][i]);
                    demo.setGrowthRate(values[1][i]);
                    demo.setCataclysmTime(values[3][i]);
                    demo.setSpikeFactor(values[2][i]);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            }
            Variate xData = new Variate.Double();
            Variate yDataMean = new Variate.Double();
            Variate yDataMedian = new Variate.Double();
            Variate yDataUpper = new Variate.Double();
            Variate yDataLower = new Variate.Double();

            double t;
            if (ageOfYoungest > 0.0) {
                t = maxTime;
            } else {
                t = minTime;
            }
            for (Variate bin : bins) {
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

            frame.addDemographic(title + ": " + traceList.getName(), xData,
                    yDataMean, yDataMedian,
                    yDataUpper, yDataLower,
                    timeMean, timeMedian,
                    timeUpper, timeLower);

            return null;
        }

        private void addDemographic(Variate[] bins, int binCount, double maxHeight, double delta, DemographicFunction demo) {
            double height;
            if (ageOfYoungest > 0.0) {
                height = ageOfYoungest - maxTime;
            } else {
                height = ageOfYoungest;
            }

            for (int k = 0; k < binCount; k++) {
                if (height >= 0.0 && height <= maxHeight) {
                    bins[k].add(demo.getDemographic(height));
                }
                height += delta;
            }
            current++;
        }
    }
}