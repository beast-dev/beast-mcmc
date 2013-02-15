/*
 * DemographicDialog.java
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
import dr.evolution.coalescent.TwoEpochDemographic;
import dr.evolution.util.Units;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import jam.framework.DocumentFrame;
import jam.panels.OptionsPanel;
import jebl.evolution.coalescent.*;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.List;

public class DemographicDialog {

    private JFrame frame;

    private JComboBox demographicCombo;
    private WholeNumberField binCountField;

    public static String[] demographicModels = {
            "Constant Population",
            "Exponential Growth (Growth Rate)",
            "Exponential Growth (Doubling Time)",
            "Logistic Growth (Growth Rate)",
            "Logistic Growth (Doubling Time)",
            "Expansion (Growth Rate)",
            "Expansion (Doubling Time)",
            "Constant-Exponential",
            "Constant-Logistic",
            "Constant-Exponential-Constant",
            "Exponential-Logistic",
            "Boom-Bust",
            "Two Epoch"
    };

    private String[][] argumentGuesses = {
            {"populationsize", "population", "popsize", "n0", "size", "pop"},
            {"ancestralsize", "ancestralproportion", "ancpopsize", "proportion", "ancestral", "n1"},
            {"exponentialgrowthrate", "exponentialrate", "growthrate", "expgrowth", "growth", "rate", "r"},
            {"doublingtime", "doubling", "time", "t"},
            {"logisticshape", "halflife", "t50", "time50", "logt50", "shape"},
            {"spikefactor", "spike", "factor", "f"},
            {"cataclysmtime", "cataclysm", "time", "t"},
            {"transitiontime", "time1", "time", "t1", "t"},
            {"logisticgrowthrate", "logisticgrowth", "loggrowth", "logisticrate"},
            {"populationsize2", "population2", "popsize2", "n02", "size2", "pop2"},
            {"exponentialgrowthrate2", "exponentialrate2", "growthrate2", "expgrowth2", "growth2", "rate2", "r2"}
    };

    private String[] argumentNames = new String[]{
            "Population Size",
            "Ancestral Proportion",
            "Growth Rate",
            "Doubling Time",
            "Logistic Shape",
            "Spike Factor",
            "Spike Time",
            "Transition Time",
            "Logistic Growth Rate",
            "Population Size 2",
            "Growth Rate 2"
    };

    private int[][] argumentIndices = {
            {0},            // const
            {0, 2},         // exp
            {0, 3},         // exp doubling time
            {0, 2, 4},      // logistic
            {0, 3, 4},      // logistic doubling time
            {0, 1, 2},      // expansion
            {0, 1, 3},      // expansion doubling time
            {0, 7, 2},      // const-exp
            {0, 1, 2, 4},   // const-log
            {0, 1, 2, 7},   // const-exp-const
            {0, 2, 4, 7, 8},// exp-logistic
            {0, 2, 5, 6},    // boom bust
            {0, 2, 9, 10, 7}    // Two Epoch
    };

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
//            int n = traceList.getStateCount();
            current = 0;

            int[] argIndices = argumentIndices[demographicCombo.getSelectedIndex()];
            ArrayList<ArrayList> values = new ArrayList<ArrayList>();

            Variate.D[] bins = new Variate.D[binCount];
            for (int k = 0; k < binCount; k++) {
                bins[k] = new Variate.D();
            }

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
            String title = "";

            for (int j = 0; j < argIndices.length; j++) {
                int index = traceList.getTraceIndex(argumentTraces[argIndices[j]]);
                values.add(new ArrayList(traceList.getValues(index)));
            }

            if (demographicCombo.getSelectedIndex() == 0) { // Constant Size
                title = "Constant Population Size";
                ConstantPopulation demo = new ConstantPopulation();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 1) { // Exponential Growth (Growth Rate)
                title = "Exponential Growth";
                ExponentialGrowth demo = new ExponentialGrowth();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setGrowthRate((Double) values.get(1).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);

                    current++;
                }


            } else if (demographicCombo.getSelectedIndex() == 2) { // Exponential Growth (Doubling Time)
                title = "Exponential Growth";
                ExponentialGrowth demo = new ExponentialGrowth();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setDoublingTime((Double) values.get(1).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 3) { // Logistic Growth (Growth Rate)
                title = "Logistic Growth";
                LogisticGrowth demo = new LogisticGrowth();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setGrowthRate((Double) values.get(1).get(i));
                    demo.setTime50((Double) values.get(2).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 4) { // Logistic Growth (Doubling Time)
                title = "Logistic Growth";
                LogisticGrowth demo = new LogisticGrowth();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setDoublingTime((Double) values.get(1).get(i));
                    demo.setTime50((Double) values.get(2).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 5) { // Expansion (Growth Rate)
                title = "Expansion";
                Expansion demo = new Expansion();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setProportion((Double) values.get(1).get(i));
                    demo.setGrowthRate((Double) values.get(2).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 6) { // Expansion (Doubling Time)
                title = "Expansion";
                Expansion demo = new Expansion();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setProportion((Double) values.get(1).get(i));
                    demo.setDoublingTime((Double) values.get(2).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 7) { // ConstExponential Growth
                title = "Constant-Exponential Growth";
                ConstExponential demo = new ConstExponential();

                for (int i = 0; i < values.get(0).size(); i++) {

                    double N0 = (Double) values.get(0).get(i);
                    double time = (Double) values.get(1).get(i);
                    double r = (Double) values.get(2).get(i);

                    double N1 = N0 * Math.exp(-r * time);

                    demo.setN0(N0);
                    demo.setN1(N1);
                    demo.setGrowthRate(r);

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 8) { // ConstLogistic Growth
                title = "Constant-Logistic Growth";
                ConstLogistic demo = new ConstLogistic();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setN1((Double) values.get(1).get(i));
                    demo.setGrowthRate((Double) values.get(2).get(i));
                    demo.setTime50((Double) values.get(3).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 9) { // ConstExpConst
//                title = "Constant-Exponential-Constant";
//                ConstExpConst demo = new ConstExpConst();
//                for (int i = 0; i < values[0].size(); i++) {
//                    demo.setN0(values[0].get(i));
//	                demo.setN1(values[1].get(i));
//                    demo.setGrowthRate(values[2].get(i));
//                    //demo.setTime50(values[3].get(i));
//
//                    addDemographic(bins, binCount, maxHeight, delta, demo);
//                    current++;
//                }

            } else if (demographicCombo.getSelectedIndex() == 10) { // ExpLogistic Growth
                title = "Exponential-Logistic Growth";
                ExponentialLogistic demo = new ExponentialLogistic();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setR2((Double) values.get(1).get(i));
                    demo.setTime50((Double) values.get(2).get(i));
                    demo.setTime((Double) values.get(3).get(i));
                    demo.setGrowthRate((Double) values.get(4).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 11) { // Cataclysm
                title = "Boom-Bust";
                CataclysmicDemographic demo = new CataclysmicDemographic();
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo.setN0((Double) values.get(0).get(i));
                    demo.setGrowthRate((Double) values.get(1).get(i));
                    demo.setCataclysmTime((Double) values.get(3).get(i));
                    demo.setSpikeFactor((Double) values.get(2).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);
                    current++;
                }

            } else if (demographicCombo.getSelectedIndex() == 12) { // Two Epoch
                title = "Two Epoch";
                dr.evolution.coalescent.ExponentialGrowth demo1 = new dr.evolution.coalescent.ExponentialGrowth(Units.Type.SUBSTITUTIONS);
                dr.evolution.coalescent.ExponentialGrowth demo2 = new dr.evolution.coalescent.ExponentialGrowth(Units.Type.SUBSTITUTIONS);
                TwoEpochDemographic demo = new TwoEpochDemographic(demo1, demo2, Units.Type.SUBSTITUTIONS);
                for (int i = 0; i < values.get(0).size(); i++) {
                    demo1.setN0((Double) values.get(0).get(i));
                    demo1.setGrowthRate((Double) values.get(1).get(i));
                    demo2.setN0((Double) values.get(2).get(i));
                    demo2.setGrowthRate((Double) values.get(3).get(i));
                    demo.setTransitionTime((Double) values.get(4).get(i));

                    addDemographic(bins, binCount, maxHeight, delta, demo);

                    current++;
                }
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

        private void addDemographic(Variate[] bins, int binCount, double maxHeight, double delta, dr.evolution.coalescent.DemographicFunction demo) {
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