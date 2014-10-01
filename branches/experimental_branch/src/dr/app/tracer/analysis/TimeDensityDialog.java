/*
 * TimeDensityDialog.java
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
import dr.app.gui.util.LongTask;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

public class TimeDensityDialog {

    private JFrame frame;

    private JComboBox traceCombo;
    private String timeTrace = "None selected";

    private RealNumberField ageOfYoungestField = new RealNumberField();

    private OptionsPanel optionPanel;

    public TimeDensityDialog(JFrame frame) {
        this.frame = frame;

        traceCombo = new JComboBox();

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

        setArguments();

        for (int j = 0; j < traceList.getTraceCount(); j++) {
            String statistic = traceList.getTraceName(j);
            traceCombo.addItem(statistic);
        }
        int index = findArgument(traceCombo, timeTrace);
        if (index == -1) index = findArgument(traceCombo, "root");
        if (index == -1) index = findArgument(traceCombo, "height");
        if (index == -1) index = findArgument(traceCombo, "time");
        if (index == -1) index = 0;
        traceCombo.setSelectedIndex(index);

        final JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Time Density Analysis");
        dialog.pack();

        int result = JOptionPane.CANCEL_OPTION;

        dialog.setVisible(true);

        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }


        return result;
    }

    private void setArguments() {
        optionPanel.removeAll();

        optionPanel.addComponentWithLabel("Select the trace:", traceCombo);

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

    public void addToTemporalAnalysis(TraceList traceList, TemporalAnalysisFrame frame) {

        final AnalyseTimeDensityTask analyseTask = new AnalyseTimeDensityTask(traceList, frame);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Analysing Time-Density",
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

    class AnalyseTimeDensityTask extends LongTask {

        TraceList traceList;
        TemporalAnalysisFrame frame;
        double ageOfYoungest;
        int binCount;
        double minTime;
        double maxTime;

        int stateCount;

        private int lengthOfTask = 0;
        private int current = 0;

        public AnalyseTimeDensityTask(TraceList traceList, TemporalAnalysisFrame frame) {
            this.traceList = traceList;
            this.frame = frame;

            this.binCount = frame.getBinCount();

            lengthOfTask = traceList.getStateCount() + binCount;

            ageOfYoungest = ageOfYoungestField.getValue();

            stateCount = traceList.getStateCount();

        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Calculating Density...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {

            List<Double> times = traceList.getValues(traceList.getTraceIndex((String) traceCombo.getSelectedItem()));

            minTime = frame.getMinTime();
            maxTime = frame.getMaxTime();

            double delta = (maxTime - minTime) / (binCount - 1);

            FrequencyDistribution frequency = new FrequencyDistribution(minTime, binCount, delta);

            for (int i = 0; i < times.size(); i++) {
                frequency.addValue(getTime(times.get(i)));
            }

            Variate.D xData = new Variate.D();
            Variate.D yData = new Variate.D();

            double x = frequency.getLowerBound() - frequency.getBinSize();
            xData.add(x + (frequency.getBinSize() / 2.0));
            yData.add(0.0);
            x += frequency.getBinSize();

            for (int i = 0; i < frequency.getBinCount(); i++) {
                xData.add(x + (frequency.getBinSize() / 2.0));
                double density = frequency.getFrequency(i) / frequency.getBinSize() / times.size();
                yData.add(density);
                x += frequency.getBinSize();
            }

            xData.add(x + (frequency.getBinSize() / 2.0));
            yData.add(0.0);

            frame.addDensity("Density: " + traceList.getName(), xData, yData);

            return null;
        }

        private double getTime(double height) {
            if (ageOfYoungest > 0.0) {
                return ageOfYoungest - height;
            } else {
                return ageOfYoungest;
            }
        }


    }
}