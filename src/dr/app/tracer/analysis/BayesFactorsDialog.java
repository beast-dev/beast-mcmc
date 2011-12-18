/*
 * BayesFactorsDialog.java
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

import dr.app.gui.components.WholeNumberField;
import dr.app.gui.util.LongTask;
import dr.inference.trace.MarginalLikelihoodAnalysis;
import dr.inference.trace.TraceList;
import dr.util.TaskListener;
import jam.framework.DocumentFrame;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class BayesFactorsDialog {

    private JFrame frame;

    private JComboBox likelihoodCombo;
    private JComboBox analysisTypeCombo;
    private WholeNumberField bootstrapCountField;

    private String[] likelihoodGuesses = {
            "likelihood", "treelikelihood", "lnl", "lik"
    };

    private String likelihoodTrace = "None selected";

    private OptionsPanel optionPanel;

    public BayesFactorsDialog(JFrame frame) {
        this.frame = frame;

        likelihoodCombo = new JComboBox();
        analysisTypeCombo = new JComboBox();
        bootstrapCountField = new WholeNumberField(0, Integer.MAX_VALUE);
        bootstrapCountField.setValue(1000);

        optionPanel = new OptionsPanel(12, 12);
    }

    private int findArgument(JComboBox comboBox, String argument) {
        for (int i = 0; i < comboBox.getItemCount(); i++) {
            String item = ((String) comboBox.getItemAt(i)).toLowerCase();
            if (item.indexOf(argument) != -1) return i;
        }
        return -1;
    }

    public String getLikelihoodTrace() {
        return likelihoodTrace;
    }

    public int showDialog(List<TraceList> traceLists) {

        setArguments();

        likelihoodCombo.removeAllItems();

        List<String> statistics = new ArrayList<String>();
        TraceList tl = traceLists.get(0);
        for (int j = 0; j < tl.getTraceCount(); j++) {
            statistics.add(tl.getTraceName(j));
        }

        for (int i = 1; i < traceLists.size(); i++) {
            tl = traceLists.get(i);
            Set<String> statistics2 = new HashSet<String>();
            for (int j = 0; j < tl.getTraceCount(); j++) {
                statistics2.add(tl.getTraceName(j));
            }

            statistics.retainAll(statistics2);
        }

        if (statistics.size() == 0) {
            JOptionPane.showMessageDialog(frame,
                    "These trace files don't seem to have any traces in common.",
                    "Unable to perform analysis",
                    JOptionPane.ERROR_MESSAGE);
            return JOptionPane.CANCEL_OPTION;
        }

        for (String statistic : statistics) {
            likelihoodCombo.addItem(statistic);
        }

        int index = -1;
        for (String guess : likelihoodGuesses) {
            if (index != -1) break;

            index = findArgument(likelihoodCombo, guess);
        }
        if (index == -1) index = 0;
        likelihoodCombo.setSelectedIndex(index);

        analysisTypeCombo.removeAllItems();
        analysisTypeCombo.addItem("aicm");
        analysisTypeCombo.addItem("harmonic mean");

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Model Comparison");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        if (result == JOptionPane.OK_OPTION) {
            likelihoodTrace = (String) likelihoodCombo.getSelectedItem();
        }

        return result;
    }

    private void setArguments() {
        optionPanel.removeAll();

        optionPanel.addLabel("Select the trace to analyse:");

        optionPanel.addComponents(new JLabel("Likelihood trace:"), likelihoodCombo);

        optionPanel.addComponents(new JLabel("Analysis type:"), analysisTypeCombo);

        bootstrapCountField.setColumns(12);
        optionPanel.addComponents(new JLabel("Bootstrap replicates:"), bootstrapCountField);

    }

    Timer timer = null;

    public void createBayesFactorsFrame(List<TraceList> traceLists, DocumentFrame parent) {

        String analysisType = analysisTypeCombo.getSelectedItem().toString();
        int bootstrapLength = bootstrapCountField.getValue();

        boolean isAICM = false;
        if (analysisType.equals("aicm")) {
            isAICM = true;
        }

        BayesFactorsFrame frame = new BayesFactorsFrame(parent, "Model Comparison", bootstrapLength > 1, isAICM);
        frame.initialize();
        frame.setVisible(true);

        final MarginalLikelihoodTask analyseTask = new MarginalLikelihoodTask(frame, traceLists, analysisType, bootstrapLength);

        final ProgressMonitor progressMonitor = new ProgressMonitor(frame,
                "Estimating Marginal Likelihoods",
                "", 0, analyseTask.getLengthOfTask());
        progressMonitor.setMillisToPopup(0);
        progressMonitor.setMillisToDecideToPopup(0);

        timer = new Timer(100, new ActionListener() {
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

    class MarginalLikelihoodTask extends LongTask {

        List<TraceList> traceLists;
        BayesFactorsFrame frame;
        String analysisType;
        int bootstrapLength;

        private int lengthOfTask = 0;
        private int current = 0;

        public MarginalLikelihoodTask(BayesFactorsFrame frame, List<TraceList> traceLists, String analysisType, int bootstrapLength) {
            this.traceLists = traceLists;
            this.analysisType = analysisType;
            this.bootstrapLength = bootstrapLength;
            this.frame = frame;
            lengthOfTask = traceLists.size() * 100;
        }

        public int getCurrent() {
            return current;
        }

        public int getLengthOfTask() {
            return lengthOfTask;
        }

        public String getDescription() {
            return "Estimating marginal likelihoods...";
        }

        public String getMessage() {
            return null;
        }

        public Object doWork() {

            current = 0;
            int i = 0;

            for (TraceList traceList : traceLists) {
                final int offset = i * 100;

                int index = traceList.getTraceIndex(likelihoodTrace);
                List<Double> likelihoods = traceList.getValues(index);

                final MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(likelihoods,
                        traceList.getName(), traceList.getBurnIn(), analysisType, bootstrapLength);

                analysis.setTaskListener(new TaskListener() {
                    public void progress(double progress) {
                        current = offset + (int) (progress * 100);
                    }
                });

                // call to force calculation...
                analysis.calculate();

                EventQueue.invokeLater(
                        new Runnable() {
                            public void run() {
                                frame.addMarginalLikelihood(analysis);
                            }
                        });

                i++;
            }

            return null;
        }

    }

}