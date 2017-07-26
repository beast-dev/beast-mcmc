/*
 * MLEDialog.java
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

package dr.app.beauti.components.marginalLikelihoodEstimation;

import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.WholeNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;


/**
 * @author Guy Baele
 */
public class MLEDialog {

    private JFrame frame;

    private final OptionsPanel optionsPanel;

    private JLabel labelPathSteps, labelChainLength, labelLogEvery, labelLogFileName, labelResultFileName, labelStepDistribution;

    private WholeNumberField pathStepsField = new WholeNumberField(1, Integer.MAX_VALUE);
    private WholeNumberField chainLengthField = new WholeNumberField(1, Integer.MAX_VALUE);
    private WholeNumberField logEveryField = new WholeNumberField(1, Integer.MAX_VALUE);

    private JTextArea logFileNameField = new JTextArea("mle.log");
    private JTextArea resultFileNameField = new JTextArea("mle.result.log");

    JCheckBox operatorAnalysis = new JCheckBox("Print operator analysis");

    private JComboBox stepDistribution = new JComboBox();

    private MarginalLikelihoodEstimationOptions options;

    private String description = "Settings for marginal likelihood estimation using PS/SS";

    public MLEDialog(final JFrame frame, final MarginalLikelihoodEstimationOptions options) {
        this.frame = frame;
        this.options = options;

        optionsPanel = new OptionsPanel(12, 12);

        optionsPanel.setOpaque(false);

        JTextArea mleInfo = new JTextArea("Set the options to perform marginal likelihood " +
                "estimation (MLE) using path sampling (PS) / stepping-stone sampling (SS).");
        mleInfo.setColumns(56);
        PanelUtils.setupComponent(mleInfo);
        optionsPanel.addSpanningComponent(mleInfo);

        pathStepsField.setValue(options.pathSteps);
        pathStepsField.setColumns(16);
        pathStepsField.setMinimumSize(pathStepsField.getPreferredSize());
        labelPathSteps = optionsPanel.addComponentWithLabel("Number of path steps:", pathStepsField);
        pathStepsField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.pathSteps = pathStepsField.getValue();
            }
        });

        chainLengthField.setValue(options.mleChainLength);
        chainLengthField.setColumns(16);
        chainLengthField.setMinimumSize(chainLengthField.getPreferredSize());
        labelChainLength = optionsPanel.addComponentWithLabel("Length of chains:", chainLengthField);
        chainLengthField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.mleChainLength = chainLengthField.getValue();
            }
        });

        optionsPanel.addSeparator();

        logEveryField.setValue(options.mleLogEvery);
        logEveryField.setColumns(16);
        logEveryField.setMinimumSize(logEveryField.getPreferredSize());
        labelLogEvery = optionsPanel.addComponentWithLabel("Log likelihood every:", logEveryField);
        logEveryField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.mleLogEvery = logEveryField.getValue();
            }
        });

        optionsPanel.addSeparator();

        logFileNameField.setColumns(32);
        logFileNameField.setEditable(false);
        logFileNameField.setMinimumSize(logFileNameField.getPreferredSize());
        labelLogFileName = optionsPanel.addComponentWithLabel("Log file name:", logFileNameField);
        logFileNameField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.mleFileName = logFileNameField.getText();
            }
        });

        optionsPanel.addSeparator();

        resultFileNameField.setColumns(32);
        resultFileNameField.setEditable(false);
        resultFileNameField.setMinimumSize(resultFileNameField.getPreferredSize());
        labelResultFileName = optionsPanel.addComponentWithLabel("Results file name:", resultFileNameField);
        resultFileNameField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.mleResultFileName = resultFileNameField.getText();
            }
        });

        optionsPanel.addSeparator();

        JTextArea betaInfo = new JTextArea("By default, the power posteriors are determined according to " +
                "evenly spaced quantiles of a Beta(0.3, 1.0) distribution, thereby estimating " +
                "more power posteriors close to the prior.");
        betaInfo.setColumns(56);
        PanelUtils.setupComponent(betaInfo);
        optionsPanel.addSpanningComponent(betaInfo);

        stepDistribution.addItem("Beta");
        labelStepDistribution = optionsPanel.addComponentWithLabel("Path step distribution:", stepDistribution);

        optionsPanel.addSeparator();

        operatorAnalysis.setSelected(options.printOperatorAnalysis);
        optionsPanel.addComponent(operatorAnalysis);

        operatorAnalysis.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (operatorAnalysis.isSelected()) {
                    options.printOperatorAnalysis = true;
                } else {
                    options.printOperatorAnalysis = false;
                }
            }
        });

        optionsPanel.addSeparator();

        JTextArea mleTutorial = new JTextArea("Additional information on marginal likelihood estimation in BEAST " +
                "can be found on http://beast.community/model_selection");
        mleTutorial.setColumns(56);
        PanelUtils.setupComponent(mleTutorial);
        optionsPanel.addSpanningComponent(mleTutorial);

        JTextArea citationText = new JTextArea("Baele G, Lemey P, Bedford T, Rambaut A, Suchard MA, Alekseyenko AV (2012)\n" +
                "Mol Biol Evol 29(9), 2157-2167 [Advantages of PS/SS].\n" +
                "Baele G, Li WLS, Drummond AJ, Suchard MA, Lemey P (2013)\nMol Biol Evol 30(2), 239-243 " +
                "[Importance of using proper priors].");
        citationText.setColumns(45);
        optionsPanel.addComponentWithLabel("Citation:", citationText);

        optionsPanel.addSeparator();
        
        /*JPanel panel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel.add(optionsPanel, BorderLayout.CENTER);
        panel.setOpaque(false);
        
        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        add(scrollPane, BorderLayout.CENTER);*/

    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(optionsPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, description);
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

	/*public JComponent getExportableComponent() {
		return optionsPanel;
	}*/

    public void setFilenameStem(String fileNameStem, boolean addTxt) {
        logFileNameField.setText(fileNameStem + ".mle.log" + (addTxt ? ".txt" : ""));
        resultFileNameField.setText(fileNameStem + ".mle.result.log" + (addTxt ? ".txt" : ""));
        options.mleFileName = logFileNameField.getText();
        options.mleResultFileName = resultFileNameField.getText();
    }

    public void setOptions(MarginalLikelihoodEstimationOptions options) {
        //this.options = options;

        /*System.err.println("setOptions: " + options);
        System.err.println("options.pathSteps: " + options.pathSteps);
        System.err.println("options.mleChainLength: " + options.mleChainLength);
        System.err.println("options.mleLogEvery: " + options.mleLogEvery);*/

        pathStepsField.setValue(options.pathSteps);
        chainLengthField.setValue(options.mleChainLength);
        logEveryField.setValue(options.mleLogEvery);

        logFileNameField.setText(options.mleFileName);
        resultFileNameField.setText(options.mleResultFileName);

        operatorAnalysis.setSelected(options.printOperatorAnalysis);

        optionsPanel.validate();
        optionsPanel.repaint();
    }

    public void getOptions(MarginalLikelihoodEstimationOptions options) {
        //this.options = options;

        options.pathSteps = pathStepsField.getValue();
        options.mleChainLength = chainLengthField.getValue();
        options.mleLogEvery = logEveryField.getValue();

        options.printOperatorAnalysis = operatorAnalysis.isSelected();

        options.mleFileName = logFileNameField.getText();
        options.mleResultFileName = resultFileNameField.getText();

        /*System.err.println("getOptions: " + options);
        System.err.println("options.pathSteps: " + options.pathSteps);
        System.err.println("options.mleChainLength: " + options.mleChainLength);
        System.err.println("options.mleLogEvery: " + options.mleLogEvery);*/

    }

}
