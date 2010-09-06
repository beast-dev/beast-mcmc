/*
 * PriorDialog.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.oldbeauti;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.chart.*;
import dr.util.NumberFormatter;
import dr.math.distributions.*;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorDialog {

    private JFrame frame;

    public static PriorType[] priors = {
            PriorType.UNIFORM_PRIOR,
            PriorType.EXPONENTIAL_PRIOR,
            PriorType.NORMAL_PRIOR,
            PriorType.LOGNORMAL_PRIOR,
            PriorType.GAMMA_PRIOR,
            PriorType.JEFFREYS_PRIOR,
    };

    public static PriorType[] rootHeightPriors = {
            PriorType.NONE,
            PriorType.UNIFORM_PRIOR,
            PriorType.EXPONENTIAL_PRIOR,
            PriorType.NORMAL_PRIOR,
            PriorType.LOGNORMAL_PRIOR,
            PriorType.GAMMA_PRIOR,
            PriorType.JEFFREYS_PRIOR,
    };

    private String[] argumentNames = new String[]{
            "Lower Bound", "Upper Bound", "Exponential Mean", "Zero Offset", "Normal Mean", "Normal Stdev",
            "LogNormal Mean", "LogNormal Stdev", "Zero Offset",
            "Gamma Shape (alpha)", "Gamma Scale (beta)", "Zero Offset",
    };

    private JComboBox priorCombo;
    private JComboBox rootHeightPriorCombo;
    private int[][] argumentIndices = {{0, 1}, {2, 3}, {4, 5}, {6, 7, 8}, {9, 10, 11}, {}, {}, {4, 5, 0, 1}};
    private RealNumberField initialField = new RealNumberField();
    private RealNumberField[] argumentFields = new RealNumberField[argumentNames.length];
    private OptionsPanel optionPanel;
    private JChart chart;
    private JLabel quantileLabels;
    private JTextArea quantileText;

    private JCheckBox truncatedCheck;
    private boolean isTruncated;

    private BeautiOptions.Parameter parameter;

    public PriorDialog(JFrame frame) {
        this.frame = frame;

        priorCombo = new JComboBox(priors);
        rootHeightPriorCombo = new JComboBox(rootHeightPriors);

        truncatedCheck = new JCheckBox("Use a truncated normal distribution");
        truncatedCheck.setOpaque(false);

        initialField.setColumns(8);
        for (int i = 0; i < argumentNames.length; i++) {
            argumentFields[i] = new RealNumberField();
            argumentFields[i].setColumns(8);
        }

        optionPanel = new OptionsPanel(12, 12);

        chart = new JChart(new LinearAxis(Axis.AT_MINOR_TICK, Axis.AT_MINOR_TICK),
                new LinearAxis(Axis.AT_ZERO, Axis.AT_DATA));

        quantileLabels = new JLabel();
        quantileLabels.setFont(quantileLabels.getFont().deriveFont(10.0f));
        quantileLabels.setOpaque(false);
        quantileLabels.setText("<html><p align=\"right\">Quantiles: 2.5%:<br>5%:<br>Median:<br>95%:<br>97.5%:</p></html>");
        quantileText = new JTextArea(0, 5);
        quantileText.setFont(quantileText.getFont().deriveFont(10.0f));
        quantileText.setOpaque(false);
        quantileText.setEditable(false);
    }

    public int showDialog(final BeautiOptions.Parameter parameter) {

        PriorType priorType;

        this.parameter = parameter;

        priorType = parameter.priorType;

        if (parameter.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
            isTruncated = true;
            truncatedCheck.setSelected(isTruncated);
            priorType = PriorType.NORMAL_PRIOR;
        }

        if (parameter.isNodeHeight) {
            if (priorType != PriorType.NONE) {
                rootHeightPriorCombo.setSelectedItem(priorType);
            } else {
                rootHeightPriorCombo.setSelectedIndex(0);
            }
        } else {
            priorCombo.setSelectedItem(priorType);
        }

        if (!parameter.isStatistic) {
            initialField.setRange(parameter.lower, parameter.upper);
            initialField.setValue(parameter.initial);
        }


        setArguments();
        setupComponents();

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Prior for Parameter");
        dialog.pack();

        priorCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                setupComponents();
                dialog.pack();
                dialog.repaint();
            }
        });

        rootHeightPriorCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                setupComponents();
                dialog.pack();
                dialog.repaint();
            }
        });

        truncatedCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                isTruncated = truncatedCheck.isSelected();
                setupComponents();
                dialog.pack();
                dialog.repaint();
            }
        });

        KeyListener listener = new KeyAdapter() {
            public void keyReleased(KeyEvent e) {
                setupChart();
                dialog.repaint();
            }
        };
        for (int i = 0; i < argumentNames.length; i++) {
            argumentFields[i].addKeyListener(listener);
        }

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        if (result == JOptionPane.OK_OPTION) {
            getArguments();
        }

        return result;
    }

    private void setArguments() {
        argumentFields[0].setRange(parameter.lower, parameter.upper);
        argumentFields[0].setValue(parameter.uniformLower);
        argumentFields[1].setRange(parameter.lower, parameter.upper);
        argumentFields[1].setValue(parameter.uniformUpper);

        argumentFields[2].setRange(0.0, Double.MAX_VALUE);
        argumentFields[2].setValue(parameter.exponentialMean);
        argumentFields[3].setValue(parameter.exponentialOffset);

        argumentFields[4].setValue(parameter.normalMean);
        argumentFields[5].setRange(0.0, Double.MAX_VALUE);
        argumentFields[5].setValue(parameter.normalStdev);

        argumentFields[6].setValue(parameter.logNormalMean);
        argumentFields[7].setValue(parameter.logNormalStdev);
        argumentFields[8].setValue(parameter.logNormalOffset);

        argumentFields[9].setRange(0.0, Double.MAX_VALUE);
        argumentFields[9].setValue(parameter.gammaAlpha);
        argumentFields[10].setRange(0.0, Double.MAX_VALUE);
        argumentFields[10].setValue(parameter.gammaBeta);
        argumentFields[11].setValue(parameter.gammaOffset);
    }

    private void getArguments() {
        if (parameter.isNodeHeight) {
            if (rootHeightPriorCombo.getSelectedIndex() == 0) {
                parameter.priorType = PriorType.NONE;
                parameter.initial = Double.NaN;
                return;
            } else {
                parameter.priorType = (PriorType) rootHeightPriorCombo.getSelectedItem();
            }
        } else {
            parameter.priorType = (PriorType) priorCombo.getSelectedItem();
        }

        if (parameter.priorType == PriorType.NORMAL_PRIOR && isTruncated)
            parameter.priorType = PriorType.TRUNC_NORMAL_PRIOR;

        if (initialField.getValue() != null) parameter.initial = initialField.getValue();

        switch (parameter.priorType) {
            case UNIFORM_PRIOR:
                if (argumentFields[0].getValue() != null) parameter.uniformLower = argumentFields[0].getValue();
                if (argumentFields[1].getValue() != null) parameter.uniformUpper = argumentFields[1].getValue();
                break;
            case EXPONENTIAL_PRIOR:
                if (argumentFields[2].getValue() != null) parameter.exponentialMean = argumentFields[2].getValue();
                if (argumentFields[3].getValue() != null) parameter.exponentialOffset = argumentFields[3].getValue();
                break;
            case NORMAL_PRIOR:
                if (argumentFields[4].getValue() != null) parameter.normalMean = argumentFields[4].getValue();
                if (argumentFields[5].getValue() != null) parameter.normalStdev = argumentFields[5].getValue();
                break;
            case LOGNORMAL_PRIOR:
                if (argumentFields[6].getValue() != null) parameter.logNormalMean = argumentFields[6].getValue();
                if (argumentFields[7].getValue() != null) parameter.logNormalStdev = argumentFields[7].getValue();
                if (argumentFields[8].getValue() != null) parameter.logNormalOffset = argumentFields[8].getValue();
                break;
            case GAMMA_PRIOR:
                if (argumentFields[9].getValue() != null) parameter.gammaAlpha = argumentFields[9].getValue();
                if (argumentFields[10].getValue() != null) parameter.gammaBeta = argumentFields[10].getValue();
                if (argumentFields[11].getValue() != null) parameter.gammaOffset = argumentFields[11].getValue();
                break;
            case JEFFREYS_PRIOR:
                break;
            case TRUNC_NORMAL_PRIOR:
                if (argumentFields[0].getValue() != null) parameter.uniformLower = argumentFields[0].getValue();
                if (argumentFields[1].getValue() != null) parameter.uniformUpper = argumentFields[1].getValue();
                if (argumentFields[4].getValue() != null) parameter.normalMean = argumentFields[4].getValue();
                if (argumentFields[5].getValue() != null) parameter.normalStdev = argumentFields[5].getValue();
                break;
            default:
                throw new IllegalArgumentException("Unknown prior index");
        }
    }

    private void setupComponents() {
        optionPanel.removeAll();

        optionPanel.addSpanningComponent(new JLabel("Select prior distribution for " + parameter.getName()));

        PriorType priorType;
        if (parameter.isNodeHeight) {
            optionPanel.addComponents(new JLabel("Prior Distribution:"), rootHeightPriorCombo);
            if (rootHeightPriorCombo.getSelectedIndex() == 0) {
                return;
            } else {
                priorType = (PriorType) rootHeightPriorCombo.getSelectedItem();
            }
        } else {
            optionPanel.addComponents(new JLabel("Prior Distribution:"), priorCombo);
            priorType = (PriorType) priorCombo.getSelectedItem();
        }

        if (priorType == PriorType.NORMAL_PRIOR && isTruncated)
            priorType = PriorType.TRUNC_NORMAL_PRIOR;

        if (priorType != PriorType.JEFFREYS_PRIOR) {
            optionPanel.addSeparator();

            if (priorType == PriorType.NORMAL_PRIOR || priorType == PriorType.TRUNC_NORMAL_PRIOR) {
                optionPanel.addComponent(truncatedCheck);
            }

            for (int i = 0; i < argumentIndices[priorType.ordinal() - 1].length; i++) {
                int k = argumentIndices[priorType.ordinal() - 1][i];
                optionPanel.addComponentWithLabel(argumentNames[k] + ":", argumentFields[k]);
            }
        }

        if (!parameter.isStatistic) {
            optionPanel.addSeparator();
            optionPanel.addComponents(new JLabel("Initial Value:"), initialField);
        }

        if (priorType != PriorType.UNIFORM_PRIOR && priorType != PriorType.JEFFREYS_PRIOR) {
            optionPanel.addSeparator();

            setupChart();
            chart.setPreferredSize(new Dimension(300, 200));
            chart.setFontSize(8);
            optionPanel.addSpanningComponent(chart);
            optionPanel.addComponents(quantileLabels, quantileText);
        }
    }

    NumberFormatter formatter = new NumberFormatter(4);

    private void setupChart() {
        chart.removeAllPlots();

        PriorType priorType;
        if (parameter.isNodeHeight) {
            if (rootHeightPriorCombo.getSelectedIndex() == 0) {
                return;
            } else {
                priorType = (PriorType) rootHeightPriorCombo.getSelectedItem();
            }
        } else {
            priorType = (PriorType) priorCombo.getSelectedItem();
        }

        if (priorType == PriorType.NORMAL_PRIOR && isTruncated)
            priorType = PriorType.TRUNC_NORMAL_PRIOR;

        if (priorType == PriorType.TRUNC_NORMAL_PRIOR && !isTruncated)
            priorType = PriorType.NORMAL_PRIOR;

        Distribution distribution = null;
        double offset = 0.0;
        switch (priorType) {
            case UNIFORM_PRIOR:
                return;
            case EXPONENTIAL_PRIOR:
                double exponentialMean = getValue(argumentFields[2].getValue(), 1.0);
                offset = getValue(argumentFields[3].getValue(), 0.0);
                distribution = new ExponentialDistribution(1.0 / exponentialMean);
                break;
            case NORMAL_PRIOR:
                double normalMean = getValue(argumentFields[4].getValue(), 0.0);
                double normalStdev = getValue(argumentFields[5].getValue(), 1.0);
                distribution = new NormalDistribution(normalMean, normalStdev);
                break;
            case LOGNORMAL_PRIOR:
                double logNormalMean = getValue(argumentFields[6].getValue(), 0.0);
                double logNormalStdev = getValue(argumentFields[7].getValue(), 1.0);
                offset = getValue(argumentFields[8].getValue(), 0.0);
                distribution = new LogNormalDistribution(logNormalMean, logNormalStdev);
                break;
            case GAMMA_PRIOR:
                double gammaAlpha = getValue(argumentFields[9].getValue(), 1.0);
                double gammaBeta = getValue(argumentFields[10].getValue(), 1.0);
                offset = getValue(argumentFields[11].getValue(), 0.0);
                distribution = new GammaDistribution(gammaAlpha, gammaBeta);
                break;
            case JEFFREYS_PRIOR:
                break;
            case TRUNC_NORMAL_PRIOR:
                double truncNormalMean = getValue(argumentFields[4].getValue(), 0.0);
                double truncNormalStdev = getValue(argumentFields[5].getValue(), 1.0);
                double truncLower = getValue(argumentFields[0].getValue(), 0.0);
                double truncUpper = getValue(argumentFields[1].getValue(), 1.0);
                distribution = new TruncatedNormalDistribution(truncNormalMean, truncNormalStdev, truncLower, truncUpper);
                break;
            default:
                throw new IllegalArgumentException("Unknown prior index");

        }
        chart.addPlot(new PDFPlot(distribution, offset));
        if (distribution != null) {
            quantileText.setText(formatter.format(distribution.quantile(0.025)) +
                    "\n" + formatter.format(distribution.quantile(0.05)) +
                    "\n" + formatter.format(distribution.quantile(0.5)) +
                    "\n" + formatter.format(distribution.quantile(0.95)) +
                    "\n" + formatter.format(distribution.quantile(0.975)));
        }
    }

    private double getValue(Double field, double defaultValue) {
        if (field != null) {
            return field;
        }
        return defaultValue;
    }
}
