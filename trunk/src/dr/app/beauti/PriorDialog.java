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

package dr.app.beauti;

import dr.gui.chart.*;
import dr.math.*;
import dr.util.NumberFormatter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorDialog {

	private JFrame frame;

	public static String[] priors = { "Uniform",
			"Exponential",
			"Normal",
			"LogNormal",
			"Gamma",
			"Jeffreys"};

	public static String[] rootHeightPriors = { "Tree Prior",
			"Uniform",
			"Exponential",
			"Normal",
			"LogNormal",
			"Gamma",
			"Jeffreys"};

	private String[] argumentNames = new String[] {
			"Lower Bound", "Upper Bound", "Exponential Mean", "Zero Offset", "Normal Mean", "Normal Stdev",
			"LogNormal Mean", "LogNormal Stdev", "Zero Offset",
			"Gamma Shape (alpha)", "Gamma Scale (beta)", "Zero Offset"
	};

	private JComboBox priorCombo;
	private JComboBox rootHeightPriorCombo;
	private int[][] argumentIndices = { {0, 1}, {2, 3}, {4, 5}, {6, 7, 8}, {9, 10, 11}, {} };
	private RealNumberField initialField = new RealNumberField();
	private RealNumberField[] argumentFields = new RealNumberField[argumentNames.length];
	private OptionsPanel optionPanel;
	private JChart chart;
	private JLabel quantileLabels;
	private JTextArea quantileText;

	private BeautiOptions.Parameter parameter;

	public PriorDialog(JFrame frame) {
		this.frame = frame;

		priorCombo = new JComboBox(priors);
		rootHeightPriorCombo = new JComboBox(rootHeightPriors);

		initialField.setColumns(8);
		for (int i = 0; i < argumentNames.length; i++) {
			argumentFields[i] = new RealNumberField();
			argumentFields[i].setColumns(8);
		}

		optionPanel = new OptionsPanel(12,12);

		chart = new JChart(new LinearAxis(), new LinearAxis());

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

		this.parameter = parameter;

		if (parameter.isNodeHeight) {
			if (parameter.priorType != BeautiOptions.NONE) {
				rootHeightPriorCombo.setSelectedIndex(parameter.priorType + 1);
			} else {
				rootHeightPriorCombo.setSelectedIndex(0);
			}
		} else {
			priorCombo.setSelectedIndex(parameter.priorType);
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
			}});

		rootHeightPriorCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				setupComponents();
				dialog.pack();
				dialog.repaint();
			}});

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
		Integer value = (Integer)optionPane.getValue();
		if (value != null && value.intValue() != -1) {
			result = value.intValue();
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
				parameter.priorType = BeautiOptions.NONE;
				parameter.initial = Double.NaN;
				return;
			} else {
				parameter.priorType = rootHeightPriorCombo.getSelectedIndex() - 1;
			}
		} else {
			parameter.priorType = priorCombo.getSelectedIndex();
		}

		if (initialField.getValue() != null) parameter.initial = initialField.getValue().doubleValue();

		switch (parameter.priorType) {
			case BeautiOptions.UNIFORM_PRIOR:
				if (argumentFields[0].getValue() != null) parameter.uniformLower = argumentFields[0].getValue().doubleValue();
				if (argumentFields[1].getValue() != null) parameter.uniformUpper = argumentFields[1].getValue().doubleValue();
				break;
			case BeautiOptions.EXPONENTIAL_PRIOR:
				if (argumentFields[2].getValue() != null) parameter.exponentialMean = argumentFields[2].getValue().doubleValue();
				if (argumentFields[3].getValue() != null) parameter.exponentialOffset = argumentFields[3].getValue().doubleValue();
				break;
			case BeautiOptions.NORMAL_PRIOR:
				if (argumentFields[4].getValue() != null) parameter.normalMean = argumentFields[4].getValue().doubleValue();
				if (argumentFields[5].getValue() != null) parameter.normalStdev = argumentFields[5].getValue().doubleValue();
				break;
			case BeautiOptions.LOG_NORMAL_PRIOR:
				if (argumentFields[6].getValue() != null) parameter.logNormalMean = argumentFields[6].getValue().doubleValue();
				if (argumentFields[7].getValue() != null) parameter.logNormalStdev = argumentFields[7].getValue().doubleValue();
				if (argumentFields[8].getValue() != null) parameter.logNormalOffset = argumentFields[8].getValue().doubleValue();
				break;
			case BeautiOptions.GAMMA_PRIOR:
				if (argumentFields[9].getValue() != null) parameter.gammaAlpha = argumentFields[9].getValue().doubleValue();
				if (argumentFields[10].getValue() != null) parameter.gammaBeta = argumentFields[10].getValue().doubleValue();
				if (argumentFields[11].getValue() != null) parameter.gammaOffset = argumentFields[11].getValue().doubleValue();
				break;
			case BeautiOptions.JEFFREYS_PRIOR:
				break;
			default: throw new IllegalArgumentException("Unknown prior index");
		}
	}

	private void setupComponents() {
		optionPanel.removeAll();

		optionPanel.addSpanningComponent(new JLabel("Select prior distribution for " + parameter.getName()));

		int priorType;
		if (parameter.isNodeHeight) {
			optionPanel.addComponents(new JLabel("Prior Distribution:"), rootHeightPriorCombo);
			if (rootHeightPriorCombo.getSelectedIndex() == 0) {
				return;
			} else {
				priorType = rootHeightPriorCombo.getSelectedIndex() - 1;
			}
		} else {
			optionPanel.addComponents(new JLabel("Prior Distribution:"), priorCombo);
			priorType = priorCombo.getSelectedIndex();
		}

		if (priorType != BeautiOptions.JEFFREYS_PRIOR) {
			optionPanel.addSeparator();

			for (int i = 0; i < argumentIndices[priorType].length; i++) {
				int k = argumentIndices[priorType][i];
				optionPanel.addComponentWithLabel(argumentNames[k] + ":", argumentFields[k]);
			}

		}

		if (!parameter.isStatistic) {
			optionPanel.addSeparator();
			optionPanel.addComponents(new JLabel("Initial Value:"), initialField);
		}

		if (priorType != BeautiOptions.UNIFORM_PRIOR && priorType != BeautiOptions.JEFFREYS_PRIOR) {
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

		int priorType;
		if (parameter.isNodeHeight) {
			if (rootHeightPriorCombo.getSelectedIndex() == 0) {
				return;
			} else {
				priorType = rootHeightPriorCombo.getSelectedIndex() - 1;
			}
		} else {
			priorType = priorCombo.getSelectedIndex();
		}
		Distribution distribution = null;
		double offset = 0.0;
		switch (priorType) {
			case BeautiOptions.UNIFORM_PRIOR:
				return;
			case BeautiOptions.EXPONENTIAL_PRIOR:
				double exponentialMean = argumentFields[2].getValue().doubleValue();
				offset = argumentFields[3].getValue().doubleValue();
				distribution = new ExponentialDistribution(1.0 / exponentialMean);
				break;
			case BeautiOptions.NORMAL_PRIOR:
				double normalMean = argumentFields[4].getValue().doubleValue();
				double normalStdev = argumentFields[5].getValue().doubleValue();
				distribution = new NormalDistribution(normalMean, normalStdev);
				break;
			case BeautiOptions.LOG_NORMAL_PRIOR:
				double logNormalMean = argumentFields[6].getValue().doubleValue();
				double logNormalStdev = argumentFields[7].getValue().doubleValue();
				offset = argumentFields[8].getValue().doubleValue();
				distribution = new LogNormalDistribution(logNormalMean, logNormalStdev);
				break;
			case BeautiOptions.GAMMA_PRIOR:
				double gammaAlpha = argumentFields[9].getValue().doubleValue();
				double gammaBeta = argumentFields[10].getValue().doubleValue();
				offset = argumentFields[11].getValue().doubleValue();
				distribution = new GammaDistribution(gammaAlpha, gammaBeta);
				break;
			case BeautiOptions.JEFFREYS_PRIOR:
				break;
			default: throw new IllegalArgumentException("Unknown prior index");
		}
		chart.addPlot(new PDFPlot(distribution, offset));
		quantileText.setText(formatter.format(distribution.quantile(0.025)) +
				"\n" +formatter.format(distribution.quantile(0.05)) +
				"\n" +formatter.format(distribution.quantile(0.5)) +
				"\n" +formatter.format(distribution.quantile(0.95)) +
				"\n" +formatter.format(distribution.quantile(0.975)));
	}
}