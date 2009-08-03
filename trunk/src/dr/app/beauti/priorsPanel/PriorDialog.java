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

package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.Parameter;
import dr.gui.chart.Axis;
import dr.gui.chart.JChart;
import dr.gui.chart.LinearAxis;
import dr.gui.chart.PDFPlot;
import dr.math.distributions.*;
import dr.util.NumberFormatter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.HashMap;
import java.util.Map;

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
			PriorType.LAPLACE_PRIOR,
			PriorType.NORMAL_PRIOR,
			PriorType.TRUNC_NORMAL_PRIOR,
			PriorType.LOGNORMAL_PRIOR,
			PriorType.GAMMA_PRIOR,
			PriorType.JEFFREYS_PRIOR,
	};

	public static PriorType[] rootHeightPriors = {
			PriorType.NONE,
			PriorType.UNIFORM_PRIOR,
			PriorType.EXPONENTIAL_PRIOR,
			PriorType.LAPLACE_PRIOR,
			PriorType.NORMAL_PRIOR,
			PriorType.TRUNC_NORMAL_PRIOR,
			PriorType.LOGNORMAL_PRIOR,
			PriorType.GAMMA_PRIOR,
			PriorType.JEFFREYS_PRIOR,
	};

	private Map<PriorType, PriorOptionsPanel> optionsPanels = new HashMap<PriorType, PriorOptionsPanel>();

	private JComboBox priorCombo;
	private JComboBox rootHeightPriorCombo;
	private RealNumberField initialField = new RealNumberField();

	private OptionsPanel optionPanel;
	private JChart chart;
	private JLabel quantileLabels;
	private JTextArea quantileText;

	private Parameter parameter;

	public PriorDialog(JFrame frame) {
		this.frame = frame;

		priorCombo = new JComboBox(priors);
		rootHeightPriorCombo = new JComboBox(rootHeightPriors);

		initialField.setColumns(8);

		optionsPanels.put(PriorType.UNIFORM_PRIOR, new UniformOptionsPanel());
		optionsPanels.put(PriorType.LAPLACE_PRIOR, new LaplaceOptionsPanel());
		optionsPanels.put(PriorType.NORMAL_PRIOR, new NormalOptionsPanel());
		optionsPanels.put(PriorType.TRUNC_NORMAL_PRIOR, new TruncatedNormalOptionsPanel());
		optionsPanels.put(PriorType.LOGNORMAL_PRIOR, new LogNormalOptionsPanel());
		optionsPanels.put(PriorType.GAMMA_PRIOR, new GammaOptionsPanel());
		optionsPanels.put(PriorType.EXPONENTIAL_PRIOR, new ExponentialOptionsPanel());
//        optionsPanels.put(PriorType.GMRF_PRIOR, new GMRFOptionsPanel());

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

	public int showDialog(final Parameter parameter) {

		this.parameter = parameter;
		PriorType priorType = parameter.priorType;

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

		KeyListener listener = new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				setupChart();
				dialog.repaint();
			}
		};

		for (PriorOptionsPanel optionsPanel : optionsPanels.values()) {
			for (RealNumberField field : optionsPanel.getFields()) {
				field.addKeyListener(listener);
			}
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

		optionsPanels.get(PriorType.UNIFORM_PRIOR).getField(0).setRange(parameter.lower, parameter.upper);
		optionsPanels.get(PriorType.UNIFORM_PRIOR).getField(0).setValue(parameter.uniformLower);
		optionsPanels.get(PriorType.UNIFORM_PRIOR).getField(1).setRange(parameter.lower, parameter.upper);
		optionsPanels.get(PriorType.UNIFORM_PRIOR).getField(1).setValue(parameter.uniformUpper);

		optionsPanels.get(PriorType.EXPONENTIAL_PRIOR).getField(0).setRange(0.0, Double.MAX_VALUE);
		optionsPanels.get(PriorType.EXPONENTIAL_PRIOR).getField(0).setValue(parameter.exponentialMean);
		optionsPanels.get(PriorType.EXPONENTIAL_PRIOR).getField(1).setValue(parameter.exponentialOffset);

		optionsPanels.get(PriorType.NORMAL_PRIOR).getField(0).setValue(parameter.normalMean);
		optionsPanels.get(PriorType.NORMAL_PRIOR).getField(1).setRange(0.0, Double.MAX_VALUE);
		optionsPanels.get(PriorType.NORMAL_PRIOR).getField(1).setValue(parameter.normalStdev);

		optionsPanels.get(PriorType.LOGNORMAL_PRIOR).getField(0).setValue(parameter.logNormalMean);
		optionsPanels.get(PriorType.LOGNORMAL_PRIOR).getField(1).setValue(parameter.logNormalStdev);
		optionsPanels.get(PriorType.LOGNORMAL_PRIOR).getField(2).setValue(parameter.logNormalOffset);

		optionsPanels.get(PriorType.GAMMA_PRIOR).getField(0).setValue(parameter.gammaAlpha);
		optionsPanels.get(PriorType.GAMMA_PRIOR).getField(0).setRange(0.0, Double.MAX_VALUE);
		optionsPanels.get(PriorType.GAMMA_PRIOR).getField(1).setValue(parameter.gammaBeta);
		optionsPanels.get(PriorType.GAMMA_PRIOR).getField(1).setRange(0.0, Double.MAX_VALUE);
		optionsPanels.get(PriorType.GAMMA_PRIOR).getField(2).setValue(parameter.gammaOffset);
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

		if (initialField.getValue() != null) parameter.initial = initialField.getValue();

		if (parameter.priorType != PriorType.JEFFREYS_PRIOR) optionsPanels.get(parameter.priorType).setParameterPrior(parameter);
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
			if (!parameter.priorFixed) {
				optionPanel.addComponents(new JLabel("Prior Distribution:"), priorCombo);
				priorType = (PriorType) priorCombo.getSelectedItem();
			} else {
				priorType = (PriorType) priorCombo.getSelectedItem();
				optionPanel.addComponents(new JLabel("Prior Distribution: "), new JLabel(priorType.toString()));
			}
		}

		if (priorType != PriorType.JEFFREYS_PRIOR) {
			optionPanel.addSeparator();

			optionPanel.addComponent(optionsPanels.get(priorType));
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

		double offset = 0.0;
		Distribution distribution = optionsPanels.get(priorType).getDistribution();

		chart.addPlot(new PDFPlot(distribution, offset));
		if (distribution != null) {
			quantileText.setText(formatter.format(distribution.quantile(0.025)) +
					"\n" + formatter.format(distribution.quantile(0.05)) +
					"\n" + formatter.format(distribution.quantile(0.5)) +
					"\n" + formatter.format(distribution.quantile(0.95)) +
					"\n" + formatter.format(distribution.quantile(0.975)));
		}
	}

	// options panels

	class LaplaceOptionsPanel extends PriorOptionsPanel {

		public LaplaceOptionsPanel() {
			addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
		}

		public Distribution getDistribution() {
			return new LaplaceDistribution(getValue(0), getValue(1));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.laplaceMean = getValue(0);
			parameter.laplaceStdev = getValue(1);
		}
	}

	class UniformOptionsPanel extends PriorOptionsPanel {

		public UniformOptionsPanel() {
			addField("Lower", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			addField("Upper", 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
		}

		public Distribution getDistribution() {
			return new UniformDistribution(getValue(0), getValue(1));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.uniformLower = getValue(0);
			parameter.uniformUpper = getValue(1);
		}

	}

	class ExponentialOptionsPanel extends PriorOptionsPanel {

		public ExponentialOptionsPanel() {
			addField("Mean", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
			addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
		}

		public Distribution getDistribution() {
			return new OffsetPositiveDistribution(
					new ExponentialDistribution(1.0 / getValue(0)), getValue(1));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.exponentialMean = getValue(0);
			parameter.exponentialOffset = getValue(1);
		}
	}

	class NormalOptionsPanel extends PriorOptionsPanel {

		public NormalOptionsPanel() {

			addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			addField("Stdev", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
		}

		public Distribution getDistribution() {
			return new NormalDistribution(getValue(0), getValue(1));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.normalMean = getValue(0);
			parameter.normalStdev = getValue(1);
		}
	}

	class TruncatedNormalOptionsPanel extends PriorOptionsPanel {

		public TruncatedNormalOptionsPanel() {

			addField("Mean", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			addField("Stdev", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
			addField("Lower", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			addField("Upper", 1.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
		}

		public Distribution getDistribution() {
			return new TruncatedNormalDistribution(getValue(0), getValue(1), getValue(2), getValue(3));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.normalMean = getValue(0);
			parameter.normalStdev = getValue(1);
			parameter.uniformLower = getValue(2);
			parameter.uniformUpper = getValue(3);
		}
	}

	class LogNormalOptionsPanel extends PriorOptionsPanel {

		public LogNormalOptionsPanel() {
			addField("M", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
			addField("S", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
			addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
		}

		public Distribution getDistribution() {
			return new OffsetPositiveDistribution(
					new LogNormalDistribution(getValue(0), getValue(1)), getValue(2));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.logNormalMean = getValue(0);
			parameter.logNormalStdev = getValue(1);
			parameter.logNormalOffset = getValue(2);
		}

	}

	class GammaOptionsPanel extends PriorOptionsPanel {

		public GammaOptionsPanel() {
			addField("Shape", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
			addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
			addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
		}

		public Distribution getDistribution() {
			return new OffsetPositiveDistribution(
					new GammaDistribution(getValue(0), getValue(1)), getValue(2));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.gammaAlpha = getValue(0);
			parameter.gammaBeta = getValue(1);
			parameter.gammaOffset = getValue(2);
		}
	}
}