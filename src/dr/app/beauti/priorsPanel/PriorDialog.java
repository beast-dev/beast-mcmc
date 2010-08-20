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

import dr.app.beauti.enumTypes.PriorType;
import dr.app.beauti.options.Parameter;
import dr.app.util.OSType;
import dr.gui.chart.Axis;
import dr.gui.chart.JChart;
import dr.gui.chart.LinearAxis;
import dr.gui.chart.PDFPlot;
import dr.math.distributions.*;
import dr.util.NumberFormatter;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorDialog {

	private JFrame frame;

	private Map<PriorType, PriorOptionsPanel> optionsPanels = new HashMap<PriorType, PriorOptionsPanel>();

	private JComboBox priorCombo = new JComboBox(EnumSet.range(PriorType.UNIFORM_PRIOR, PriorType.TRUNC_NORMAL_PRIOR).toArray());
	private JComboBox rootHeightPriorCombo = new JComboBox(EnumSet.range(PriorType.NONE, PriorType.TRUNC_NORMAL_PRIOR).toArray());
    private JCheckBox meanInRealSpaceCheck = new JCheckBox();
	private RealNumberField initialField = new RealNumberField();
    private RealNumberField selectedField;

	private JPanel panel;

    private final SpecialNumberPanel specialNumberPanel;
	private JChart chart;
    private JPanel quantilePanel;
    private JTextArea quantileText;

	private Parameter parameter;

	public PriorDialog(JFrame frame) {
		this.frame = frame;
        
		initialField.setColumns(8);

		optionsPanels.put(PriorType.UNIFORM_PRIOR, new UniformOptionsPanel());
		optionsPanels.put(PriorType.LAPLACE_PRIOR, new LaplaceOptionsPanel());
		optionsPanels.put(PriorType.NORMAL_PRIOR, new NormalOptionsPanel());
		optionsPanels.put(PriorType.EXPONENTIAL_PRIOR, new ExponentialOptionsPanel());
		optionsPanels.put(PriorType.LOGNORMAL_PRIOR, new LogNormalOptionsPanel());
		optionsPanels.put(PriorType.GAMMA_PRIOR, new GammaOptionsPanel());
        optionsPanels.put(PriorType.INVERSE_GAMMA_PRIOR, new InverseGammaOptionsPanel());
		optionsPanels.put(PriorType.TRUNC_NORMAL_PRIOR, new TruncatedNormalOptionsPanel());
//        optionsPanels.put(PriorType.GMRF_PRIOR, new GMRFOptionsPanel());

		chart = new JChart(new LinearAxis(Axis.AT_MINOR_TICK, Axis.AT_MINOR_TICK),
				new LinearAxis(Axis.AT_ZERO, Axis.AT_DATA));

		JLabel quantileLabels = new JLabel();
		quantileLabels.setFont(quantileLabels.getFont().deriveFont(10.0f));
		quantileLabels.setOpaque(false);
		quantileLabels.setText("<html><p align=\"right\">Quantiles: 2.5%:<br>5%:<br>Median:<br>95%:<br>97.5%:</p></html>");

        quantileText = new JTextArea(0, 5);
		quantileText.setFont(quantileText.getFont().deriveFont(10.0f));
		quantileText.setOpaque(false);
		quantileText.setEditable(false);
        quantileLabels.setHorizontalAlignment(JLabel.LEFT);

        quantilePanel = new JPanel();
        quantilePanel.add(quantileLabels);
        quantilePanel.add(quantileText);

        specialNumberPanel = new SpecialNumberPanel(this);
        specialNumberPanel.setEnabled(false);
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

        panel = new JPanel(new GridBagLayout());

        setArguments(priorType);
        setupComponents();

        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

		JOptionPane optionPane = new JOptionPane(scrollPane,
				JOptionPane.PLAIN_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				null,
				null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		final JDialog dialog = optionPane.createDialog(frame, "Prior for Parameter");

		priorCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				setupComponents();
                dialog.validate();
				dialog.repaint();
//                dialog.pack();
			}
		});

		rootHeightPriorCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				setupComponents();
                dialog.validate();
				dialog.repaint();
//                dialog.pack();
			}
		});


        meanInRealSpaceCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                PriorOptionsPanel currentPanel = optionsPanels.get(PriorType.LOGNORMAL_PRIOR);

                if (meanInRealSpaceCheck.isSelected()) {
                    currentPanel.replaceFieldName(0, "Mean");
                    if (currentPanel.getValue(0) <= 0) {
                        currentPanel.getField(0).setValue(0.01);
                    }
                    currentPanel.getField(0).setRange(0.0, Double.POSITIVE_INFINITY);
                } else {
                    currentPanel.replaceFieldName(0, "Log(Mean)");
                    currentPanel.getField(0).setRange(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                }

                setupChart();
                dialog.validate();
                dialog.repaint();
//                dialog.pack();
            }
        });

        KeyListener listener = new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				setupChart();
				dialog.repaint();
			}
		};
        FocusListener flistener = new FocusAdapter() {
			public void focusGained(FocusEvent e) {
                if (e.getComponent() instanceof RealNumberField) {
                    selectedField = (RealNumberField) e.getComponent();
                    specialNumberPanel.setEnabled(true);
                }
			}

            public void focusLost(FocusEvent e) {
                selectedField = null;
                specialNumberPanel.setEnabled(false);
            }
        };

		for (PriorOptionsPanel optionsPanel : optionsPanels.values()) {
			for (JComponent component : optionsPanel.getJComponents()) {
                if (component instanceof RealNumberField) {
                    component.addKeyListener(listener);
                    component.addFocusListener(flistener);
                }
			}
		}

//        Toolkit tk = Toolkit.getDefaultToolkit();
//        Dimension d = tk.getScreenSize();
//        if (d.height < 700 && optionsPanel.getHeight() > 450) {
//            dialog.setSize(new java.awt.Dimension(optionsPanel.getWidth() + 100, 550));
//        } else {
//            // setSize because optionsPanel is shrunk in dialog
//            dialog.setSize(new java.awt.Dimension(optionsPanel.getWidth() + 100, optionsPanel.getHeight() + 100));
//        }

//        System.out.println("panel width = " + optionsPanel.getWidth());
//        System.out.println("panel height = " + optionsPanel.getHeight());
        
        dialog.pack();
        dialog.setMinimumSize(new Dimension(dialog.getBounds().width, 300));

        dialog.setResizable(true);
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

	private void setArguments(PriorType priorType) {
        PriorOptionsPanel panel;
        switch (priorType) {
            case UNIFORM_PRIOR:
                panel = optionsPanels.get(priorType);
//                panel.getField(0).setRange(parameter.lower, parameter.upper);
                panel.getField(0).setValue(parameter.lower);
//                panel.getField(1).setRange(parameter.lower, parameter.upper);
                panel.getField(1).setValue(parameter.upper);
                break;

            case EXPONENTIAL_PRIOR:
                panel = optionsPanels.get(priorType);
//                panel.getField(0).setRange(0.0, Double.MAX_VALUE);
                if (parameter.mean != 0) {// ExponentialDistribution(1.0 / mean)
                    panel.getField(0).setValue(parameter.mean);
                }
                panel.getField(1).setValue(parameter.offset);
                break;

            case NORMAL_PRIOR:
                panel = optionsPanels.get(priorType);
                panel.getField(0).setValue(parameter.mean);
//                panel.getField(1).setRange(0.0, Double.MAX_VALUE);
                panel.getField(1).setValue(parameter.stdev);
                break;

            case LOGNORMAL_PRIOR:
                panel = optionsPanels.get(priorType);
                if (parameter.isMeanInRealSpace() && parameter.mean <= 0) {// if LOGNORMAL && meanInRealSpace = true, then mean > 0
                    panel.getField(0).setValue(0.01);
                } else {
                    panel.getField(0).setValue(parameter.mean);
                }
                panel.getField(1).setValue(parameter.stdev);
                panel.getField(2).setValue(parameter.offset);
                meanInRealSpaceCheck.setSelected(parameter.isMeanInRealSpace());
                break;

            case LAPLACE_PRIOR:
                panel = optionsPanels.get(priorType);
                panel.getField(0).setValue(parameter.mean);
                panel.getField(1).setValue(parameter.stdev);
                break;

            case GAMMA_PRIOR:
                panel = optionsPanels.get(priorType);
                panel.getField(0).setValue(parameter.shape);
//                panel.getField(0).setRange(0.0, Double.MAX_VALUE);
                panel.getField(1).setValue(parameter.scale);
//                panel.getField(1).setRange(0.0, Double.MAX_VALUE);
                panel.getField(2).setValue(parameter.offset);
                break;

            case INVERSE_GAMMA_PRIOR:
                panel = optionsPanels.get(priorType);
                panel.getField(0).setValue(parameter.shape);
                panel.getField(1).setValue(parameter.scale);
                panel.getField(2).setValue(parameter.offset);
                break;

            case TRUNC_NORMAL_PRIOR:
                panel = optionsPanels.get(priorType);
                panel.getField(0).setValue(parameter.mean);
                panel.getField(1).setValue(parameter.stdev);
                panel.getField(2).setValue(parameter.lower);
                panel.getField(3).setValue(parameter.upper);
                break;
        }
       
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

		if (!parameter.isStatistic && initialField.getValue() != null) parameter.initial = initialField.getValue();

		if (parameter.priorType != PriorType.JEFFREYS_PRIOR) optionsPanels.get(parameter.priorType).setParameterPrior(parameter);
	}

	private void setupComponents() {
        panel.removeAll();
        
        OptionsPanel optionsPanel = new OptionsPanel(12, (OSType.isMac() ? 6 : 24));

		optionsPanel.addSpanningComponent(new JLabel("Select prior distribution for " + parameter.getName()));

		PriorType priorType;
		if (parameter.isNodeHeight) {
			optionsPanel.addComponentWithLabel("Prior Distribution: ", rootHeightPriorCombo);
			if (rootHeightPriorCombo.getSelectedIndex() == 0) {
				return; // PriorType.NONE
			} else {
				priorType = (PriorType) rootHeightPriorCombo.getSelectedItem();
			}
		} else {
			if (!parameter.priorFixed) {
				optionsPanel.addComponentWithLabel("Prior Distribution: ", priorCombo);
				priorType = (PriorType) priorCombo.getSelectedItem();
			} else {
				priorType = (PriorType) priorCombo.getSelectedItem();
				optionsPanel.addComponentWithLabel("Prior Distribution: ", new JLabel(priorType.toString()));
			}
		}

        if (!parameter.isStatistic) {
            optionsPanel.addSeparator();
            optionsPanel.addComponentWithLabel("Initial Value: ", initialField);
        }
        
		if (priorType != PriorType.JEFFREYS_PRIOR) {
			optionsPanel.addSpanningComponent(optionsPanels.get(priorType));
		}

        if (priorType == PriorType.UNIFORM_PRIOR || priorType == PriorType.TRUNC_NORMAL_PRIOR) {
            optionsPanel.addSeparator();
            optionsPanel.addSpanningComponent(specialNumberPanel);
        }

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 0.0;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.PAGE_START;
        gbc.gridwidth = GridBagConstraints.REMAINDER;

        JPanel panel1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel1.add(optionsPanel);

        panel.add(panel1, gbc);

        // UNIFORM_PRIOR and JEFFREYS_PRIOR have no chart
		if (priorType != PriorType.UNIFORM_PRIOR && priorType != PriorType.JEFFREYS_PRIOR) {
			optionsPanel.addSeparator();

			setupChart();
			chart.setPreferredSize(new Dimension(300, 200));
			chart.setFontSize(8);

            gbc.gridy = 1;
            gbc.weighty = 1.0;
            gbc.fill = GridBagConstraints.BOTH;

            panel.add(chart, gbc);

            gbc.gridy = 2;
            gbc.weighty = 0.0;
            gbc.anchor = GridBagConstraints.PAGE_END;
            gbc.fill = GridBagConstraints.HORIZONTAL;

            panel.add(quantilePanel, gbc);

        }
	}

    public void setSelectedField(RealNumberField selectedField) {
        this.selectedField = selectedField;
    }
    
    public RealNumberField getSelectedField() {
        return selectedField;
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
		// ExponentialDistribution(1.0 / mean)
//        if (priorType == PriorType.EXPONENTIAL_PRIOR && parameter.mean == 0) parameter.mean = 1; 
        
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
			addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE); //TODO Beta?
		}

		public Distribution getDistribution() {
			return new LaplaceDistribution(getValue(0), getValue(1));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.mean = getValue(0);
			parameter.stdev = getValue(1); //TODO  stdev or  Scale
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
			parameter.lower = getValue(0);
			parameter.upper = getValue(1);
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
			parameter.mean = getValue(0);
			parameter.offset = getValue(1);
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
			parameter.mean = getValue(0);
			parameter.stdev = getValue(1);
		}
	}

	class LogNormalOptionsPanel extends PriorOptionsPanel {

		public LogNormalOptionsPanel() {
            if (meanInRealSpaceCheck.isSelected()) {
                addField("Mean", 0.01, 0.0, Double.POSITIVE_INFINITY);
            } else {
                addField("Log(Mean)", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
            }
            addField("Log(Stdev)", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
			addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
            addCheckBox("Mean In Real Space", meanInRealSpaceCheck);            
		}

		public Distribution getDistribution() {
            double mean = getValue(0);
            if (meanInRealSpaceCheck.isSelected()) {
                if (mean <= 0) {
                    throw new IllegalArgumentException("meanInRealSpace works only for a positive mean");
                }
                mean = Math.log(getValue(0)) - 0.5 * getValue(1) * getValue(1);
            }
			return new OffsetPositiveDistribution(
					new LogNormalDistribution(mean, getValue(1)), getValue(2));
		}

		public void setParameterPrior(Parameter parameter) {
			parameter.mean = getValue(0);
			parameter.stdev = getValue(1);
			parameter.offset = getValue(2);
            parameter.setMeanInRealSpace(meanInRealSpaceCheck.isSelected());
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
			parameter.shape = getValue(0);
			parameter.scale = getValue(1);
			parameter.offset = getValue(2);
		}
	}

    class InverseGammaOptionsPanel extends PriorOptionsPanel {

        public InverseGammaOptionsPanel() {
            addField("Shape", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Scale", 1.0, Double.MIN_VALUE, Double.MAX_VALUE);
            addField("Offset", 0.0, 0.0, Double.MAX_VALUE);
        }

        public Distribution getDistribution() {
            return new OffsetPositiveDistribution(
                    new InverseGammaDistribution(getValue(0), getValue(1)), getValue(2));
        }

        public void setParameterPrior(Parameter parameter) {
            parameter.shape = getValue(0);
            parameter.scale = getValue(1);
            parameter.offset = getValue(2);
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
            parameter.mean = getValue(0);
            parameter.stdev = getValue(1);
            parameter.lower = getValue(2);
            parameter.upper = getValue(3);
        }
    }
    
}