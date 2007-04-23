/*
 * PriorsPanel.java
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

import dr.util.NumberFormatter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorsPanel extends JPanel implements Exportable {

	/**
	 *
	 */
	private static final long serialVersionUID = -2936049032365493416L;
	JScrollPane scrollPane = new JScrollPane();
	JTable priorTable = null;
	PriorTableModel priorTableModel = null;

	OptionsPanel treePriorPanel = new OptionsPanel();
	JComboBox treePriorCombo;
	JComboBox parameterizationCombo = new JComboBox(new String[] {
			"Growth Rate", "Doubling Time"});
	JComboBox bayesianSkylineCombo = new JComboBox(new String[] {
			"Constant", "Linear"});
	WholeNumberField groupCountField = new WholeNumberField(2, Integer.MAX_VALUE);

	RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

	JCheckBox upgmaStartingTreeCheck = new JCheckBox("Use UPGMA to construct a starting tree");

	public ArrayList parameters = new ArrayList();

	BeautiFrame frame = null;

	public PriorsPanel(BeautiFrame parent) {

		this.frame = parent;

		priorTableModel = new PriorTableModel();
		priorTable = new JTable(priorTableModel);

		priorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
		priorTable.getTableHeader().setReorderingAllowed(false);
		priorTable.getTableHeader().setDefaultRenderer(
				new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

		priorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_LAST_COLUMN);
		priorTable.getColumnModel().getColumn(0).setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
		priorTable.getColumnModel().getColumn(0).setPreferredWidth(160);

		priorTable.getColumnModel().getColumn(1).setCellRenderer(
				new ButtonRenderer(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
		priorTable.getColumnModel().getColumn(1).setCellEditor(
				new ButtonEditor(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
		priorTable.getColumnModel().getColumn(1).setPreferredWidth(260);

		priorTable.getColumnModel().getColumn(2).setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
		priorTable.getColumnModel().getColumn(2).setPreferredWidth(400);

		scrollPane = new JScrollPane(priorTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		scrollPane.setOpaque(false);

		java.awt.event.ItemListener listener = new java.awt.event.ItemListener() {
			public void itemStateChanged(java.awt.event.ItemEvent ev) {
				if (!settingOptions) frame.priorsChanged();
			}
		};

		if (BeautiApp.developer) {
			treePriorCombo = new JComboBox(new String[] {
					"Coalescent: Constant Size",
					"Coalescent: Exponential Growth",
					"Coalescent: Logistic Growth",
					"Coalescent: Expansion Growth",
					"Coalescent: Bayesian Skyline",
					"Speciation: Yule Process",
					"Speciation: Birth-Death Process"
			});
		} else {
			treePriorCombo = new JComboBox(new String[] {
					"Coalescent: Constant Size",
					"Coalescent: Exponential Growth",
					"Coalescent: Logistic Growth",
					"Coalescent: Expansion Growth",
					"Coalescent: Bayesian Skyline",
					"Speciation: Yule Process"
					// Until we have tested the Birth-Death process properly, I have hidden this option
					//"Speciation: Birth-Death Process"
			});
		}
		treePriorCombo.setOpaque(false);
		treePriorCombo.addItemListener(
				new java.awt.event.ItemListener() {
					public void itemStateChanged(java.awt.event.ItemEvent ev) {
						if (!settingOptions) frame.priorsChanged();
						setupPanel();
					}
				}
		);
		groupCountField.addKeyListener(new java.awt.event.KeyAdapter() {
			public void keyTyped(java.awt.event.KeyEvent ev) {
				if (!settingOptions) frame.priorsChanged();
			}});
		samplingProportionField.addKeyListener(new java.awt.event.KeyAdapter() {
			public void keyTyped(java.awt.event.KeyEvent ev) {
				if (!settingOptions) frame.priorsChanged();
			}});

		parameterizationCombo.setOpaque(false);
		parameterizationCombo.addItemListener(listener);

		bayesianSkylineCombo.setOpaque(false);
		bayesianSkylineCombo.addItemListener(listener);

		upgmaStartingTreeCheck.setOpaque(false);

		setOpaque(false);
		setLayout(new BorderLayout(0,0));
		setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

		JPanel panel = new JPanel(new BorderLayout(0,0));
		panel.setOpaque(false);
		panel.add(new JLabel("Priors for model parameters:"), BorderLayout.NORTH);
		panel.add(scrollPane, BorderLayout.CENTER);
		panel.add(new JLabel("* Marked parameters currently have a default prior distribution. " +
				"You could check that these are appropriate."), BorderLayout.SOUTH);

		treePriorPanel.setBorder(null);
		add(treePriorPanel, BorderLayout.NORTH);
		add(panel, BorderLayout.CENTER);
	}

	private void setupPanel() {

		treePriorPanel.removeAll();

		treePriorPanel.addComponentWithLabel("Tree Prior:", treePriorCombo);
		if (treePriorCombo.getSelectedIndex() == 1 || // exponential
				treePriorCombo.getSelectedIndex() == 2 || // logistic
				treePriorCombo.getSelectedIndex() == 3 ) { // expansion
			treePriorPanel.addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
		} else if (treePriorCombo.getSelectedIndex() == 4 ) { // bayesian skyline
			groupCountField.setColumns(6);
			treePriorPanel.addComponentWithLabel("Number of groups:", groupCountField);
			treePriorPanel.addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
		} else if (treePriorCombo.getSelectedIndex() == 6 ) { // birth-death
			samplingProportionField.setColumns(8);
			treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
		}

		treePriorPanel.addComponent(upgmaStartingTreeCheck);

		validate();
		repaint();
	}

	private boolean settingOptions = false;

	public void setOptions(BeautiOptions options) {
		settingOptions = true;
		parameters = options.selectParameters();
		priorTableModel.fireTableDataChanged();

		if (options.nodeHeightPrior == BeautiOptions.CONSTANT) {
			treePriorCombo.setSelectedIndex(0);
		} else if (options.nodeHeightPrior == BeautiOptions.EXPONENTIAL) {
			treePriorCombo.setSelectedIndex(1);
		} else if (options.nodeHeightPrior == BeautiOptions.LOGISTIC) {
			treePriorCombo.setSelectedIndex(2);
		} else if (options.nodeHeightPrior == BeautiOptions.EXPANSION) {
			treePriorCombo.setSelectedIndex(3);
		} else if (options.nodeHeightPrior == BeautiOptions.SKYLINE) {
			treePriorCombo.setSelectedIndex(4);
		} else if (options.nodeHeightPrior == BeautiOptions.YULE) {
			treePriorCombo.setSelectedIndex(5);
		} else if (options.nodeHeightPrior == BeautiOptions.BIRTH_DEATH) {
			treePriorCombo.setSelectedIndex(6);
		}
		groupCountField.setValue(options.skylineGroupCount);
		samplingProportionField.setValue(options.birthDeathSamplingProportion);

		parameterizationCombo.setSelectedIndex(options.parameterization);
		bayesianSkylineCombo.setSelectedIndex(options.skylineModel);

		upgmaStartingTreeCheck.setSelected(options.upgmaStartingTree);

		setupPanel();

		settingOptions = false;

		validate();
		repaint();
	}

	private PriorDialog priorDialog = null;

	private void priorButtonPressed(int row) {
		if (priorDialog == null) {
			priorDialog = new PriorDialog(frame);
		}

		BeautiOptions.Parameter param = (BeautiOptions.Parameter)parameters.get(row);

		if (priorDialog.showDialog(param) == JOptionPane.CANCEL_OPTION) {
			return;
		}

		param.priorEdited = true;

		priorTableModel.fireTableDataChanged();
	}

	public void getOptions(BeautiOptions options) {
		if (settingOptions) return;

		if (treePriorCombo.getSelectedIndex() == 0) {
			options.nodeHeightPrior = BeautiOptions.CONSTANT;
		} else if (treePriorCombo.getSelectedIndex() == 1) {
			options.nodeHeightPrior = BeautiOptions.EXPONENTIAL;
		} else if (treePriorCombo.getSelectedIndex() == 2) {
			options.nodeHeightPrior = BeautiOptions.LOGISTIC;
		} else if (treePriorCombo.getSelectedIndex() == 3) {
			options.nodeHeightPrior = BeautiOptions.EXPANSION;
		} else if (treePriorCombo.getSelectedIndex() == 4) {
			options.nodeHeightPrior = BeautiOptions.SKYLINE;
			Integer groupCount = groupCountField.getValue();
			if (groupCount != null) {
				options.skylineGroupCount = groupCount.intValue();
			} else {
				options.skylineGroupCount = 5;
			}
		} else if (treePriorCombo.getSelectedIndex() == 5) {
			options.nodeHeightPrior = BeautiOptions.YULE;
		} else if (treePriorCombo.getSelectedIndex() == 6) {
			options.nodeHeightPrior = BeautiOptions.BIRTH_DEATH;
			Double samplingProportion = samplingProportionField.getValue();
			if (samplingProportion != null) {
				options.birthDeathSamplingProportion = samplingProportion.doubleValue();
			} else {
				options.birthDeathSamplingProportion = 1.0;
			}
		} else {
			throw new RuntimeException("Unexpected value from treePriorCombo");
		}

		options.parameterization = parameterizationCombo.getSelectedIndex();
		options.skylineModel = bayesianSkylineCombo.getSelectedIndex();

		options.upgmaStartingTree = upgmaStartingTreeCheck.isSelected();
	}

	public JComponent getExportableComponent() {
		return priorTable;
	}

	NumberFormatter formatter = new NumberFormatter(4);

	class PriorTableModel extends AbstractTableModel {

		/**
		 *
		 */
		private static final long serialVersionUID = -8864178122484971872L;
		String[] columnNames = { "Parameter", "Prior", "Description" };

		public PriorTableModel() {
		}

		public int getColumnCount() {
			return columnNames.length;
		}

		public int getRowCount() {
			return parameters.size();
		}

		public Object getValueAt(int row, int col) {
			BeastGenerator.Parameter param = (BeastGenerator.Parameter)parameters.get(row);
			switch (col) {
				case 0: return param.getName();
				case 1: return getPriorString(param);
				case 2: return param.getDescription();
			}
			return null;
		}

		public String getPriorString(BeastGenerator.Parameter param) {
			StringBuffer buffer = new StringBuffer();

			if (!param.priorEdited) {
				buffer.append("* ");
			}
			switch (param.priorType) {
				case BeautiOptions.NONE:
					buffer.append("Using Tree Prior");
					break;
				case BeautiOptions.UNIFORM_PRIOR:
					buffer.append("Uniform [");
					buffer.append(formatter.format(param.uniformLower));
					buffer.append(", ");
					buffer.append(formatter.format(param.uniformUpper));
					buffer.append("]");
					break;
				case BeautiOptions.EXPONENTIAL_PRIOR:
					buffer.append("Exponential [");
					buffer.append(formatter.format(param.exponentialMean));
					buffer.append("]");
					break;
				case BeautiOptions.NORMAL_PRIOR:
					buffer.append("Normal [");
					buffer.append(formatter.format(param.normalMean));
					buffer.append(", ");
					buffer.append(formatter.format(param.normalStdev));
					buffer.append("]");
					break;
				case BeautiOptions.LOG_NORMAL_PRIOR:
					buffer.append("LogNormal [");
					buffer.append(formatter.format(param.logNormalMean));
					buffer.append(", ");
					buffer.append(formatter.format(param.logNormalStdev));
					buffer.append("]");
					break;
				case BeautiOptions.GAMMA_PRIOR:
					buffer.append("Gamma [");
					buffer.append(formatter.format(param.gammaAlpha));
					buffer.append(", ");
					buffer.append(formatter.format(param.gammaBeta));
					buffer.append("]");
					break;
				case BeautiOptions.JEFFREYS_PRIOR:
					buffer.append("Jeffreys");
					break;
			}
			if (param.priorType != BeautiOptions.NONE && !param.isStatistic) {
				buffer.append(", initial=" + param.initial);
			}

			return buffer.toString();
		}

		public String getColumnName(int column) {
			return columnNames[column];
		}

		public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}

		public boolean isCellEditable(int row, int col) {
			if (col == 1) {
				return true;
			}
			return false;
		}

		public String toString() {
			StringBuffer buffer = new StringBuffer();

			buffer.append(getColumnName(0));
			for (int j = 1; j < getColumnCount(); j++) {
				buffer.append("\t");
				buffer.append(getColumnName(j));
			}
			buffer.append("\n");

			for (int i = 0; i < getRowCount(); i++) {
				buffer.append(getValueAt(i, 0));
				for (int j = 1; j < getColumnCount(); j++) {
					buffer.append("\t");
					buffer.append(getValueAt(i, j));
				}
				buffer.append("\n");
			}

			return buffer.toString();
		}
	};

	class DoubleRenderer extends TableRenderer {

		/**
		 *
		 */
		private static final long serialVersionUID = -2614341608257369805L;

		public DoubleRenderer(int alignment, Insets insets) {

			super(true, alignment, insets);
		}

		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
		                                               boolean hasFocus, int row, int column) {

			String s;
			if (((Double)value).isNaN()) {
				s = "random";
			} else {
				s = formatter.format(((Double)value).doubleValue());
			}
			return super.getTableCellRendererComponent(table, s, isSelected, hasFocus, row, column);

		}
	};

	public class ButtonRenderer extends JButton implements TableCellRenderer {

		/**
		 *
		 */
		private static final long serialVersionUID = -2416184092883649169L;

		public ButtonRenderer(int alignment, Insets insets) {
			setOpaque(true);
			setHorizontalAlignment(alignment);
			setMargin(insets);
		}

		public Component getTableCellRendererComponent(JTable table, Object value,
		                                               boolean isSelected, boolean hasFocus, int row, int column) {
			setEnabled(table.isEnabled());
			setFont(table.getFont());
			if (isSelected) {
				setForeground(table.getSelectionForeground());
				setBackground(table.getSelectionBackground());
			} else{
				setForeground(table.getForeground());
				setBackground(UIManager.getColor("Button.background"));
			}
			setText( (value ==null) ? "" : value.toString() );
			return this;
		}
	};

	public class ButtonEditor extends DefaultCellEditor {
		/**
		 *
		 */
		private static final long serialVersionUID = 6372738480075411674L;
		protected JButton button;
		private String label;
		private boolean isPushed;
		private int row;

		public ButtonEditor(int alignment, Insets insets) {
			super(new JCheckBox());
			button = new JButton();
			button.setOpaque(true);
			button.setHorizontalAlignment(alignment);
			button.setMargin(insets);
			button.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					fireEditingStopped();
				}
			});
		}

		public Component getTableCellEditorComponent(JTable table, Object value,
		                                             boolean isSelected, int row, int column) {
			button.setEnabled(table.isEnabled());
			button.setFont(table.getFont());
			if (isSelected) {
				button.setForeground(table.getSelectionForeground());
				button.setBackground(table.getSelectionBackground());
			} else{
				button.setForeground(table.getForeground());
				button.setBackground(table.getBackground());
			}
			label = (value ==null) ? "" : value.toString();
			button.setText( label );
			isPushed = true;
			this.row = row;
			return button;
		}

		public Object getCellEditorValue() {
			if (isPushed)  {
				priorButtonPressed(row);
			}
			isPushed = false;
			return new String( label ) ;
		}

		public boolean stopCellEditing() {
			isPushed = false;
			return super.stopCellEditing();
		}

		protected void fireEditingStopped() {
			super.fireEditingStopped();
		}
	}


}
