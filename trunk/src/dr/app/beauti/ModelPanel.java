/*
 * ModelPanel.java
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

import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: ModelPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class ModelPanel extends OptionsPanel implements Exportable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 2778103564318492601L;
	JComboBox nucSubstCombo = new JComboBox(new String[] {"HKY", "GTR"});
	JComboBox aaSubstCombo = new JComboBox(new String[] {"Blosum62", "Dayhoff", "JTT", "mtREV", "cpREV", "WAG"});
	JComboBox heteroCombo = new JComboBox(new String[] {"None", "Gamma", "Invariant Sites", "Gamma + Invariant Sites"});

	JComboBox gammaCatCombo = new JComboBox(new String[] {"4", "5", "6", "7", "8", "9", "10"});
	JLabel gammaCatLabel;

	JCheckBox codingCheck = new JCheckBox("Partition into codon positions");
	JCheckBox substUnlinkCheck = new JCheckBox("Unlink substitution model across codon positions");
	JCheckBox heteroUnlinkCheck = new JCheckBox("Unlink rate heterogeneity model across codon positions");
	JCheckBox freqsUnlinkCheck = new JCheckBox("Unlink base frequencies across codon positions");

	JCheckBox fixedSubstitutionRateCheck = new JCheckBox("Fix mean substitution rate");
	JLabel substitutionRateLabel = new JLabel("Mean substitution rate:");
	RealNumberField substitutionRateField = new RealNumberField(Double.MIN_VALUE, Double.POSITIVE_INFINITY);

	JComboBox clockModelCombo = new JComboBox(new String[] {
	    "Strict Clock",
	    "Relaxed Clock: Uncorrelated Exponential",
		"Relaxed Clock: Uncorrelated Lognormal" });

	BeautiFrame frame = null;

	boolean warningShown = false;
	boolean hasSetFixedSubstitutionRate = false;

	boolean settingOptions = false;

	boolean hasAlignment = false;
	boolean isNucleotides = true;

	public ModelPanel(BeautiFrame parent) {

		super(12, 18);

		this.frame = parent;

		setOpaque(false);

		substUnlinkCheck.setOpaque(false);
		substUnlinkCheck.setEnabled(false);
		heteroUnlinkCheck.setOpaque(false);
		heteroUnlinkCheck.setEnabled(false);
		freqsUnlinkCheck.setOpaque(false);
		freqsUnlinkCheck.setEnabled(false);

		java.awt.event.ItemListener listener = new java.awt.event.ItemListener() {
			public void itemStateChanged(java.awt.event.ItemEvent ev) {
				frame.modelChanged();
			}
		};

		nucSubstCombo.setOpaque(false);
		nucSubstCombo.addItemListener(listener);
		aaSubstCombo.setOpaque(false);
		aaSubstCombo.addItemListener(listener);

		heteroCombo.setOpaque(false);
		heteroCombo.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {

					frame.modelChanged();

					if (heteroCombo.getSelectedIndex() == 1 || heteroCombo.getSelectedIndex() == 3) {
						gammaCatLabel.setEnabled(true);
						gammaCatCombo.setEnabled(true);
					} else {
						gammaCatLabel.setEnabled(false);
						gammaCatCombo.setEnabled(false);
					}
				}
			}
		);

		gammaCatCombo.setOpaque(false);
		gammaCatCombo.addItemListener(listener);

		codingCheck.setOpaque(false);
		codingCheck.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {

					frame.modelChanged();

					if (codingCheck.isSelected()) {
						substUnlinkCheck.setEnabled(true);
						heteroUnlinkCheck.setEnabled(true);
						freqsUnlinkCheck.setEnabled(true);
					} else {
						substUnlinkCheck.setEnabled(false);
						substUnlinkCheck.setSelected(false);
						heteroUnlinkCheck.setEnabled(false);
						heteroUnlinkCheck.setSelected(false);
						freqsUnlinkCheck.setEnabled(false);
						freqsUnlinkCheck.setSelected(false);
					}
				}
			}
		);

		substUnlinkCheck.addItemListener(listener);
		heteroUnlinkCheck.addItemListener(listener);
		freqsUnlinkCheck.addItemListener(listener);

		fixedSubstitutionRateCheck.setOpaque(false);
		fixedSubstitutionRateCheck.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					boolean fixed = fixedSubstitutionRateCheck.isSelected();
					substitutionRateLabel.setEnabled(fixed);
					substitutionRateField.setEnabled(fixed);
					hasSetFixedSubstitutionRate = true;
					frame.modelChanged();
				}
			}
		);

		substitutionRateField.addKeyListener(new java.awt.event.KeyAdapter() {
			public void keyTyped(java.awt.event.KeyEvent ev) {
				frame.mcmcChanged();
			}});

        clockModelCombo.setOpaque(false);
        clockModelCombo.addItemListener(listener);

		setupPanel();
	}

	private void setupPanel() {

		removeAll();

		if (hasAlignment) {
			if (isNucleotides) {
				addComponentWithLabel("Substitution Model:", nucSubstCombo);
			} else {
				addComponentWithLabel("Substitution Model:", aaSubstCombo);
			}

			addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
			gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);

			if (isNucleotides) {
				addSeparator();

				addComponent(codingCheck);

				addComponent(substUnlinkCheck);
				addComponent(heteroUnlinkCheck);
		//		addComponent(freqsUnlinkCheck);
			}

			addSeparator();

			addComponent(fixedSubstitutionRateCheck);
			addComponents(substitutionRateLabel, substitutionRateField);
			substitutionRateField.setColumns(10);

			addSeparator();
        }

		addComponentWithLabel("Molecular Clock Model:", clockModelCombo);
        validate();
        repaint();
	}

	public void setOptions(BeautiOptions options) {

		settingOptions = true;

		if (options.alignment != null) {
			hasAlignment = true;
			if (options.alignment.getDataType() == dr.evolution.datatype.Nucleotides.INSTANCE) {

				if (options.nucSubstitutionModel == BeautiOptions.GTR) {
					nucSubstCombo.setSelectedIndex(1);
				} else {
					nucSubstCombo.setSelectedIndex(0);
				}

				isNucleotides = true;
			} else {

				aaSubstCombo.setSelectedIndex(options.aaSubstitutionModel);
				isNucleotides = false;
			}
		} else {
			hasAlignment = false;
		}

		if (options.gammaHetero && !options.invarHetero) {
			heteroCombo.setSelectedIndex(1);
		} else if (!options.gammaHetero && options.invarHetero) {
			heteroCombo.setSelectedIndex(2);
		} else if (options.gammaHetero && options.invarHetero) {
			heteroCombo.setSelectedIndex(3);
		} else {
			heteroCombo.setSelectedIndex(0);
		}

		gammaCatCombo.setSelectedIndex(options.gammaCategories - 4);

		codingCheck.setSelected(options.codonHetero);
		substUnlinkCheck.setSelected(options.unlinkedSubstitutionModel);
		heteroUnlinkCheck.setSelected(options.unlinkedHeterogeneityModel);
		freqsUnlinkCheck.setSelected(options.unlinkedFrequencyModel);

		hasSetFixedSubstitutionRate = options.hasSetFixedSubstitutionRate;
		if (!hasSetFixedSubstitutionRate) {
			if (options.maximumTipHeight > 0.0) {
				options.meanSubstitutionRate = 0.001;
				options.fixedSubstitutionRate = false;
			} else {
				options.meanSubstitutionRate = 1.0;
				options.fixedSubstitutionRate = true;
			}
		}

		fixedSubstitutionRateCheck.setSelected(options.fixedSubstitutionRate);
		substitutionRateField.setValue(options.meanSubstitutionRate);
		substitutionRateField.setEnabled(options.fixedSubstitutionRate);

		clockModelCombo.setSelectedIndex(options.clockModel);

		setupPanel();

		settingOptions = false;

		validate();
		repaint();
	}

	public void getOptions(BeautiOptions options) {

		// This prevents options be overwritten due to listeners calling
		// this function (indirectly through modelChanged()) whilst in the
		// middle of the setOptions() method.
		if (settingOptions) return;

		if (nucSubstCombo.getSelectedIndex() == 1) {
			options.nucSubstitutionModel = BeautiOptions.GTR;
		} else {
			options.nucSubstitutionModel = BeautiOptions.HKY;
		}
		options.aaSubstitutionModel = aaSubstCombo.getSelectedIndex();

		if (heteroCombo.getSelectedIndex() == 1 || heteroCombo.getSelectedIndex() == 3) {
			options.gammaHetero = true;
		} else {
			options.gammaHetero = false;
		}

		if (heteroCombo.getSelectedIndex() == 2 || heteroCombo.getSelectedIndex() == 3) {
			options.invarHetero = true;
		} else {
			options.invarHetero = false;
		}

		options.gammaCategories = gammaCatCombo.getSelectedIndex() + 4;

		options.codonHetero = codingCheck.isSelected();
		options.unlinkedSubstitutionModel = substUnlinkCheck.isSelected();
		options.unlinkedHeterogeneityModel = heteroUnlinkCheck.isSelected();
		options.unlinkedFrequencyModel = freqsUnlinkCheck.isSelected();

		options.hasSetFixedSubstitutionRate = hasSetFixedSubstitutionRate;
		options.fixedSubstitutionRate = fixedSubstitutionRateCheck.isSelected();
		options.meanSubstitutionRate = substitutionRateField.getValue().doubleValue();

		boolean fixed = fixedSubstitutionRateCheck.isSelected();
		if (!warningShown && !fixed && options.maximumTipHeight == 0.0) {
			JOptionPane.showMessageDialog(frame,
			        "You have chosen to sample substitution rates but all \n"+
			        "the sequences have the same date. In order for this to \n"+
			        "work, a strong prior is required on the substitution\n"+
			        "rate or the root of the tree.",
			        "Warning",
			        JOptionPane.WARNING_MESSAGE);
			warningShown = true;
		}

		options.clockModel = clockModelCombo.getSelectedIndex();
	}

    public JComponent getExportableComponent() {
		return this;
	}
}
