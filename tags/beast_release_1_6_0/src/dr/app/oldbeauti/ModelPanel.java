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

package dr.app.oldbeauti;

import dr.app.gui.components.RealNumberField;
import dr.evolution.datatype.DataType;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

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
    JComboBox binarySubstCombo = new JComboBox(new String[] {"Simple", "Covarion"});

    JComboBox frequencyCombo =new JComboBox(new String[] {"Estimated", "Empirical", "All equal"});

    JComboBox heteroCombo = new JComboBox(new String[] {"None", "Gamma", "Invariant Sites", "Gamma + Invariant Sites"});

    JComboBox gammaCatCombo = new JComboBox(new String[] {"4", "5", "6", "7", "8", "9", "10"});
    JLabel gammaCatLabel;

    JComboBox codingCombo = new JComboBox(new String[] {
            "Off",
            "2 partitions: codon positions (1 + 2), 3",
            "3 partitions: codon positions 1, 2, 3"});
    JCheckBox substUnlinkCheck = new JCheckBox("Unlink substitution model across codon positions");
    JCheckBox heteroUnlinkCheck = new JCheckBox("Unlink rate heterogeneity model across codon positions");
    JCheckBox freqsUnlinkCheck = new JCheckBox("Unlink base frequencies across codon positions");

    JButton setSRD06Button;

    JCheckBox fixedSubstitutionRateCheck = new JCheckBox("Fix mean substitution rate:");
    JLabel substitutionRateLabel = new JLabel("Mean substitution rate:");
    RealNumberField substitutionRateField = new RealNumberField(Double.MIN_VALUE, Double.POSITIVE_INFINITY);

    JComboBox clockModelCombo = new JComboBox(new String[] {
            "Strict Clock",
            "Random Local Clock",
            "Relaxed Clock: Uncorrelated Lognormal",
            "Relaxed Clock: Uncorrelated Exponential" } );

    BeautiFrame frame = null;

    boolean warningShown = false;
    boolean hasSetFixedSubstitutionRate = false;

    boolean settingOptions = false;

    boolean hasAlignment = false;

    int dataType = DataType.NUCLEOTIDES;

    public ModelPanel(BeautiFrame parent) {

        super(12, 18);

        this.frame = parent;

        setOpaque(false);

        setupComponent(substUnlinkCheck);
        substUnlinkCheck.setEnabled(false);
        substUnlinkCheck.setToolTipText("" +
                "<html>Gives each codon position partition different<br>" +
                "substitution model parameters.</html>");

        setupComponent(heteroUnlinkCheck);
        heteroUnlinkCheck.setEnabled(false);
        heteroUnlinkCheck.setToolTipText("<html>Gives each codon position partition different<br>rate heterogeneity model parameters.</html>");

        setupComponent(freqsUnlinkCheck);
        freqsUnlinkCheck.setEnabled(false);
        freqsUnlinkCheck.setToolTipText("<html>Gives each codon position partition different<br>nucleotide frequency parameters.</html>");

        java.awt.event.ItemListener listener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                frame.modelChanged();
            }
        };

        setupComponent(nucSubstCombo);
        nucSubstCombo.addItemListener(listener);
        nucSubstCombo.setToolTipText("<html>Select the type of nucleotide substitution model.</html>");

        setupComponent(aaSubstCombo);
        aaSubstCombo.addItemListener(listener);
        aaSubstCombo.setToolTipText("<html>Select the type of amino acid substitution model.</html>");

        setupComponent(binarySubstCombo);
        binarySubstCombo.addItemListener(listener);
        binarySubstCombo.setToolTipText("<html>Select the type of binay substitution model.</html>");

        setupComponent(frequencyCombo);
        frequencyCombo.addItemListener(listener);
        frequencyCombo.setToolTipText("<html>Select the policy for determining the base frequencies.</html>");

        setupComponent(heteroCombo);
        heteroCombo.setToolTipText("<html>Select the type of site-specific rate<br>heterogeneity model.</html>");
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

        setupComponent(gammaCatCombo);
        gammaCatCombo.setToolTipText("<html>Select the number of categories to use for<br>the discrete gamma rate heterogeneity model.</html>");
        gammaCatCombo.addItemListener(listener);

        setupComponent(codingCombo);
        codingCombo.setToolTipText("<html>Select how to partition the codon positions.</html>");
        codingCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {

                        frame.modelChanged();

                        if (codingCombo.getSelectedIndex() != 0) { // codon position partitioning
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

        setSRD06Button = new JButton(setSRD06Action);
        setupComponent(setSRD06Button);
        setSRD06Button.setToolTipText("<html>Sets the SRD06 model as described in<br>" +
                "Shapiro, Rambaut & Drummond (2006) <i>MBE</i> <b>23</b>: 7-9.</html>");

        setupComponent(fixedSubstitutionRateCheck);
        fixedSubstitutionRateCheck.setToolTipText(
                "<html>Select this option to fix the substitution rate<br>" +
                        "rather than try to infer it. If this option is<br>" +
                        "turned off then either the sequences should have<br>" +
                        "dates or the tree should have sufficient calibration<br>" +
                        "informations specified as priors.</html>");
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

        setupComponent(substitutionRateField);
        substitutionRateField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
                frame.mcmcChanged();
            }});
        substitutionRateField.setToolTipText("<html>Enter the substitution rate here.</html>");

        setupComponent(clockModelCombo);
        clockModelCombo.setToolTipText("<html>Select either a strict molecular clock or<br>or a relaxed clock model.</html>");
        clockModelCombo.addItemListener(listener);

        setupPanel();
    }

    private void setupComponent(JComponent comp) {
        comp.setOpaque(false);

        //comp.setFont(UIManager.getFont("SmallSystemFont"));
        //comp.putClientProperty("JComponent.sizeVariant", "small");
        if (comp instanceof JButton) {
            comp.putClientProperty("JButton.buttonType", "roundRect");
        }
        if (comp instanceof JComboBox) {
            comp.putClientProperty("JComboBox.isSquare", Boolean.TRUE);
        }
    }

    private void setupPanel() {

        removeAll();

        if (hasAlignment) {

            switch (dataType){
                case DataType.NUCLEOTIDES:
                    addComponentWithLabel("Substitution Model:", nucSubstCombo, true);
                    addComponentWithLabel("Base frequencies:", frequencyCombo);
                    addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                    gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);

                    addSeparator();

                    JPanel panel = new JPanel(new BorderLayout(6,6));
                    panel.setOpaque(false);
                    panel.add(codingCombo, BorderLayout.CENTER);
                    panel.add(setSRD06Button, BorderLayout.EAST);
                    addComponentWithLabel("Partition into codon positions:", panel);

                    panel = new JPanel();
                    panel.setOpaque(false);
                    panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));
                    panel.setBorder(BorderFactory.createTitledBorder("Link/Unlink parameters:"));
                    panel.add(substUnlinkCheck);
                    panel.add(heteroUnlinkCheck);
                    panel.add(freqsUnlinkCheck);

                    addComponent(panel);
                    break;

                case DataType.AMINO_ACIDS:
                    addComponentWithLabel("Substitution Model:", aaSubstCombo);
                    addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                    gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);

                    break;

                case DataType.TWO_STATES:
                case DataType.COVARION:
                    addComponentWithLabel("Substitution Model:", binarySubstCombo);
                    addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                    gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);

                    break;

                default:
                    throw new IllegalArgumentException("Unknown data type");

            }

            addSeparator();

            //addComponent(fixedSubstitutionRateCheck);
            substitutionRateField.setColumns(10);
            addComponents(fixedSubstitutionRateCheck, substitutionRateField);

            addSeparator();
        }

        addComponentWithLabel("Molecular Clock Model:", clockModelCombo);
        validate();
        repaint();
    }

    private void setSRD06Model() {
        nucSubstCombo.setSelectedIndex(0);
        heteroCombo.setSelectedIndex(1);
        codingCombo.setSelectedIndex(1);
        substUnlinkCheck.setSelected(true);
        heteroUnlinkCheck.setSelected(true);
    }

    public void setOptions(BeautiOptions options) {

        settingOptions = true;

        if (options.alignment != null) {
            hasAlignment = true;

            dataType=options.dataType;
            switch(dataType){
                case DataType.NUCLEOTIDES:
                    if (options.nucSubstitutionModel == BeautiOptions.GTR) {
                        nucSubstCombo.setSelectedIndex(1);
                    } else {
                        nucSubstCombo.setSelectedIndex(0);
                    }

                    frequencyCombo.setSelectedIndex(options.frequencyPolicy);

                    break;

                case DataType.AMINO_ACIDS:
                    aaSubstCombo.setSelectedIndex(options.aaSubstitutionModel);
                    break;

                case DataType.TWO_STATES:
                case DataType.COVARION:
                    binarySubstCombo.setSelectedIndex(options.binarySubstitutionModel);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown data type");
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

        if (options.codonHeteroPattern == null) {
            codingCombo.setSelectedIndex(0);
        } else if (options.codonHeteroPattern.equals("112")) {
            codingCombo.setSelectedIndex(1);
        } else {
            codingCombo.setSelectedIndex(2);
        }

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

        switch (options.clockModel) {
            case BeautiOptions.STRICT_CLOCK:
                clockModelCombo.setSelectedIndex(0); break;
            case BeautiOptions.RANDOM_LOCAL_CLOCK:
                clockModelCombo.setSelectedIndex(1); break;
            case BeautiOptions.UNCORRELATED_LOGNORMAL:
                clockModelCombo.setSelectedIndex(2); break;
            case BeautiOptions.UNCORRELATED_EXPONENTIAL:
                clockModelCombo.setSelectedIndex(3); break;
            default:
                throw new IllegalArgumentException("Unknown option for clock model");
        }
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

        options.binarySubstitutionModel = binarySubstCombo.getSelectedIndex();

        options.frequencyPolicy = frequencyCombo.getSelectedIndex();

        options.gammaHetero = heteroCombo.getSelectedIndex() == 1 || heteroCombo.getSelectedIndex() == 3;

        options.invarHetero = heteroCombo.getSelectedIndex() == 2 || heteroCombo.getSelectedIndex() == 3;

        options.gammaCategories = gammaCatCombo.getSelectedIndex() + 4;

        if (codingCombo.getSelectedIndex() == 0) {
            options.codonHeteroPattern = null;
        } else if (codingCombo.getSelectedIndex() == 1) {
            options.codonHeteroPattern = "112";
        } else {
            options.codonHeteroPattern = "123";
        }

        options.unlinkedSubstitutionModel = substUnlinkCheck.isSelected();
        options.unlinkedHeterogeneityModel = heteroUnlinkCheck.isSelected();
        options.unlinkedFrequencyModel = freqsUnlinkCheck.isSelected();

        options.hasSetFixedSubstitutionRate = hasSetFixedSubstitutionRate;
        options.fixedSubstitutionRate = fixedSubstitutionRateCheck.isSelected();
        options.meanSubstitutionRate = substitutionRateField.getValue();

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

        switch (clockModelCombo.getSelectedIndex()) {
            case 0:
                options.clockModel = BeautiOptions.STRICT_CLOCK; break;
            case 1:
                options.clockModel = BeautiOptions.RANDOM_LOCAL_CLOCK; break;
            case 2:
                options.clockModel = BeautiOptions.UNCORRELATED_LOGNORMAL; break;
            case 3:
                options.clockModel = BeautiOptions.UNCORRELATED_EXPONENTIAL; break;
            default:
                throw new IllegalArgumentException("Unknown option for clock model");
        }
    }

    public JComponent getExportableComponent() {

        return this;
    }

    private Action setSRD06Action = new AbstractAction("Use SRD06 Model") {
        public void actionPerformed(ActionEvent actionEvent) {
            setSRD06Model();
        }
    };

}
