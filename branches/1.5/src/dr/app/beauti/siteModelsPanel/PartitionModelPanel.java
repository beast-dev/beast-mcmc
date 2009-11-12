/*
 * PartitionModelPanel.java
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

package dr.app.beauti.siteModelsPanel;

import dr.app.beauti.BeautiApp;
import dr.evomodel.substmodel.AminoAcidModelType;
import dr.app.beauti.enumTypes.BinaryModelType;
import dr.app.beauti.enumTypes.FrequencyPolicyType;
import dr.evomodel.substmodel.NucModelType;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.util.PanelUtils;
import dr.evolution.datatype.DataType;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;
import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionModelPanel extends OptionsPanel {

    // Components
	private static final long serialVersionUID = -1645661616353099424L;
	
	private JComboBox nucSubstCombo = new JComboBox(EnumSet.range(NucModelType.HKY, NucModelType.TN93).toArray());
    private JComboBox aaSubstCombo = new JComboBox(AminoAcidModelType.values());
    private JComboBox binarySubstCombo = new JComboBox(BinaryModelType.values());
    private JCheckBox useAmbiguitiesTreeLikelihoodCheck
            = new JCheckBox("Use ambiguities in those treeLikelihood associating this substitution model");

    private JComboBox frequencyCombo = new JComboBox(FrequencyPolicyType.values());

    private JComboBox heteroCombo = new JComboBox(
            new String[]{"None", "Gamma", "Invariant Sites", "Gamma + Invariant Sites"});

    private JComboBox gammaCatCombo = new JComboBox(new String[]{"4", "5", "6", "7", "8", "9", "10"});
    private JLabel gammaCatLabel;

    private JComboBox codingCombo = new JComboBox(new String[]{
            "Off",
            "2 partitions: codon positions (1 + 2), 3",
            "3 partitions: codon positions 1, 2, 3"});

    private JCheckBox substUnlinkCheck = new JCheckBox("Unlink substitution rate parameters across codon positions");
    private JCheckBox heteroUnlinkCheck =
            new JCheckBox("Unlink rate heterogeneity model across codon positions");
    private JCheckBox freqsUnlinkCheck = new JCheckBox("Unlink base frequencies across codon positions");

    private JButton setSRD06Button;

    private JCheckBox dolloCheck = new JCheckBox("Use Stochastic Dollo Model");
    // private JComboBox dolloCombo = new JComboBox(new String[]{"Analytical", "Sample"});

    PartitionSubstitutionModel model;

    public PartitionModelPanel(PartitionSubstitutionModel partitionModel) {

        super(12, 30);

        this.model = partitionModel;

        initCodonPartitionComponents();

        PanelUtils.setupComponent(nucSubstCombo);
        nucSubstCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setNucSubstitutionModel( (NucModelType) nucSubstCombo.getSelectedItem());                
            }
        });
        nucSubstCombo.setToolTipText("<html>Select the type of nucleotide substitution model.</html>");

        PanelUtils.setupComponent(aaSubstCombo);
        aaSubstCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setAaSubstitutionModel((AminoAcidModelType) aaSubstCombo.getSelectedItem());
            }
        });
        aaSubstCombo.setToolTipText("<html>Select the type of amino acid substitution model.</html>");

        PanelUtils.setupComponent(binarySubstCombo);
        binarySubstCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setBinarySubstitutionModel((BinaryModelType) binarySubstCombo.getSelectedItem());
                useAmbiguitiesTreeLikelihoodCheck.setSelected(binarySubstCombo.getSelectedItem() == BinaryModelType.BIN_COVARION);
                useAmbiguitiesTreeLikelihoodCheck.setEnabled(binarySubstCombo.getSelectedItem() != BinaryModelType.BIN_COVARION);
            }
        });
        binarySubstCombo.setToolTipText("<html>Select the type of binary substitution model.</html>");

        PanelUtils.setupComponent(useAmbiguitiesTreeLikelihoodCheck);
        useAmbiguitiesTreeLikelihoodCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
            	model.setUseAmbiguitiesTreeLikelihood(useAmbiguitiesTreeLikelihoodCheck.isSelected());
            }
        });
        useAmbiguitiesTreeLikelihoodCheck.setToolTipText("<html>Detemine useAmbiguities in &lt treeLikelihood &gt .</html>");

        PanelUtils.setupComponent(frequencyCombo);
        frequencyCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setFrequencyPolicy((FrequencyPolicyType) frequencyCombo.getSelectedItem());
            }
        });
        frequencyCombo.setToolTipText("<html>Select the policy for determining the base frequencies.</html>");

        PanelUtils.setupComponent(heteroCombo);
        heteroCombo.setToolTipText("<html>Select the type of site-specific rate<br>heterogeneity model.</html>");
        heteroCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {

                        boolean gammaHetero = heteroCombo.getSelectedIndex() == 1 || heteroCombo.getSelectedIndex() == 3;

                        model.setGammaHetero(gammaHetero);
                        model.setInvarHetero(heteroCombo.getSelectedIndex() == 2 || heteroCombo.getSelectedIndex() == 3);

                        if (gammaHetero) {
                            gammaCatLabel.setEnabled(true);
                            gammaCatCombo.setEnabled(true);
                        } else {
                            gammaCatLabel.setEnabled(false);
                            gammaCatCombo.setEnabled(false);
                        }                        

                        if (codingCombo.getSelectedIndex() != 0) {
                            heteroUnlinkCheck.setEnabled(heteroCombo.getSelectedIndex() != 0);
                            heteroUnlinkCheck.setSelected(heteroCombo.getSelectedIndex() != 0);
                        }
                    }
                }
        );

        PanelUtils.setupComponent(gammaCatCombo);
        gammaCatCombo.setToolTipText("<html>Select the number of categories to use for<br>the discrete gamma rate heterogeneity model.</html>");
        gammaCatCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {

                model.setGammaCategories(gammaCatCombo.getSelectedIndex() + 4);
            }
        });

        Action setSRD06Action = new AbstractAction("Use SRD06 Model") {
            public void actionPerformed(ActionEvent actionEvent) {
                setSRD06Model();
            }
        };
        setSRD06Button = new JButton(setSRD06Action);
        PanelUtils.setupComponent(setSRD06Button);
        setSRD06Button.setToolTipText("<html>Sets the SRD06 model as described in<br>" +
                "Shapiro, Rambaut & Drummond (2006) <i>MBE</i> <b>23</b>: 7-9.</html>");


        PanelUtils.setupComponent(dolloCheck);
        dolloCheck.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                model.setDolloModel(true);
            }
        });
        dolloCheck.setEnabled(true);
        dolloCheck.setToolTipText("<html>Activates a stochastic dollo model - Alekseyenko, Lee & Suchard(2008) Syst Biol 57: 772-784.</html>");

        setupPanel();
        setOpaque(false);
    }

    /**
     * Sets the components up according to the partition model - but does not layout the top
     * level options panel.
     */
    private void setOptions() {

        if (SiteModelsPanel.DEBUG) {
            String modelName = (model == null) ? "null" : model.getName();
            Logger.getLogger("dr.app.beauti").info("ModelsPanel.setModelOptions(" + modelName + ")");
        }

        if (model == null) {
            return;
        }

        int dataType = model.getDataType().getType();
        switch (dataType) {
            case DataType.NUCLEOTIDES:                
                nucSubstCombo.setSelectedItem(model.getNucSubstitutionModel());  
                frequencyCombo.setSelectedItem(model.getFrequencyPolicy());

                break;

            case DataType.AMINO_ACIDS:
                aaSubstCombo.setSelectedItem(model.getAaSubstitutionModel());
                
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                binarySubstCombo.setSelectedItem(model.getBinarySubstitutionModel());
                useAmbiguitiesTreeLikelihoodCheck.setSelected(model.isUseAmbiguitiesTreeLikelihood());
                
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");
        }

        if (model.isGammaHetero() && !model.isInvarHetero()) {
            heteroCombo.setSelectedIndex(1);
        } else if (!model.isGammaHetero() && model.isInvarHetero()) {
            heteroCombo.setSelectedIndex(2);
        } else if (model.isGammaHetero() && model.isInvarHetero()) {
            heteroCombo.setSelectedIndex(3);
        } else {
            heteroCombo.setSelectedIndex(0);
        }

        gammaCatCombo.setSelectedIndex(model.getGammaCategories() - 4);

        if (model.getCodonHeteroPattern() == null) {
            codingCombo.setSelectedIndex(0);
        } else if (model.getCodonHeteroPattern().equals("112")) {
            codingCombo.setSelectedIndex(1);
        } else {
            codingCombo.setSelectedIndex(2);
        }

        substUnlinkCheck.setSelected(model.isUnlinkedSubstitutionModel());
        heteroUnlinkCheck.setSelected(model.isUnlinkedHeterogeneityModel());
        freqsUnlinkCheck.setSelected(model.isUnlinkedFrequencyModel());

        dolloCheck.setSelected(model.isDolloModel());
    }


    /**
     * Configure this panel for the Shapiro, Rambaut and Drummond 2006 codon position model
     */
    private void setSRD06Model() {
        nucSubstCombo.setSelectedIndex(0);
        heteroCombo.setSelectedIndex(1);
        codingCombo.setSelectedIndex(1);
        substUnlinkCheck.setSelected(true);
        heteroUnlinkCheck.setSelected(true);
    }

    /**
     * Lays out the appropriate components in the panel for this partition model.
     */
    private void setupPanel() {

        switch (model.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                addComponentWithLabel("Substitution Model:", nucSubstCombo);
                addComponentWithLabel("Base frequencies:", frequencyCombo);
                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                heteroCombo.setSelectedIndex(0);
                gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);
                gammaCatCombo.setEnabled(false);

                addSeparator();

                JPanel panel1 = new JPanel(new BorderLayout(6, 6));
                panel1.setOpaque(false);
                panel1.add(codingCombo, BorderLayout.CENTER);
                panel1.add(setSRD06Button, BorderLayout.EAST);
                addComponentWithLabel("Partition into codon positions:", panel1);

                JPanel panel2 = new JPanel();
                panel2.setOpaque(false);
                panel2.setLayout(new BoxLayout(panel2, BoxLayout.PAGE_AXIS));
                panel2.setBorder(BorderFactory.createTitledBorder("Link/Unlink parameters:"));
                panel2.add(substUnlinkCheck);
                panel2.add(heteroUnlinkCheck);
                panel2.add(freqsUnlinkCheck);

                addComponent(panel2);
                break;

            case DataType.AMINO_ACIDS:
                addComponentWithLabel("Substitution Model:", aaSubstCombo);
                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                heteroCombo.setSelectedIndex(0);
                gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);
                gammaCatCombo.setEnabled(false);

                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                addComponentWithLabel("Substitution Model:", binarySubstCombo);
                addComponentWithLabel("Base frequencies:", frequencyCombo);
                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                heteroCombo.setSelectedIndex(0);
                gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);
                gammaCatCombo.setEnabled(false);

                addSeparator();
                
                addComponentWithLabel("", useAmbiguitiesTreeLikelihoodCheck);

                break;

            default:
                throw new IllegalArgumentException("Unknown data type");

        }

        if (BeautiApp.advanced) {
            addSeparator();
            addComponent(dolloCheck);
        }

        setOptions();
    }

    /**
     * Initializes and binds the components related to modeling codon positions.
     */
    private void initCodonPartitionComponents() {

        PanelUtils.setupComponent(substUnlinkCheck);

        substUnlinkCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setUnlinkedSubstitutionModel(substUnlinkCheck.isSelected());
            }
        });
        substUnlinkCheck.setEnabled(false);
        substUnlinkCheck.setToolTipText("" +
                "<html>Gives each codon position partition different<br>" +
                "substitution model parameters.</html>");

        PanelUtils.setupComponent(heteroUnlinkCheck);
        heteroUnlinkCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setUnlinkedHeterogeneityModel(heteroUnlinkCheck.isSelected());
            }
        });
        heteroUnlinkCheck.setEnabled(heteroCombo.getSelectedIndex() != 0);
        heteroUnlinkCheck.setToolTipText("<html>Gives each codon position partition different<br>rate heterogeneity model parameters.</html>");

        PanelUtils.setupComponent(freqsUnlinkCheck);
        freqsUnlinkCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setUnlinkedFrequencyModel(freqsUnlinkCheck.isSelected());
            }
        });
        freqsUnlinkCheck.setEnabled(false);
        freqsUnlinkCheck.setToolTipText("<html>Gives each codon position partition different<br>nucleotide frequency parameters.</html>");

        PanelUtils.setupComponent(codingCombo);
        codingCombo.setToolTipText("<html>Select how to partition the codon positions.</html>");
        codingCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {

                switch (codingCombo.getSelectedIndex()) {
                    case 0:
                        model.setCodonHeteroPattern(null);
                        break;
                    case 1:
                        model.setCodonHeteroPattern("112");
                        break;
                    default:
                        model.setCodonHeteroPattern("123");
                        break;

                }

                if (codingCombo.getSelectedIndex() != 0) { // codon position partitioning
                    substUnlinkCheck.setEnabled(true);
                    heteroUnlinkCheck.setEnabled(heteroCombo.getSelectedIndex() != 0);
                    freqsUnlinkCheck.setEnabled(true);
                    substUnlinkCheck.setSelected(true);
                    heteroUnlinkCheck.setSelected(heteroCombo.getSelectedIndex() != 0);
                    freqsUnlinkCheck.setSelected(true);
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
    }
   
}
