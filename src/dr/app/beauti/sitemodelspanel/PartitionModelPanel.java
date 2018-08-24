/*
 * PartitionModelPanel.java
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

package dr.app.beauti.sitemodelspanel;

import dr.app.beauti.BeautiFrame;
import dr.evomodel.substmodel.aminoacid.AminoAcidModelType;
import dr.evomodel.substmodel.nucleotide.NucModelType;
import dr.app.beauti.components.continuous.ContinuousComponentOptions;
import dr.app.beauti.components.continuous.ContinuousSubstModelType;
import dr.app.beauti.components.discrete.DiscreteSubstModelType;
import dr.app.beauti.components.dollo.DolloComponentOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.BinaryModelType;
import dr.app.beauti.types.FrequencyPolicyType;
import dr.app.beauti.types.MicroSatModelType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.event.*;
import java.util.EnumSet;
import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PartitionModelPanel extends OptionsPanel {

    private static final boolean ENABLE_STOCHASTIC_DOLLO = true;

    // Components
    private static final long serialVersionUID = -1645661616353099424L;

    private JComboBox nucSubstCombo = new JComboBox(EnumSet.range(
            NucModelType.JC, NucModelType.TN93).toArray());
    private JComboBox aaSubstCombo = new JComboBox(AminoAcidModelType.values());
    private JComboBox binarySubstCombo = new JComboBox(
            new BinaryModelType[] { BinaryModelType.BIN_SIMPLE, BinaryModelType.BIN_COVARION });
    private JCheckBox useAmbiguitiesTreeLikelihoodCheck = new JCheckBox(
            "Use ambiguities in the tree likelihood associated with this model");

    private JComboBox frequencyCombo = new JComboBox(FrequencyPolicyType
            .values());

    private JComboBox heteroCombo = new JComboBox(new String[] { "None",
            "Gamma", "Invariant Sites", "Gamma + Invariant Sites" });

    private JComboBox gammaCatCombo = new JComboBox(new String[] { "4", "5",
            "6", "7", "8", "9", "10" });
    private JLabel gammaCatLabel;

    private JComboBox codingCombo = new JComboBox(new String[] { "Off",
            "2 partitions: positions (1 + 2), 3",
            "3 partitions: positions 1, 2, 3" });

    private JCheckBox substUnlinkCheck = new JCheckBox(
            "Unlink substitution rate parameters across codon positions");
    private JCheckBox heteroUnlinkCheck = new JCheckBox(
            "Unlink rate heterogeneity model across codon positions");
    private JCheckBox freqsUnlinkCheck = new JCheckBox(
            "Unlink base frequencies across codon positions");

    private JButton setYang96Button;
    private JButton setSRD06Button;

    private JCheckBox dolloCheck = new JCheckBox("Use stochastic Dollo model");
    // private JComboBox dolloCombo = new JComboBox(new String[]{"Analytical",
    // "Sample"});

    private JComboBox discreteTraitSiteModelCombo = new JComboBox(
            DiscreteSubstModelType.values());
    private JCheckBox activateBSSVS = new JCheckBox(
            // "Activate BSSVS"
            "Infer social network with BSSVS");
    private JButton setupGLMButton;
    private GLMSettingsDialog glmSettingsDialog = null;

    private JComboBox continuousTraitSiteModelCombo = new JComboBox(
            ContinuousSubstModelType.values());

    private JCheckBox latLongCheck = new JCheckBox(
            "Bivariate trait represents latitude and longitude");

    private JCheckBox useLambdaCheck = new JCheckBox(
            "Estimate phylogenetic signal using tree transform");

    private JCheckBox addJitterCheck = new JCheckBox(
            "Add random jitter to tips");
    private JLabel jitterWindowLabel = new JLabel("Jitter window size:");
    private RealNumberField jitterWindowText = new RealNumberField(0, Double.POSITIVE_INFINITY);

    private JTextArea citationText;

    // =========== micro sat ===========
    private JTextField microsatName = new JTextField();
    private WholeNumberField microsatMax = new WholeNumberField(2,
            Integer.MAX_VALUE);
    private WholeNumberField microsatMin = new WholeNumberField(1,
            Integer.MAX_VALUE);
    private JComboBox rateProportionCombo = new JComboBox(
            MicroSatModelType.RateProportionality.values());
    private JComboBox mutationBiasCombo = new JComboBox(
            MicroSatModelType.MutationalBias.values());
    private JComboBox phaseCombo = new JComboBox(MicroSatModelType.Phase
            .values());
    JCheckBox shareMicroSatCheck = new JCheckBox(
            "Share one microsatellite among all substitution model(s)");

    protected final PartitionSubstitutionModel model;

    final BeautiFrame frame;

    public PartitionModelPanel(final BeautiFrame frame, final PartitionSubstitutionModel partitionModel) {

        super(12, (OSType.isMac() ? 6 : 24));

        this.frame = frame;

        this.model = partitionModel;

        initCodonPartitionComponents();

        PanelUtils.setupComponent(nucSubstCombo);
        nucSubstCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setNucSubstitutionModel((NucModelType) nucSubstCombo
                        .getSelectedItem());
                if (model.getNucSubstitutionModel() == NucModelType.JC) {
                    frequencyCombo.getSelectedItem();

                    frequencyCombo.setSelectedItem(FrequencyPolicyType.ALLEQUAL);
                }
                frequencyCombo.setEnabled(model.getNucSubstitutionModel() != NucModelType.JC);
            }
        });
        nucSubstCombo
                .setToolTipText("<html>Select the type of nucleotide substitution model.</html>");

        PanelUtils.setupComponent(aaSubstCombo);
        aaSubstCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                AminoAcidModelType type = (AminoAcidModelType) aaSubstCombo.getSelectedItem();
                model.setAaSubstitutionModel(type);
                citationText.setText(type.getCitation().toString());
            }
        });
        aaSubstCombo
                .setToolTipText("<html>Select the type of amino acid substitution model.</html>");

        PanelUtils.setupComponent(binarySubstCombo);
        binarySubstCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model
                        .setBinarySubstitutionModel((BinaryModelType) binarySubstCombo
                                .getSelectedItem());
                useAmbiguitiesTreeLikelihoodCheck.setSelected(binarySubstCombo
                        .getSelectedItem() == BinaryModelType.BIN_COVARION);
                useAmbiguitiesTreeLikelihoodCheck.setEnabled(binarySubstCombo
                        .getSelectedItem() != BinaryModelType.BIN_COVARION);
            }
        });
        binarySubstCombo
                .setToolTipText("<html>Select the type of binary substitution model.</html>");

        PanelUtils.setupComponent(useAmbiguitiesTreeLikelihoodCheck);
        useAmbiguitiesTreeLikelihoodCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model
                        .setUseAmbiguitiesTreeLikelihood(useAmbiguitiesTreeLikelihoodCheck
                                .isSelected());
            }
        });
        useAmbiguitiesTreeLikelihoodCheck
                .setToolTipText("<html>Detemine useAmbiguities in &lt treeLikelihood &gt .</html>");

        PanelUtils.setupComponent(frequencyCombo);
        frequencyCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setFrequencyPolicy((FrequencyPolicyType) frequencyCombo
                        .getSelectedItem());
            }
        });
        frequencyCombo
                .setToolTipText("<html>Select the policy for determining the base frequencies.</html>");

        PanelUtils.setupComponent(heteroCombo);
        heteroCombo
                .setToolTipText("<html>Select the type of site-specific rate<br>heterogeneity model.</html>");
        heteroCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {

                boolean gammaHetero = heteroCombo.getSelectedIndex() == 1
                        || heteroCombo.getSelectedIndex() == 3;

                model.setGammaHetero(gammaHetero);
                model.setInvarHetero(heteroCombo.getSelectedIndex() == 2
                        || heteroCombo.getSelectedIndex() == 3);

                if (gammaHetero) {
                    gammaCatLabel.setEnabled(true);
                    gammaCatCombo.setEnabled(true);
                } else {
                    gammaCatLabel.setEnabled(false);
                    gammaCatCombo.setEnabled(false);
                }

                if (codingCombo.getSelectedIndex() != 0) {
                    heteroUnlinkCheck
                            .setEnabled(heteroCombo.getSelectedIndex() != 0);
                    heteroUnlinkCheck.setSelected(heteroCombo
                            .getSelectedIndex() != 0);
                }
            }
        });

        PanelUtils.setupComponent(gammaCatCombo);
        gammaCatCombo
                .setToolTipText("<html>Select the number of categories to use for<br>the discrete gamma rate heterogeneity model.</html>");
        gammaCatCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {

                model.setGammaCategories(gammaCatCombo.getSelectedIndex() + 4);
            }
        });

        setYang96Button = new JButton("Use Yang96 model");
        setYang96Button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                setYang96Model();
            }
        });
        PanelUtils.setupComponent(setYang96Button);
        setYang96Button
                .setToolTipText("<html>Sets a 3 codon-position model with independent GTR and Gamma as described in<br>"
                        + "Yang (1996) <i>J Mol Evol</i> <b>42</b>: 587-596. This model is named 3' in this paper.</html>");

        setSRD06Button = new JButton("Use SRD06 model");
        setSRD06Button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                setSRD06Model();
            }
        });
        PanelUtils.setupComponent(setSRD06Button);
        setSRD06Button
                .setToolTipText("<html>Sets the SRD06 model as described in<br>"
                        + "Shapiro, Rambaut & Drummond (2006) <i>MBE</i> <b>23</b>: 7-9.</html>");

        citationText = new JTextArea(1, 40);
        citationText.setLineWrap(true);
        citationText.setWrapStyleWord(true);
        citationText.setEditable(false);
        citationText.setFont(this.getFont());
        citationText.setOpaque(false);
        AminoAcidModelType type = (AminoAcidModelType) aaSubstCombo.getSelectedItem();
        citationText.setText(type.getCitation().toString());

        dolloCheck.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                if (dolloCheck.isSelected()) {
                    binarySubstCombo.setSelectedIndex(0);
                    binarySubstCombo.setEnabled(false);
                    useAmbiguitiesTreeLikelihoodCheck.setSelected(true);
                    useAmbiguitiesTreeLikelihoodCheck.setEnabled(false);
                    frequencyCombo.setEnabled(false);
                    frequencyCombo.setSelectedItem(FrequencyPolicyType.EMPIRICAL);
                    heteroCombo.setSelectedIndex(0);
                    heteroCombo.setEnabled(false);
                    model.setBinarySubstitutionModel(BinaryModelType.BIN_DOLLO);
                    model.setDolloModel(true);
                    DolloComponentOptions comp = (DolloComponentOptions)
                            model.getOptions().getComponentOptions(DolloComponentOptions.class);
                    comp.createParameters(model.getOptions());
                    comp.setActive(true);

                } else {
                    binarySubstCombo.setEnabled(true);
                    useAmbiguitiesTreeLikelihoodCheck.setEnabled(true);
                    frequencyCombo.setEnabled(true);
                    heteroCombo.setEnabled(true);
                    model.setBinarySubstitutionModel((BinaryModelType) binarySubstCombo.getSelectedItem());
                    model.setDolloModel(false);
                }
            }
        });

        PanelUtils.setupComponent(dolloCheck);
//        dolloCheck.addChangeListener(new ChangeListener() {
//            public void stateChanged(ChangeEvent e) {
//                model.setDolloModel(true);
//            }
//        });
        dolloCheck.setEnabled(true);
        dolloCheck
                .setToolTipText("<html>Activates a Stochastic Dollo model as described in<br>"
                        + "Alekseyenko, Lee & Suchard (2008) <i>Syst Biol</i> <b>57</b>: 772-784.</html>");

        PanelUtils.setupComponent(discreteTraitSiteModelCombo);
        discreteTraitSiteModelCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model
                        .setDiscreteSubstType((DiscreteSubstModelType) discreteTraitSiteModelCombo
                                .getSelectedItem());
                activateBSSVS.setEnabled(model.getDiscreteSubstType() != DiscreteSubstModelType.GLM_SUBST);
                setupGLMButton.setEnabled(model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST);
            }
        });

        PanelUtils.setupComponent(continuousTraitSiteModelCombo);
        continuousTraitSiteModelCombo
                .setToolTipText("<html>Select the model of continuous random walk, either homogenous<br>" +
                        "or relaxed random walk (RRW) with a choice of distributions.</html>");
        continuousTraitSiteModelCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model
                        .setContinuousSubstModelType((ContinuousSubstModelType) continuousTraitSiteModelCombo
                                .getSelectedItem());
            }
        });

        PanelUtils.setupComponent(latLongCheck);
        latLongCheck
                .setToolTipText("<html>Specify whether this is a geographical trait representing <br>"
                        + "latitude and longitude. Provides appropriate statistics to log file.</html>");

        latLongCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setIsLatitudeLongitude(latLongCheck.isSelected());
            }
        });
        latLongCheck.setEnabled(false);

        PanelUtils.setupComponent(addJitterCheck);
        addJitterCheck
                .setToolTipText("<html>Specify if the tip values should have some added random<br>" +
                        "noise. This can be useful if some tips have precisely<br>" +
                        "the same location.</html>");

        jitterWindowText.setValue(model.getJitterWindow());
        jitterWindowText.setColumns(10);

        addJitterCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                jitterWindowLabel.setEnabled(addJitterCheck.isSelected());
                jitterWindowText.setEnabled(addJitterCheck.isSelected());
                model.setJitterWindow(addJitterCheck.isSelected() ? jitterWindowText.getValue() : 0.0);
            }
        });
        jitterWindowText.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                model.setJitterWindow(addJitterCheck.isSelected() ? jitterWindowText.getValue() : 0.0);
            }
        });

        PanelUtils.setupComponent(useLambdaCheck);
        useLambdaCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                ContinuousComponentOptions component = (ContinuousComponentOptions) model.getOptions()
                        .getComponentOptions(ContinuousComponentOptions.class);
                component.setUseLambda(model, useLambdaCheck.isSelected());
            }
        });
        useLambdaCheck
                .setToolTipText("<html>Estimate degree of phylogenetic correlation in continuous traits using <br>"
                        + "a tree transform. Inspired by Pagel (1999), described in Lemey et al (2013) <i>in prep</i></html>");

        PanelUtils.setupComponent(activateBSSVS);
        activateBSSVS.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setActivateBSSVS(activateBSSVS.isSelected());
            }
        });
        activateBSSVS
                .setToolTipText("<html>Activates Bayesian stochastic search variable selection on the rates as described in<br>"
                        + "Lemey, Rambaut, Drummond & Suchard (2009) <i>PLoS Computational Biology</i> <b>5</b>: e1000520</html>");
        activateBSSVS.setEnabled(model.getDiscreteSubstType() != DiscreteSubstModelType.GLM_SUBST);

        setupGLMButton = new JButton("Setup GLM");
        setupGLMButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                setupGLM();
            }
        });
        PanelUtils.setupComponent(setupGLMButton);
        setupGLMButton
                .setToolTipText("<html>Set-up design of phylogenetic GLM.</html>");
        setupGLMButton.setEnabled(model.getDiscreteSubstType() == DiscreteSubstModelType.GLM_SUBST);

        // ============ micro-sat ================
        microsatName.setColumns(30);
        microsatName.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                model.getMicrosatellite().setName(microsatName.getText());
            }
        });
        microsatMax.setColumns(10);
        microsatMax.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                model.getMicrosatellite().setMax(
                        Integer.parseInt(microsatMax.getText()));
            }
        });
        microsatMin.setColumns(10);
        microsatMin.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                model.getMicrosatellite().setMin(
                        Integer.parseInt(microsatMin.getText()));
            }
        });

        PanelUtils.setupComponent(shareMicroSatCheck);
        shareMicroSatCheck.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                model.getOptions().shareMicroSat = shareMicroSatCheck
                        .isSelected();
                if (shareMicroSatCheck.isSelected()) {
                    model.getOptions().shareMicroSat();
                } else {
                    model.getOptions().unshareMicroSat();
                }
                setOptions();
            }
        });

        PanelUtils.setupComponent(rateProportionCombo);
        rateProportionCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model
                        .setRatePorportion((MicroSatModelType.RateProportionality) rateProportionCombo
                                .getSelectedItem());
            }
        });
        // rateProportionCombo.setToolTipText("<html>Select the type of microsatellite substitution model.</html>");
        PanelUtils.setupComponent(mutationBiasCombo);
        mutationBiasCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model
                        .setMutationBias((MicroSatModelType.MutationalBias) mutationBiasCombo
                                .getSelectedItem());
            }
        });
        PanelUtils.setupComponent(phaseCombo);
        phaseCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setPhase((MicroSatModelType.Phase) phaseCombo
                        .getSelectedItem());
            }
        });

        setupPanel();
        setOpaque(false);
    }

    /**
     * Sets the components up according to the partition model - but does not
     * layout the top level options panel.
     */
    public void setOptions() {

        if (SiteModelsPanel.DEBUG) {
            String modelName = (model == null) ? "null" : model.getName();
            Logger.getLogger("dr.app.beauti").info(
                    "ModelsPanel.setModelOptions(" + modelName + ")");
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
                binarySubstCombo
                        .setSelectedItem(model.getBinarySubstitutionModel());
                useAmbiguitiesTreeLikelihoodCheck.setSelected(model
                        .isUseAmbiguitiesTreeLikelihood());

                break;

            case DataType.GENERAL:
                discreteTraitSiteModelCombo.setSelectedItem(model
                        .getDiscreteSubstType());
                activateBSSVS.setSelected(model.getDiscreteSubstType() != DiscreteSubstModelType.GLM_SUBST ?
                        model.isActivateBSSVS() : false);
                break;

            case DataType.CONTINUOUS:
                continuousTraitSiteModelCombo.setSelectedItem(model
                        .getContinuousSubstModelType());

                ContinuousComponentOptions component = (ContinuousComponentOptions) model.getOptions()
                        .getComponentOptions(ContinuousComponentOptions.class);

                latLongCheck.setSelected(model.isLatitudeLongitude());
                latLongCheck.setEnabled(model.getContinuousTraitCount() == 2);
                useLambdaCheck.setSelected(component.useLambda(model));
                break;
            case DataType.MICRO_SAT:
                microsatName.setText(model.getMicrosatellite().getName());
                microsatMax.setText(Integer.toString(model.getMicrosatellite()
                        .getMax()));
                microsatMin.setText(Integer.toString(model.getMicrosatellite()
                        .getMin()));
                shareMicroSatCheck.setSelected(model.getOptions().shareMicroSat);
                rateProportionCombo.setSelectedItem(model.getRatePorportion());
                mutationBiasCombo.setSelectedItem(model.getMutationBias());
                phaseCombo.setSelectedItem(model.getPhase());
                shareMicroSatCheck.setEnabled(model.getOptions()
                        .getPartitionSubstitutionModels(Microsatellite.INSTANCE)
                        .size() > 1);
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
     * Configure this panel for the Yang 96 codon
     * position model
     */
    private void setYang96Model() {
        nucSubstCombo.setSelectedIndex(2);
        heteroCombo.setSelectedIndex(1);
        codingCombo.setSelectedIndex(2);
        substUnlinkCheck.setSelected(true);
        heteroUnlinkCheck.setSelected(true);
        freqsUnlinkCheck.setSelected(true);
    }

    /**
     * Configure this panel for the Shapiro, Rambaut and Drummond 2006 codon
     * position model
     */
    private void setSRD06Model() {
        nucSubstCombo.setSelectedIndex(1);
        heteroCombo.setSelectedIndex(1);
        codingCombo.setSelectedIndex(1);
        substUnlinkCheck.setSelected(true);
        heteroUnlinkCheck.setSelected(true);
        freqsUnlinkCheck.setSelected(false);
    }

    private void setupGLM() {
        if (glmSettingsDialog == null) {
            glmSettingsDialog = new GLMSettingsDialog(frame);
        }

        glmSettingsDialog.setTrait(model.getTraitData());

        int result = glmSettingsDialog.showDialog();

        if (result == JOptionPane.OK_OPTION) {
            // Only do this if OK button is pressed (not cancel):

            frame.setAllOptions();
        }
    }

    /**
     * Lays out the appropriate components in the panel for this partition
     * model.
     */
    private void setupPanel() {

        switch (model.getDataType().getType()) {
            case DataType.NUCLEOTIDES:
                addComponentWithLabel("Substitution Model:", nucSubstCombo);
                addComponentWithLabel("Base frequencies:", frequencyCombo);
                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                heteroCombo.setSelectedIndex(0);
                gammaCatLabel = addComponentWithLabel(
                        "Number of Gamma Categories:", gammaCatCombo);
                gammaCatCombo.setEnabled(false);

                addSeparator();

                addComponentWithLabel("Partition into codon positions:",
                        codingCombo);

                JPanel panel2 = new JPanel();
                panel2.setOpaque(false);
                panel2.setLayout(new BoxLayout(panel2, BoxLayout.PAGE_AXIS));
                panel2.setBorder(new TitledBorder(null, "Link/Unlink parameters:", TitledBorder.DEFAULT_JUSTIFICATION, TitledBorder.ABOVE_TOP));
                panel2.add(substUnlinkCheck);
                panel2.add(heteroUnlinkCheck);
                panel2.add(freqsUnlinkCheck);

                addComponent(panel2);

                addComponent(setYang96Button);
                addComponent(setSRD06Button);

                break;

            case DataType.AMINO_ACIDS:
                addComponentWithLabel("Substitution Model:", aaSubstCombo);
                addComponentWithLabel("Citation:", citationText);

                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                heteroCombo.setSelectedIndex(0);
                gammaCatLabel = addComponentWithLabel(
                        "Number of Gamma Categories:", gammaCatCombo);
                gammaCatCombo.setEnabled(false);

                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                addComponentWithLabel("Substitution Model:", binarySubstCombo);
                addComponentWithLabel("Base frequencies:", frequencyCombo);
                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                heteroCombo.setSelectedIndex(0);
                gammaCatLabel = addComponentWithLabel(
                        "Number of Gamma Categories:", gammaCatCombo);
                gammaCatCombo.setEnabled(false);

                addSeparator();

                addComponentWithLabel("", useAmbiguitiesTreeLikelihoodCheck);

                // Easy XML specification is currently only available for binary models
                if (ENABLE_STOCHASTIC_DOLLO) {
                    addSeparator();
                    addComponent(dolloCheck);
                }

                break;

            case DataType.GENERAL:
                addComponentWithLabel("Discrete Trait Substitution Model:",
                        discreteTraitSiteModelCombo);
                addComponent(activateBSSVS);
                addComponent(setupGLMButton);
                break;

            case DataType.CONTINUOUS:
                addComponentWithLabel("Continuous Trait Model:",
                        continuousTraitSiteModelCombo);
                addComponent(latLongCheck);
                addSeparator();
                addComponent(addJitterCheck);
                OptionsPanel panel = new OptionsPanel();
                panel.addComponents(jitterWindowLabel, jitterWindowText);
                addComponent(panel);
                jitterWindowLabel.setEnabled(false);
                jitterWindowText.setEnabled(false);

                addSeparator();
                addComponent(useLambdaCheck);
                break;

            case DataType.MICRO_SAT:
                addComponentWithLabel("Microsatellite Name:", microsatName);
                addComponentWithLabel("Max of Length:", microsatMax);
                addComponentWithLabel("Min of Length:", microsatMin);
                addComponent(shareMicroSatCheck);

                addSeparator();

                addComponentWithLabel("Rate Proportionality:", rateProportionCombo);
                addComponentWithLabel("Mutational Bias:", mutationBiasCombo);
                addComponentWithLabel("Phase:", phaseCombo);
                break;

            default:
                throw new IllegalArgumentException("Unknown data type");

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
                model.setUnlinkedSubstitutionModel(substUnlinkCheck
                        .isSelected());
            }
        });
        substUnlinkCheck.setEnabled(false);
        substUnlinkCheck.setToolTipText(""
                + "<html>Gives each codon position partition different<br>"
                + "substitution model parameters.</html>");

        PanelUtils.setupComponent(heteroUnlinkCheck);
        heteroUnlinkCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setUnlinkedHeterogeneityModel(heteroUnlinkCheck
                        .isSelected());
            }
        });
        heteroUnlinkCheck.setEnabled(heteroCombo.getSelectedIndex() != 0);
        heteroUnlinkCheck
                .setToolTipText("<html>Gives each codon position partition different<br>rate heterogeneity model parameters.</html>");

        PanelUtils.setupComponent(freqsUnlinkCheck);
        freqsUnlinkCheck.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                model.setUnlinkedFrequencyModel(freqsUnlinkCheck.isSelected());
            }
        });
        freqsUnlinkCheck.setEnabled(false);
        freqsUnlinkCheck
                .setToolTipText("<html>Gives each codon position partition different<br>nucleotide frequency parameters.</html>");

        PanelUtils.setupComponent(codingCombo);
        codingCombo
                .setToolTipText("<html>Select how to partition the codon positions.</html>");
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

                if (codingCombo.getSelectedIndex() != 0) {
                    // codon position partitioning
                    substUnlinkCheck.setEnabled(true);
                    heteroUnlinkCheck
                            .setEnabled(heteroCombo.getSelectedIndex() != 3);
                    freqsUnlinkCheck.setEnabled(true);
                    substUnlinkCheck.setSelected(true);
                    heteroUnlinkCheck.setSelected(heteroCombo
                            .getSelectedIndex() != 0);
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
        });
    }
}
