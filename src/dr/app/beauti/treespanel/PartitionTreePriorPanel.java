/*
 * PartitionTreePriorPanel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.beauti.treespanel;

import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.types.PopulationSizeModelType;
import dr.app.beauti.types.TreePriorParameterizationType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.components.WholeNumberField;
import dr.app.util.OSType;
import dr.evolution.util.Taxa;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.BirthDeathSerialSamplingModelParser;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionTreePriorPanel extends OptionsPanel {

    private static final long serialVersionUID = 5016996360264782252L;

    private final JComboBox treePriorCombo = new JComboBox();

    private final JComboBox parameterizationCombo = new JComboBox(EnumSet.range(TreePriorParameterizationType.GROWTH_RATE,
            TreePriorParameterizationType.DOUBLING_TIME).toArray());
    //    private JComboBox parameterizationCombo1 = new JComboBox(EnumSet.of(TreePriorParameterizationType.DOUBLING_TIME).toArray());
    private final JComboBox gmrfBayesianSkyrideCombo = new JComboBox(EnumSet.range(TreePriorParameterizationType.UNIFORM_SKYRIDE,
            TreePriorParameterizationType.TIME_AWARE_SKYRIDE).toArray());

    private final WholeNumberField skyGridPointsField = new WholeNumberField(2, Integer.MAX_VALUE);
    private final RealNumberField skyGridInterval = new RealNumberField(0.0, Double.MAX_VALUE);

//    private JComboBox skyGridCombo = new JComboBox(EnumSet.range())

    private final JComboBox populationSizeCombo = new JComboBox(PopulationSizeModelType.values());

    private final JCheckBox subtreePriorCheckBox = new JCheckBox("Use seperate prior for subtree");
    private final JComboBox subtreeTaxonSetCombo = new JComboBox();
    private final JLabel subtreeTaxonSetComboLabel = new JLabel("Taxon set to define subtree:");
    private final JComboBox subtreePriorCombo = new JComboBox();
    private final JLabel subtreePriorComboLabel = new JLabel("Subtree prior model:");

//    private JComboBox calibrationCorrectionCombo = new JComboBox(new CalibrationPoints.CorrectionType[]
//            {CalibrationPoints.CorrectionType.EXACT, CalibrationPoints.CorrectionType.NONE});

//    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

//	private BeautiFrame frame = null;
//	private BeautiOptions options = null;

    PartitionTreePrior partitionTreePrior;
    private final TreesPanel treesPanel;

    private boolean settingOptions = false;


    public PartitionTreePriorPanel(PartitionTreePrior partitionTreePrior, final TreesPanel parent) {
        super(12, (OSType.isMac() ? 6 : 24));

        this.partitionTreePrior = partitionTreePrior;
        this.treesPanel = parent;

        setTreePriorChoices(false);
        PanelUtils.setupComponent(treePriorCombo);
        treePriorCombo.setMaximumRowCount(10); // to show Calibrated Yule
        treePriorCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {

                if (treePriorCombo.getSelectedItem() != null) {
                    setupPanel();

                    PartitionTreePriorPanel.this.partitionTreePrior.setNodeHeightPrior((TreePriorType) treePriorCombo.getSelectedItem());
                    parent.fireTreePriorsChanged();
                }
            }
        });

        PanelUtils.setupComponent(parameterizationCombo);
        parameterizationCombo.setToolTipText("<html>" +
                "Select the parameterization of growth rate to use for<br>" +
                "the parametric coalescent model. This does not affect<br>" +
                "nature of the model but may assist with mixing.<html>");
        parameterizationCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                PartitionTreePriorPanel.this.partitionTreePrior.setParameterization((TreePriorParameterizationType) parameterizationCombo.getSelectedItem());
                parent.fireTreePriorsChanged();
            }
        });

        PanelUtils.setupComponent(skyGridInterval);
        skyGridInterval.setToolTipText("<html>" +
                "This sets the time interval over which the change-points<br>" +
                "of the Skygrid are placed. It is measured relative to the<br>" +
                "most recent tip and should be commensurate with the predicted<br>" +
                "age of the tree.<html>");
        skyGridInterval.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
            }
        });

        PanelUtils.setupComponent(gmrfBayesianSkyrideCombo);
        gmrfBayesianSkyrideCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent ev) {
                PartitionTreePriorPanel.this.partitionTreePrior.setSkyrideSmoothing((TreePriorParameterizationType) gmrfBayesianSkyrideCombo.getSelectedItem());
                parent.fireTreePriorsChanged();
            }
        });

//        PanelUtils.setupComponent(skyGridPointsCombo);
//        skyGridPointsCombo.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent itemEvent) {
//                partitionTreePrior.setSkyGridCount((Integer) skyGridPointsCombo.getSelectedItem());
//                parent.fireTreePriorsChanged();
//            }
//        });
        PanelUtils.setupComponent(skyGridPointsField);
        skyGridPointsField.setToolTipText("<html>" +
                "This sets number of population size parameters for<br>" +
                "the Skygrid. This determines the number of transition-points<br>" +
                "at which the population size changes. The last one spans to<br>" +
                "infinite time.<html>");
        skyGridPointsField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
            }
        });

        PanelUtils.setupComponent(populationSizeCombo);
        populationSizeCombo.addItemListener(new ItemListener() {
                                                public void itemStateChanged(ItemEvent ev) {
                                                    PartitionTreePriorPanel.this.partitionTreePrior.setPopulationSizeModel((PopulationSizeModelType) populationSizeCombo.getSelectedItem());
                                                    parent.fireTreePriorsChanged();
                                                }
                                            }
        );

        PanelUtils.setupComponent(subtreePriorCheckBox);
        subtreePriorCheckBox.setToolTipText("<html>" +
                "Select this to specify a different tree prior model<br>" +
                "for a specific subset of the tree as defined by a<br>" +
                "taxon set. This should be specified as monophyletic.<br>" +
                "<br>" +
                "The rest of the tree will be have the prior defined above.<html>");
        subtreePriorCheckBox.addItemListener(new ItemListener() {
                                                 public void itemStateChanged(ItemEvent ev) {
                                                     boolean selected = subtreePriorCheckBox.isSelected();
                                                     subtreeTaxonSetCombo.setEnabled(selected);
                                                     subtreeTaxonSetComboLabel.setEnabled(selected);
                                                     subtreePriorCombo.setEnabled(selected);
                                                     subtreePriorComboLabel.setEnabled(selected);

                                                     if (selected &&
                                                             subtreeTaxonSetCombo.getSelectedItem() != null &&
                                                             !subtreeTaxonSetCombo.getSelectedItem().equals("No monophyletic taxon sets defined") ) {
                                                         PartitionTreePriorPanel.this.partitionTreePrior.setSubtreeTaxonSet((Taxa)subtreeTaxonSetCombo.getSelectedItem());
                                                         PartitionTreePriorPanel.this.partitionTreePrior.setSubtreePrior((TreePriorType)subtreePriorCombo.getSelectedItem());
                                                     } else {
                                                         PartitionTreePriorPanel.this.partitionTreePrior.setSubtreeTaxonSet(null);
                                                     }
                                                     parent.fireTreePriorsChanged();
                                                 }
                                             }
        );

        subtreeTaxonSetCombo.addItem("No monophyletic taxon sets defined");
        PanelUtils.setupComponent(subtreeTaxonSetCombo);
        ItemListener listener = ev -> {
            if (subtreeTaxonSetCombo.getSelectedItem() != null && !subtreeTaxonSetCombo.getSelectedItem().equals("No monophyletic taxon sets defined")) {
                PartitionTreePriorPanel.this.partitionTreePrior.setSubtreeTaxonSet((Taxa) subtreeTaxonSetCombo.getSelectedItem());
                PartitionTreePriorPanel.this.partitionTreePrior.setSubtreePrior((TreePriorType) subtreePriorCombo.getSelectedItem());
            } else {
                PartitionTreePriorPanel.this.partitionTreePrior.setSubtreeTaxonSet(null);
            }
            parent.fireTreePriorsChanged();
        };

        subtreeTaxonSetCombo.addItemListener(listener);
        subtreeTaxonSetCombo.setEnabled(false);
        subtreeTaxonSetComboLabel.setEnabled(false);

        subtreePriorCombo.addItem(TreePriorType.CONSTANT);
        subtreePriorCombo.addItem(TreePriorType.EXPONENTIAL);
        PanelUtils.setupComponent(subtreePriorCombo);
        subtreePriorCombo.addItemListener(listener);
        subtreePriorCombo.setEnabled(false);
        subtreePriorComboLabel.setEnabled(false);

        setOptions();
    }

    private void setupPanel() {

        removeAll();

        JTextArea citationText = new JTextArea(1, 40);
        citationText.setLineWrap(true);
        citationText.setWrapStyleWord(true);
        citationText.setEditable(false);
        citationText.setFont(this.getFont());
        citationText.setOpaque(false);
//        citationText.setBackground(this.getBackground());
//        JScrollPane scrollPane = new JScrollPane(citation, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
//                    JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
//        scrollPane.setOpaque(true);
        String calYule = "Heled J, Drummond AJ (2011), Syst Biol, doi: 10.1093/sysbio/syr087 [Calibrated Yule]";
        String citation;

        String citationCoalescent = "Kingman JFC (1982) Stoch Proc Appl 13, 235-248 [Constant Coalescent].";

        addComponentWithLabel("Tree Prior:", treePriorCombo);

        if (!treesPanel.linkTreePriorCheck.isEnabled()) {
            treesPanel.updateLinkTreePriorEnablility();
        }

        TreePriorType type = (TreePriorType) treePriorCombo.getSelectedItem();
        switch (type) {
            case CONSTANT:
                citation = citationCoalescent;
                break;

            case EXPONENTIAL:
            case LOGISTIC:
            case EXPANSION:
                addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
                partitionTreePrior.setParameterization((TreePriorParameterizationType) parameterizationCombo.getSelectedItem());

                citation = //citationCoalescent +  "\n" +
                        "Griffiths RC, Tavare S (1994) Phil Trans R Soc Lond B Biol Sci 344, 403-410 [Parametric Coalescent].";
//                        + "\nDrummond AJ, Rambaut A, Shapiro B, Pybus OG (2005) Mol Biol Evol 22, 1185-1192.";
                break;

//            case SKYLINE:
//                groupCountField.setColumns(6);
//                addComponentWithLabel("Number of groups:", groupCountField);
//                addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
//
//                citation = //citationCoalescent + "\n" +
//                        "Drummond AJ, Rambaut A, Shapiro B, Pybus OG (2005) Mol Biol Evol 22, 1185-1192 [Skyline Coalescent].";
//                break;

            case GMRF_SKYRIDE:
                addComponentWithLabel("Smoothing:", gmrfBayesianSkyrideCombo);
                treesPanel.linkTreePriorCheck.setSelected(true);
                treesPanel.linkTreePriorCheck.setEnabled(false);
                //For GMRF, one tree prior has to be associated to one tree model. The validation is in BeastGenerator.checkOptions()
                addLabel("<html>For the Skyride, tree model/tree prior combination not implemented by BEAST. "
                        + "The Skyride is only available for a single tree<br>model partition in this release. "
                        + "Please try the Skygrid or link all tree models." + "</html>");

                citation = //citationCoalescent + "\n" +
                        "Minin VN, Bloomquist EW, Suchard MA (2008) Mol Biol Evol 25, 1459-1471 [Skyride Coalescent].";
                break;

            case SKYGRID:
            case SKYGRID_HMC:
                skyGridPointsField.setColumns(6);
                addComponentWithLabel("Number of parameters:", skyGridPointsField);
                skyGridInterval.setColumns(6);
                addComponentWithLabel("Time at last transition point:", skyGridInterval);
                treesPanel.linkTreePriorCheck.setSelected(true);
                treesPanel.linkTreePriorCheck.setEnabled(false);
                treesPanel.updateShareSameTreePriorChanged();

                citation = //citationCoalescent + "\n" +
                        "Gill MS, Lemey P, Faria NR, Rambaut A, Shapiro B, Suchard MA (2013) Mol Biol Evol 30, 713-724 [SkyGrid Coalescent].";

                if (type == TreePriorType.SKYGRID_HMC) {
                    citation += "\nBaele G, Gill MS, Lemey P & Suchard MA (2020) Wellcome Open Res 5:53 [Hamiltonian SkyGrid]";
                }
                break;


            case YULE:
                citation = "Gernhard T (2008) J Theor Biol 253, 769-778 [Yule Process]." +
                        "\nYule GU (1925) Phil Trans R Soc Lond B Biol Sci 213, 21-87 [Yule Process].";
                break;

            case YULE_CALIBRATION:
//                    addComponentWithLabel("Calibration Correction Type:", calibrationCorrectionCombo);
                citation = calYule;
                break;

            case BIRTH_DEATH:
                citation = BirthDeathModelParser.getCitation();
                break;

            case BIRTH_DEATH_INCOMPLETE_SAMPLING:
                citation = BirthDeathModelParser.getCitationRHO();
                break;

            case BIRTH_DEATH_SERIAL_SAMPLING:
                citation = BirthDeathSerialSamplingModelParser.getCitationPsiOrg();
                break;

//            case BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER:
//                citation = BirthDeathSerialSamplingModelParser.getCitationRT();
//                break;

            default:
                throw new RuntimeException("No such tree prior has been specified so cannot refer to it");
        }

        if (treesPanel.options.maximumTipHeight > 0)
            citation = citation
//                    + "\n" +
//                    "Rodrigo AG, Felsenstein J (1999) in Molecular Evolution of HIV (Crandall K), pp. 233-272 [Serially Sampled Data]."
                    + "\n" +
                    "Drummond AJ, Nicholls GK, Rodrigo AG, Solomon W (2002) Genetics 161, 1307-1320 [Serially Sampled Data].";

//        addComponentWithLabel("Citation:", citationText);
//        citationText.setText(citation);

        for (PartitionTreeModel model : treesPanel.treeModelPanels.keySet()) {
            if (model != null) {
                treesPanel.treeModelPanels.get(model).setOptions();
                treesPanel.treeModelPanels.get(model).setupPanel();
            }
        }

        addSeparator();
        addSpanningComponent(subtreePriorCheckBox);
        addComponents(subtreeTaxonSetComboLabel, subtreeTaxonSetCombo);
        addComponents(subtreePriorComboLabel, subtreePriorCombo);

        validate();
        repaint();
    }

    public void setOptions() {

        if (settingOptions) {
            // make sure we don't go all recursive
            return;
        }

        if (partitionTreePrior == null) {
            return;
        }

        settingOptions = true;

        treePriorCombo.setSelectedItem(partitionTreePrior.getNodeHeightPrior());


        parameterizationCombo.setSelectedItem(partitionTreePrior.getParameterization());

        gmrfBayesianSkyrideCombo.setSelectedItem(partitionTreePrior.getSkyrideSmoothing());

//        skyGridPointsCombo.setSelectedItem(partitionTreePrior.getSkyGridCount());
        skyGridPointsField.setValue(partitionTreePrior.getSkyGridCount());

        double initialCutOff = partitionTreePrior.getOptions().getPartitionTreeModels().get(0).getInitialRootHeight();
        final double arbitraryScalar = 1.0;
        initialCutOff = roundToSignificantFigures(arbitraryScalar * initialCutOff, 2);

        skyGridInterval.setAllowEmpty(true);

        if (Double.isNaN(partitionTreePrior.getSkyGridInterval())) {
            if (initialCutOff > 0.0) {
                skyGridInterval.setValue(initialCutOff);
            } else {
                skyGridInterval.setText((Double)null);
            }
        } else {
            skyGridInterval.setValue(partitionTreePrior.getSkyGridInterval());
        }

        populationSizeCombo.setSelectedItem(partitionTreePrior.getPopulationSizeModel());

//        calibrationCorrectionCombo.setSelectedItem(partitionTreePrior.getCalibCorrectionType());

        subtreeTaxonSetCombo.removeAllItems();
        for (Taxa taxonSet : partitionTreePrior.getOptions().taxonSets) {
            if (taxonSet.getTaxonCount() > 1 && partitionTreePrior.getOptions().taxonSetsMono.get(taxonSet)) {
                subtreeTaxonSetCombo.addItem(taxonSet);
            }
        }
        if (subtreeTaxonSetCombo.getItemCount() == 0) {
            subtreeTaxonSetCombo.addItem("No monophyletic taxon sets defined");
        }

        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    private static double roundToSignificantFigures(double num, int n) {
        if (num == 0) {
            return 0;
        }

        final double d = Math.ceil(Math.log10(num < 0 ? -num : num));
        final int power = n - (int) d;

        final double magnitude = Math.pow(10, power);
        final long shifted = Math.round(num * magnitude);
        return shifted / magnitude;
    }

    public void getOptions() {
        if (settingOptions) return;

//        partitionTreePrior.setNodeHeightPrior((TreePriorType) treePriorCombo.getSelectedItem());

//        if (partitionTreePrior.getNodeHeightPrior() == TreePriorType.SKYLINE) {
//            Integer groupCount = groupCountField.getValue();
//            if (groupCount != null) {
//                partitionTreePrior.setSkylineGroupCount(groupCount);
//            } else {
//                partitionTreePrior.setSkylineGroupCount(5);
//            }
//        } else
        if (partitionTreePrior.getNodeHeightPrior() == TreePriorType.SKYGRID ||
                partitionTreePrior.getNodeHeightPrior() == TreePriorType.SKYGRID_HMC) {
            Double interval = skyGridInterval.getValue();
            if (interval != null) {
                partitionTreePrior.setSkyGridInterval(interval);
            } else {
                partitionTreePrior.setSkyGridInterval(Double.NaN);
            }
            partitionTreePrior.setSkyGridCount(skyGridPointsField.getValue());
        } else if (partitionTreePrior.getNodeHeightPrior() == TreePriorType.BIRTH_DEATH) {
//            Double samplingProportion = samplingProportionField.getValue();
//            if (samplingProportion != null) {
//                partitionTreePrior.birthDeathSamplingProportion = samplingProportion;
//            } else {
//                partitionTreePrior.birthDeathSamplingProportion = 1.0;
//            }
        }

    }

    public void setTreePriorChoices(boolean isTipCalibrated) {
        TreePriorType type = (TreePriorType) treePriorCombo.getSelectedItem();
        treePriorCombo.removeAllItems();


        for (TreePriorType treePriorType : EnumSet.range(TreePriorType.CONSTANT, TreePriorType.BIRTH_DEATH_SERIAL_SAMPLING)) {
            treePriorCombo.addItem(treePriorType);
            if (treePriorType == TreePriorType.EXPANSION ||
                    treePriorType == TreePriorType.GMRF_SKYRIDE) {
                treePriorCombo.addItem(new JSeparator(JSeparator.HORIZONTAL));
            }
        }

        if (isTipCalibrated) {
            // remove models that require contemporaneous tips...
            treePriorCombo.removeItem(TreePriorType.YULE);
            treePriorCombo.removeItem(TreePriorType.YULE_CALIBRATION);
            treePriorCombo.removeItem(TreePriorType.BIRTH_DEATH);
            treePriorCombo.removeItem(TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING);
        }

        // this makes sure treePriorCombo selects correct prior
        treePriorCombo.setSelectedItem(type);
        if (treePriorCombo.getSelectedItem() == null) {
            treePriorCombo.setSelectedIndex(0);
        }
    }
}

