package dr.app.beauti.modelsPanel;

import dr.app.beauti.PanelUtils;
import dr.app.beauti.BeautiApp;
import dr.app.beauti.options.AminoAcidModelType;
import dr.app.beauti.options.FrequencyPolicy;
import dr.app.beauti.options.NucModelType;
import dr.app.beauti.options.PartitionModel;
import dr.evolution.datatype.DataType;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 */
public class PartitionModelPanel extends OptionsPanel {

    PartitionModel model;

    public PartitionModelPanel(PartitionModel partitionModel) {

        super(12, 18);

        this.model = partitionModel;

        initCodonPartitionComponents();

        PanelUtils.setupComponent(nucSubstCombo);
        nucSubstCombo.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                if (nucSubstCombo.getSelectedIndex() == 0) {
                    model.setNucSubstitutionModel(NucModelType.HKY);
                } else {
                    model.setNucSubstitutionModel(NucModelType.GTR);
                }
            }
        });
        nucSubstCombo.setToolTipText("<html>Select the type of nucleotide substitution model.</html>");

        PanelUtils.setupComponent(aaSubstCombo);
        aaSubstCombo.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                model.setAaSubstitutionModel((AminoAcidModelType) aaSubstCombo.getSelectedItem());
            }
        });
        aaSubstCombo.setToolTipText("<html>Select the type of amino acid substitution model.</html>");

        PanelUtils.setupComponent(binarySubstCombo);
        binarySubstCombo.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                model.setBinarySubstitutionModel(binarySubstCombo.getSelectedIndex());
            }
        });
        binarySubstCombo.setToolTipText("<html>Select the type of binary substitution model.</html>");

        PanelUtils.setupComponent(frequencyCombo);
        frequencyCombo.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                model.setFrequencyPolicy((FrequencyPolicy) frequencyCombo.getSelectedItem());
            }
        });
        frequencyCombo.setToolTipText("<html>Select the policy for determining the base frequencies.</html>");

        PanelUtils.setupComponent(heteroCombo);
        heteroCombo.setToolTipText("<html>Select the type of site-specific rate<br>heterogeneity model.</html>");
        heteroCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {

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
                    }
                }
        );

        PanelUtils.setupComponent(gammaCatCombo);
        gammaCatCombo.setToolTipText("<html>Select the number of categories to use for<br>the discrete gamma rate heterogeneity model.</html>");
        gammaCatCombo.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {

                model.setGammaCategories(gammaCatCombo.getSelectedIndex() + 4);
            }
        });

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
    private void setupComponents() {

        if (ModelsPanel.DEBUG) {
            String modelName = (model == null) ? "null" : model.getName();
            Logger.getLogger("dr.app.beauti").info("ModelsPanel.setModelOptions(" + modelName + ")");
        }

        if (model == null) {
            return;
        }

        int dataType = model.dataType.getType();
        switch (dataType) {
            case DataType.NUCLEOTIDES:
                if (model.getNucSubstitutionModel() == NucModelType.GTR) {
                    nucSubstCombo.setSelectedIndex(1);
                } else {
                    nucSubstCombo.setSelectedIndex(0);
                }

                frequencyCombo.setSelectedItem(model.getFrequencyPolicy());

                break;

            case DataType.AMINO_ACIDS:
                aaSubstCombo.setSelectedItem(model.getAaSubstitutionModel());
                break;

            case DataType.TWO_STATES:
            case DataType.COVARION:
                binarySubstCombo.setSelectedIndex(model.getBinarySubstitutionModel());
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

        switch (model.dataType.getType()) {
            case DataType.NUCLEOTIDES:
                addComponentWithLabel("Substitution Model:", nucSubstCombo);
                addComponentWithLabel("Base frequencies:", frequencyCombo);
                addComponentWithLabel("Site Heterogeneity Model:", heteroCombo);
                gammaCatLabel = addComponentWithLabel("Number of Gamma Categories:", gammaCatCombo);

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

        if (BeautiApp.advanced) {
            addSeparator();
            addComponent(dolloCheck);
        }

        setupComponents();
    }

    /**
     * Initializes and binds the components related to modeling codon positions.
     */
    private void initCodonPartitionComponents() {

        PanelUtils.setupComponent(substUnlinkCheck);

        substUnlinkCheck.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                model.setUnlinkedSubstitutionModel(true);
            }
        });
        substUnlinkCheck.setEnabled(false);
        substUnlinkCheck.setToolTipText("" +
                "<html>Gives each codon position partition different<br>" +
                "substitution model parameters.</html>");

        PanelUtils.setupComponent(heteroUnlinkCheck);
        heteroUnlinkCheck.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                model.setUnlinkedHeterogeneityModel(true);
            }
        });
        heteroUnlinkCheck.setEnabled(false);
        heteroUnlinkCheck.setToolTipText("<html>Gives each codon position partition different<br>rate heterogeneity model parameters.</html>");

        PanelUtils.setupComponent(freqsUnlinkCheck);
        freqsUnlinkCheck.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                model.setUnlinkedFrequencyModel(true);
            }
        });
        freqsUnlinkCheck.setEnabled(false);
        freqsUnlinkCheck.setToolTipText("<html>Gives each codon position partition different<br>nucleotide frequency parameters.</html>");

        PanelUtils.setupComponent(codingCombo);
        codingCombo.setToolTipText("<html>Select how to partition the codon positions.</html>");
        codingCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {

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
    }

    // Actions

    private Action setSRD06Action = new AbstractAction("Use SRD06 Model") {
        public void actionPerformed(ActionEvent actionEvent) {
            setSRD06Model();
        }
    };

    // Components

    private JComboBox nucSubstCombo = new JComboBox(new String[]{"HKY", "GTR"});
    private JComboBox aaSubstCombo = new JComboBox(AminoAcidModelType.values());
    private JComboBox binarySubstCombo = new JComboBox(new String[]{"Simple", "Covarion"});

    private JComboBox frequencyCombo = new JComboBox(FrequencyPolicy.values());

    private JComboBox heteroCombo = new JComboBox(
            new String[]{"None", "Gamma", "Invariant Sites", "Gamma + Invariant Sites"});

    private JComboBox gammaCatCombo = new JComboBox(new String[]{"4", "5", "6", "7", "8", "9", "10"});
    private JLabel gammaCatLabel;

    private JComboBox codingCombo = new JComboBox(new String[]{
            "Off",
            "2 partitions: codon positions (1 + 2), 3",
            "3 partitions: codon positions 1, 2, 3"});

    private JCheckBox substUnlinkCheck = new JCheckBox("Unlink substitution model across codon positions");
    private JCheckBox heteroUnlinkCheck =
            new JCheckBox("Unlink rate heterogeneity model across codon positions");
    private JCheckBox freqsUnlinkCheck = new JCheckBox("Unlink base frequencies across codon positions");

    private JButton setSRD06Button;

    private JCheckBox dolloCheck = new JCheckBox("Use Stochastic Dollo Model");
    // private JComboBox dolloCombo = new JComboBox(new String[]{"Analytical", "Sample"});

}
