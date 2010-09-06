package dr.app.phylogeography.spread.layerspanel;

import dr.app.phylogeography.builder.Builder;
import dr.app.phylogeography.builder.BuilderFactory;
import dr.app.phylogeography.spread.SpreadDocument;
import dr.app.phylogeography.spread.InputFile;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class CreateBuilderDialog {

    private JFrame frame;

    private JComboBox dataFileCombo;
    private JComboBox builderCombo;

    private OptionsPanel optionPanel;

    private static final String SELECT_LAYER_TYPE = "Select layer type...";
    private static final String SELECT_INPUT_FILE = "Select input file...";

    public CreateBuilderDialog(JFrame frame) {
        this.frame = frame;

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel() {
        optionPanel.removeAll();

        if (builderCombo != null) {
            optionPanel.addComponentWithLabel("Layer type:", builderCombo);
        }

        if (dataFileCombo != null) {
            optionPanel.addComponentWithLabel("Input file:", dataFileCombo);
        }
    }

    public int showDialog(BuilderFactory[] builderFactories, SpreadDocument document) {

        builderCombo = new JComboBox();
        builderCombo.addItem(SELECT_LAYER_TYPE);
        for (Object factory : builderFactories) {
            builderCombo.addItem(factory);
        }

        dataFileCombo = new JComboBox();
        dataFileCombo.addItem(SELECT_INPUT_FILE);
        for (Object dataFile : document.getInputFiles()) {
            dataFileCombo.addItem(dataFile);
        }

        setupPanel();

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Layer");
        dialog.pack();

//        builderCombo.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent e) {
//                if (builderCombo.getSelectedIndex() != 0) {
//                    builder = ((BuilderFactory)builderCombo.getSelectedItem()).createBuilder();
//                    setupPanel(builder);
//                } else {
//                    setupPanel(null);
//                }
//                dialog.pack();
//            }
//        });

        int result;
        boolean validated = false;

        do {
            result = JOptionPane.CANCEL_OPTION;

            dialog.setVisible(true);

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }

            if (result != JOptionPane.CANCEL_OPTION) {
                if ((dataFileCombo == null || !dataFileCombo.getSelectedItem().toString().equals(SELECT_LAYER_TYPE)) &&
                        (builderCombo == null || !builderCombo.getSelectedItem().toString().equals(SELECT_INPUT_FILE))) {
                    validated = true;
                } else {
                    Toolkit.getDefaultToolkit().beep();
                }
            }

        } while (!validated && result != JOptionPane.CANCEL_OPTION);

        return result;
    }

    public Builder getBuilder() {
        Builder builder = ((BuilderFactory)builderCombo.getSelectedItem()).createBuilder();
        builder.setName("Layer " + BuilderFactory.nextCount());
        builder.setInputFile((InputFile)dataFileCombo.getSelectedItem());
        builder.setFromEditPanel();
        return builder;
    }

}