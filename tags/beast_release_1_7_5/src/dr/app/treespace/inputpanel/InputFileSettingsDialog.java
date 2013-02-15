package dr.app.treespace.inputpanel;

import dr.app.gui.components.RealNumberField;
import jam.panels.OptionsPanel;
import dr.app.gui.components.WholeNumberField;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

import dr.app.treespace.InputFile;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class InputFileSettingsDialog {

    private JFrame frame;

    private WholeNumberField burninField;
    private RealNumberField mrsdField;

    private OptionsPanel optionPanel;

    private InputFile inputFile = null;

    public InputFileSettingsDialog(JFrame frame) {
        this.frame = frame;

        burninField = new WholeNumberField(0, Integer.MAX_VALUE);
        burninField.setColumns(8);

        mrsdField = new RealNumberField(0, Double.POSITIVE_INFINITY);
        mrsdField.setColumns(20);

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel(InputFile inputFile) {
        optionPanel.removeAll();

        if (inputFile.getTree() != null) {
            if (inputFile.getType() == InputFile.Type.POSTERIOR_TREES) {
                optionPanel.addComponentWithLabel("Burnin (number of trees):", burninField);
                burninField.setValue(inputFile.getBurnin());
            }

            optionPanel.addComponentWithLabel("Most recently sampled date:", mrsdField);
            mrsdField.setValue(inputFile.getMostRecentSampleDate());
        }

    }

    public int showDialog(InputFile inputFile) {

        this.inputFile = inputFile;

        setupPanel(inputFile);

        final JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Edit Input File Settings");
        dialog.pack();


        int result = JOptionPane.CANCEL_OPTION;

        dialog.setVisible(true);

        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public InputFile getInputFile() {
        inputFile.setBurnin(burninField.getValue());
        inputFile.setMostRecentSampleDate(mrsdField.getValue());
        return inputFile;
    }

}