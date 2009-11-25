package dr.app.phylogeography.spread;

import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.components.RealNumberField;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class InputFileSettingsDialog {

    private JFrame frame;

    private RealNumberField mrsdField;

    private OptionsPanel optionPanel;

    private InputFile inputFile = null;

    public InputFileSettingsDialog(JFrame frame) {
        this.frame = frame;

        mrsdField = new RealNumberField(0, Double.POSITIVE_INFINITY);
        mrsdField.setColumns(20);

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel(InputFile inputFile) {
        optionPanel.removeAll();

        if (inputFile.getTree() != null) {
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
        inputFile.setMostRecentSampleDate(mrsdField.getValue());
        return inputFile;
    }

}