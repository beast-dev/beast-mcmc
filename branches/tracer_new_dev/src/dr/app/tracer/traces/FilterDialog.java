package dr.app.tracer.traces;

import dr.inference.trace.TraceDistribution;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: FrequencyPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class FilterDialog {
    private JFrame frame;

    JComboBox typeCombo;
    JTextField nameField;

    private OptionsPanel filterPanel;

    public FilterDialog(JFrame frame) {
        this.frame = frame;

        typeCombo = new JComboBox();
        nameField = new JTextField();
        nameField.setColumns(20);
        nameField.setEditable(false);
        typeCombo.setEditable(false);
        

    }

    public int showDialog(String traceName, TraceDistribution td) {

        typeCombo.removeAllItems();
        typeCombo.addItem(td.getTraceType());
        nameField.setText(traceName);

        OptionsPanel tittlePanel = new OptionsPanel(12, 12);
        tittlePanel.addComponentWithLabel("Trace Name : ", nameField);
        tittlePanel.addComponentWithLabel("Trace Type : ", typeCombo);        

        filterPanel = new OptionsPanel(12, 12);


        

        JPanel basePanel = new JPanel(new BorderLayout(0, 0));
        basePanel.add(tittlePanel, BorderLayout.NORTH);
        basePanel.add(filterPanel, BorderLayout.CENTER);

Object[] options = {"Apply Filter", "Remove Filter", "Cancel"};
       JOptionPane optionPane = new JOptionPane(basePanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.YES_NO_CANCEL_OPTION,
                null,
                options,
                options[2]);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create A Filter");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public String getName() {
        return nameField.getText();
    }

}