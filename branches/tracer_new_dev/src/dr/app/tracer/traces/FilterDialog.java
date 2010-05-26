package dr.app.tracer.traces;

import dr.inference.trace.FilteredTraceList;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceFactory;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: FrequencyPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class FilterDialog {
    private JFrame frame;
    List<FilteredTraceList> filteredTraceListGroup;
    String traceName;

    JComboBox treeFileCombo;
    JComboBox typeCombo;
    JTextField nameField;


    JPanel filterPanel;
    OptionsPanel tittlePanel;

    public FilterDialog(JFrame frame) {
        this.frame = frame;

        treeFileCombo = new JComboBox();
        treeFileCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                if (filteredTraceListGroup != null)
                    initComponents(filteredTraceListGroup.get(treeFileCombo.getSelectedIndex()));
            }
        });
        typeCombo = new JComboBox();
        nameField = new JTextField();
        nameField.setColumns(20);
        nameField.setEditable(false);
        typeCombo.setEditable(false);

        tittlePanel = new OptionsPanel(12, 12);
        tittlePanel.addComponentWithLabel("Tree File : ", treeFileCombo);
        tittlePanel.addComponentWithLabel("Trace Name : ", nameField);
        tittlePanel.addComponentWithLabel("Trace Type : ", typeCombo);


    }

    public int showDialog(String traceName, List<FilteredTraceList> filteredTraceListGroup) {
        this.filteredTraceListGroup = filteredTraceListGroup;
        this.traceName = traceName;

        treeFileCombo.removeAllItems();
        for (FilteredTraceList treeFile : filteredTraceListGroup) {
            treeFileCombo.addItem(treeFile.getName());
        }

        treeFileCombo.setSelectedIndex(0);

        initComponents(filteredTraceListGroup.get(0));

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

        int result = ((Integer) optionPane.getValue()).intValue();
        
        return result;
    }

    private void initComponents(FilteredTraceList filteredTraceList) {

        TraceDistribution td = filteredTraceList.getDistributionStatistics(filteredTraceList.getTraceIndex(traceName));

        typeCombo.removeAllItems();
        typeCombo.addItem(td.getTraceType());
        nameField.setText(traceName);


        if (td.getTraceType() == TraceFactory.TraceType.CONTINUOUS) {

        } else {// integer and string
            List vl = td.credSet.getValues();
            String[] sl = new String[vl.size()];
            for (int i = 0; i < vl.size(); i++) {
                sl[i] = vl.get(i).toString();
            }

            filterPanel = new FilterDiscretePanel(sl);
        }


    }

    public String getName() {
        return nameField.getText();
    }

    class FilterDiscretePanel extends JPanel {
        JList allValues;
        JList selectedValues;
        JButton selectButton;

        FilterDiscretePanel(String[] allValuesArray) {
            setLayout(new FlowLayout());

            allValues = new JList(allValuesArray);
            allValues.setVisibleRowCount(6);
            allValues.setFixedCellWidth(100);
            allValues.setFixedCellHeight(15);
            allValues.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
            add(new JScrollPane(allValues));

            selectButton = new JButton("Select >>>");
            selectButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    selectedValues.setListData(allValues.getSelectedValues());
                }
            });
            add(selectButton);

            selectedValues = new JList();
            selectedValues.setVisibleRowCount(6);
            selectedValues.setFixedCellWidth(100);
            selectedValues.setFixedCellHeight(15);
            selectedValues.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
            add(new JScrollPane(selectedValues));
        }

    }

}