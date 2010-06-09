package dr.app.tracer.traces;

import dr.inference.trace.Filter;
import dr.inference.trace.FilteredTraceList;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceFactory;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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

    JList fileList;
    JTextField typeField;
    JTextField nameField;


    FilterPanel filterPanel;
    OptionsPanel tittlePanel;

    public FilterDialog(JFrame frame) {
        this.frame = frame;

        fileList = new JList();
        fileList.setVisibleRowCount(6);
        fileList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
//        fileList.addItemListener(new ItemListener() {
//            public void itemStateChanged(ItemEvent e) {
//                if (filteredTraceListGroup != null && fileList.getSelectedIndex() > 0)
//                    initComponents(filteredTraceListGroup.get(fileList.getSelectedIndex()));
//            }
//        });
        typeField = new JTextField();
        typeField.setColumns(20);
        typeField.setEditable(false);
        nameField = new JTextField();
        nameField.setColumns(30);
        nameField.setEditable(false);

        tittlePanel = new OptionsPanel(12, 12);
        tittlePanel.addComponentWithLabel("Selected File(s) : ", fileList);
        tittlePanel.addComponentWithLabel("Trace Name : ", nameField);
        tittlePanel.addComponentWithLabel("Trace Type : ", typeField);


    }

    public String showDialog(String traceName, List<FilteredTraceList> filteredTraceListGroup, String previousMessage) {
        this.filteredTraceListGroup = filteredTraceListGroup;
        this.traceName = traceName;

        Object[] fileNames = new Object[filteredTraceListGroup.size()];
        int[] indices = new int[filteredTraceListGroup.size()];
        for (int i = 0; i < fileNames.length; i++) {
            fileNames[i] = filteredTraceListGroup.get(i).getName();
            indices[i] = i;
        }

        fileList.setListData(fileNames);
        fileList.setSelectedIndices(indices);

//        initComponents(filteredTraceListGroup.get(0));
        FilteredTraceList filteredTraceList = filteredTraceListGroup.get(fileList.getSelectedIndex()); // only pick up the 1st one
        TraceDistribution td = filteredTraceList.getDistributionStatistics(filteredTraceList.getTraceIndex(traceName));

        typeField.setText(td.getTraceType().toString());
        nameField.setText(traceName);

        Filter f = filteredTraceList.getFilter(traceName);

        String[] sel;
        if (f == null) {
            sel = null;
        } else {
            sel = f.getIn();
        }
        if (td.getTraceType() == TraceFactory.TraceType.CONTINUOUS) {
            String[] minMax = new String[]{"0.1", "10"};  // todo
            filterPanel = new FilterContinuousPanel(minMax, sel);
        } else {// integer and string
            String[] all = td.getRangeAll();
            filterPanel = new FilterDiscretePanel(all, sel);
        }


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

        Object result = optionPane.getValue();

//        FilteredTraceList filteredTraceList = filteredTraceListGroup.get(treeFileCombo.getSelectedIndex());
//        TraceDistribution td = filteredTraceList.getDistributionStatistics(filteredTraceList.getTraceIndex(traceName));

        String message = "";
        if (result.equals(options[0])) {
            for (int i = 0; i < filteredTraceListGroup.size(); i++) {
                f = filteredTraceListGroup.get(i).getFilter(traceName);

                if (f == null) {
                    f = new Filter(traceName, td.getTraceType(), filterPanel.getSelectedValues(), td.getValues());
                } else {
                    f.setIn(filterPanel.getSelectedValues(), td.getValues());
                }

                filteredTraceListGroup.get(i).setFilter(f);
                filteredTraceListGroup.get(i).createTraceFilter(f);
            }
            message = f.getStatusMessage(); // todo

        } else if (result.equals(options[1])) {
            for (int i = 0; i < filteredTraceListGroup.size(); i++) {
                filteredTraceListGroup.get(i).removeFilter(traceName);
            }
            message = "";

        } else if (result.equals(options[2])) {
            message = previousMessage;
        }

        return message;
    }

//    private void initComponents(FilteredTraceList filteredTraceList) {
//
//        TraceDistribution td = filteredTraceList.getDistributionStatistics(filteredTraceList.getTraceIndex(traceName));
//
//        typeField.setText(td.getTraceType().toString());
//        nameField.setText(traceName);
//
//        Filter f = filteredTraceList.getFilter(traceName);
//
//        if (td.getTraceType() == TraceFactory.TraceType.CONTINUOUS) {
//
//        } else {// integer and string
//            String[] all = td.getRangeAll();
//            String[] sel;
//
//            if (f == null) {
//                sel = null;
//            } else {
//                sel = f.getIn();
//            }
//
//            filterPanel = new FilterDiscretePanel(all, sel);
//        }
//
//
//    }

    public String getName() {
        return nameField.getText();
    }

    abstract class FilterPanel extends JPanel {
        abstract Object[] getSelectedValues();
    }

    class FilterDiscretePanel extends FilterPanel {
        JList allValues;
        JList selectedValues;
//        JButton selectButton;

        FilterDiscretePanel(String[] allValuesArray, String[] selectedValuesArray) {
            setLayout(new FlowLayout());

            allValues = new JList(allValuesArray);
            allValues.setVisibleRowCount(6);
            allValues.setFixedCellWidth(100);
            allValues.setFixedCellHeight(15);
            allValues.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
            add(new JScrollPane(allValues));

            if (selectedValuesArray != null) {
                int[] indices = new int[selectedValuesArray.length];
                for (int i = 0; i < indices.length; i++) {
                    for (int j = 0; j < allValuesArray.length; j++) {
                        if (selectedValuesArray[i].equals(allValuesArray[j])) {
                            indices[i] = j;
                            break;
                        }
                    }
                }

                allValues.setSelectedIndices(indices);
            }

            JButton selectButton = new JButton("Select >>>");
            selectButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    selectedValues.setListData(allValues.getSelectedValues());
                }
            });
            add(selectButton);

            if (selectedValuesArray == null) {
                selectedValues = new JList();
            } else {
                selectedValues = new JList(selectedValuesArray);
            }
            selectedValues.setVisibleRowCount(6);
            selectedValues.setFixedCellWidth(100);
            selectedValues.setFixedCellHeight(15);
            selectedValues.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
            add(new JScrollPane(selectedValues));
        }

        public Object[] getSelectedValues() {
            return allValues.getSelectedValues();
        }

    }

    class FilterContinuousPanel extends FilterPanel {
        JTextField minField;
        JTextField maxField;

        FilterContinuousPanel(String[] minMax, String[] bound) {
            setLayout(new GridLayout(2, 3, 1, 10)); // 2 by 3, gap 5 by 1

            if (bound == null) {
               bound = new String[2];
            }

            minField = new JTextField(bound[0]);
            minField.setColumns(20);
            add(new JLabel("Set Minimum for Selecting Values : "));
            add(minField);
            add(new JLabel(", which should > " + minMax[0]));

            maxField = new JTextField(bound[1]);
            maxField.setColumns(20);
            add(new JLabel("Set Maximum for Selecting Values : "));
            add(maxField);
            add(new JLabel(", which should < " + minMax[1]));
        }

        public Object[] getSelectedValues() {
            return new String[]{minField.getText(), maxField.getText()};
        }

    }

}