package dr.app.phylogeography.spread.layerspanel;

import dr.app.gui.DeleteActionResponder;
import dr.app.phylogeography.builder.*;
import dr.app.phylogeography.spread.*;
import dr.app.gui.table.MultiLineTableCellRenderer;
import jam.framework.Exportable;
import jam.panels.ActionPanel;
import dr.app.gui.table.TableEditorStopper;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.dnd.*;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class LayersPanel extends JPanel implements DeleteActionResponder, Exportable {
    public final static BuilderFactory[] builderFactories = {
            DiscreteDiffusionTreeBuilder.FACTORY,
            ContinuousDiffusionTreeBuilder.FACTORY
    };

    private JTable layerTable = null;
    private LayerTableModel layerTableModel = null;

    private SpreadFrame frame = null;

    private CreateBuilderDialog createBuilderDialog = null;
    private EditBuilderDialog editBuilderDialog = null;

    private final SpreadDocument document;

    public LayersPanel(final SpreadFrame parent, final SpreadDocument document) {

        this.frame = parent;
        this.document = document;

        layerTableModel = new LayerTableModel();
        layerTable = new JTable(layerTableModel);

        layerTable.getTableHeader().setReorderingAllowed(false);
//        layerTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        TableColumn col = layerTable.getColumnModel().getColumn(0);
        col.setCellRenderer(new MultiLineTableCellRenderer());

        layerTable.setRowHeight(layerTable.getRowHeight() * 3);

        layerTable.setDragEnabled(true);

//        layerTable.setDropMode(DropMode.INSERT);
//        layerTable.setTransferHandler(new MyListDropHandler(layerTable));

        new MyDragListener(layerTable);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(layerTable);

        layerTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        layerTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    editSelection();
                }
            }
        });

        JScrollPane scrollPane = new JScrollPane(layerTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

//        JToolBar toolBar1 = new JToolBar();
//        toolBar1.setFloatable(false);
//        toolBar1.setOpaque(false);
//        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));

//        JButton button = new JButton(unlinkModelsAction);
//        unlinkModelsAction.setEnabled(false);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);

        ActionPanel actionPanel1 = new ActionPanel(true);
        actionPanel1.setAddAction(addAction);
        actionPanel1.setRemoveAction(removeAction);
        actionPanel1.setActionAction(editAction);

        removeAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new BorderLayout());
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1, BorderLayout.WEST);

        buildAction.setEnabled(false);
        JButton generateButton = new JButton(buildAction);
        generateButton.putClientProperty("JButton.buttonType", "roundRect");
        controlPanel1.add(generateButton, BorderLayout.EAST);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));
//        add(toolBar1, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);
        add(controlPanel1, BorderLayout.SOUTH);

        document.addListener(new SpreadDocument.Listener() {
            public void dataChanged() {
            }

            public void settingsChanged() {
                LayersPanel.this.settingsChanged();
            }
        });
    }

    public void selectionChanged() {
        int[] selRows = layerTable.getSelectedRows();
        boolean hasSelection = (selRows != null && selRows.length != 0);
        removeAction.setEnabled(hasSelection);

        frame.setRemoveActionEnabled(this, hasSelection);
    }

    public JComponent getExportableComponent() {
        return layerTable;
    }

    public void delete() {
        int[] selRows = layerTable.getSelectedRows();
        Set<Builder> buildersToRemove = new HashSet<Builder>();
        for (int row : selRows) {
            buildersToRemove.add(document.getLayerBuilders().get(row));
        }
//
//        // TODO: would probably be a good idea to check if the user wants to remove the last layer
        document.getLayerBuilders().removeAll(buildersToRemove);
        document.fireSettingsChanged();
    }

    public Action getDeleteAction() {
        return removeAction;
    }

    public void editSelection() {
        int selRow = layerTable.getSelectedRow();
        if (selRow >= 0) {
            Builder builder = document.getLayerBuilders().get(selRow);
            editSettings(builder);
        }
    }

    private void editSettings(Builder builder) {
        if (editBuilderDialog == null) {
            editBuilderDialog = new EditBuilderDialog(frame);
        }

        int result = editBuilderDialog.showDialog(builder, document);

        if (result != JOptionPane.CANCEL_OPTION) {
            editBuilderDialog.getBuilder(); // force update of builder settings
            document.fireSettingsChanged();
        }
    }

    public void addLayer() {
        if (createBuilderDialog == null) {
            createBuilderDialog = new CreateBuilderDialog(frame);
        }

        int result = createBuilderDialog.showDialog(builderFactories, document);
        if (result != JOptionPane.CANCEL_OPTION) {
            Builder builder = createBuilderDialog.getBuilder();
            document.addLayerBuilder(builder);
            editSettings(builder);
        }
    }

    public void selectAll() {
        layerTable.selectAll();
    }

    private void buildAll() {
        for (Builder builder : document.getLayerBuilders()) {
            try {
                builder.build();
            } catch (BuildException be) {
                JOptionPane.showMessageDialog(frame, "For layer, " + builder.getName() + ": " + be.getMessage(),
                        "Error building layer",
                        JOptionPane.ERROR_MESSAGE);

            }
        }
        document.fireSettingsChanged();
    }

    private void settingsChanged() {
        boolean needsBuilding = false;
        for (Builder builder : document.getLayerBuilders()) {
            if (!builder.isBuilt()) {
                needsBuilding = true;
            }
        }
        buildAction.setEnabled(needsBuilding);
        layerTableModel.fireTableDataChanged();
    }

    class LayerTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        //        String[] columnNames = {"Name", "Layer Type", "Input File", "Built?"};
        String[] columnNames = {"Layers"};

        public LayerTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return document.getLayerBuilders().size();
        }

        public Object getValueAt(int row, int col) {
            Builder builder = document.getLayerBuilders().get(row);
            return builder;
//            switch (col) {
//                case 0:
//                    return builder.getName();
//                case 1:
//                    return builder.getBuilderName();
//                case 2:
//                    return builder.getInputFile();
//                case 3:
//                    return (builder.isBuilt() ? "Yes" : "No");
//                default:
//                    throw new IllegalArgumentException("unknown column, " + col);
//            }
        }

        public void setValueAt(Object aValue, int row, int col) {
//            Builder builder = document.getLayerBuilders().get(row);
//            switch (col) {
//                case 0:
//                    String name = ((String) aValue).trim();
//                    if (name.length() > 0) {
//                        builder.setName(name);
//                    }
//                    break;
//            }
//            document.fireSettingsChanged();
        }

        public boolean isCellEditable(int row, int col) {
            boolean editable = false;

//            switch (col) {
//                case 0:// name
//                    editable = true;
//                    break;
//                default:
//                    editable = false;
//            }
//
            return editable;
        }


        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            if (getRowCount() == 0) {
                return Object.class;
            }
            return getValueAt(0, c).getClass();
        }

        public String toString() {
            StringBuffer buffer = new StringBuffer();

            buffer.append(getColumnName(0));
            for (int j = 1; j < getColumnCount(); j++) {
                buffer.append("\t");
                buffer.append(getColumnName(j));
            }
            buffer.append("\n");

            for (int i = 0; i < getRowCount(); i++) {
                buffer.append(getValueAt(i, 0));
                for (int j = 1; j < getColumnCount(); j++) {
                    buffer.append("\t");
                    buffer.append(getValueAt(i, j));
                }
                buffer.append("\n");
            }

            return buffer.toString();
        }
    }

    private final Action addAction = new AddLayerAction();
    private final Action removeAction = new RemoveLayerAction();
    private final Action editAction = new EditLayerAction();
    private final Action buildAction = new BuildAction();

    public class AddLayerAction extends AbstractAction {
        public AddLayerAction() {
            super("Add");
            setToolTipText("Use this button to create a new layer");
        }

        public void actionPerformed(ActionEvent ae) {
            addLayer();
        }
    }

    public class RemoveLayerAction extends AbstractAction {
        public RemoveLayerAction() {
            super("Remove");
            setToolTipText("Use this button to remove a selected layer from the table");
        }

        public void actionPerformed(ActionEvent ae) {
            delete();
        }
    }

    public class EditLayerAction extends AbstractAction {
        public EditLayerAction() {
            super("Edit");
            setToolTipText("Use this button to edit a selected layer in the table");
        }

        public void actionPerformed(ActionEvent ae) {
            editSelection();
        }
    }


    public class BuildAction extends AbstractAction {
        public BuildAction() {
            super("Build");
            setToolTipText("Use this button to build the layers");
        }

        public void actionPerformed(ActionEvent ae) {
            buildAll();
        }
    }

    class MyDragListener implements DragSourceListener, DragGestureListener {
        JTable table;

        DragSource ds = new DragSource();

        public MyDragListener(JTable table) {
            this.table = table;
            DragGestureRecognizer dgr = ds.createDefaultDragGestureRecognizer(table, DnDConstants.ACTION_MOVE, this);
        }

        public void dragGestureRecognized(DragGestureEvent dge) {
            StringSelection transferable = new StringSelection(Integer.toString(table.getSelectedRow()));
            ds.startDrag(dge, DragSource.DefaultCopyDrop, transferable, this);
        }

        public void dragEnter(DragSourceDragEvent dsde) {
        }

        public void dragExit(DragSourceEvent dse) {
        }

        public void dragOver(DragSourceDragEvent dsde) {
        }

        public void dragDropEnd(DragSourceDropEvent dsde) {
        }

        public void dropActionChanged(DragSourceDragEvent dsde) {
        }
    }

//    class MyListDropHandler extends TransferHandler {
//        JTable table;
//
//        public MyListDropHandler(JTable table) {
//            this.table = table;
//        }
//
//        public boolean canImport(TransferHandler.TransferSupport support) {
//            if (!support.isDataFlavorSupported(DataFlavor.stringFlavor)) {
//                return false;
//            }
//            JTable.DropLocation dl = (JTable.DropLocation) support.getDropLocation();
//            return dl.getRow() != -1;
//        }
//
//        public boolean importData(TransferHandler.TransferSupport support) {
//            if (!canImport(support)) {
//                return false;
//            }
//
//            Transferable transferable = support.getTransferable();
//            String indexString;
//            try {
//                indexString = (String) transferable.getTransferData(DataFlavor.stringFlavor);
//            } catch (Exception e) {
//                return false;
//            }
//
//            int index = Integer.parseInt(indexString);
//            Builder builder = document.getLayerBuilders().get(index);
//
//            JTable.DropLocation dl = (JTable.DropLocation) support.getDropLocation();
//            int dropTargetIndex = dl.getRow();
//
//            document.getLayerBuilders().add(dropTargetIndex, builder);
//            document.getLayerBuilders().remove(builder);
//
//            document.fireSettingsChanged();
//
//            return true;
//        }
//    }
}