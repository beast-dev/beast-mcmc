package dr.app.bss;

import jam.framework.Exportable;

import java.awt.BorderLayout;
import java.util.ArrayList;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

@SuppressWarnings("serial")
public class PartitionsPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private ArrayList<BeagleSequenceSimulatorData> dataList = null;

	private JTable partitionTable = null;
	private PartitionTableModel partitionTableModel = null;
	private JScrollPane scrollPane;

	public PartitionsPanel(final BeagleSequenceSimulatorFrame frame,
			final ArrayList<BeagleSequenceSimulatorData> dataList) {

		super();

		this.frame = frame;
		this.dataList = dataList;

		partitionTableModel = new PartitionTableModel(this.dataList);
		partitionTableModel
				.addTableModelListener(new PartitionTableModelListener());
		partitionTable = new JTable(partitionTableModel);

		setLayout(new BorderLayout());
		
		scrollPane = new JScrollPane(partitionTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		RowNumberTable rowNumberTable = new RowNumberTable(partitionTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER, rowNumberTable
				.getTableHeader());
		
		
		scrollPane.getViewport().setOpaque(false);

		
		
		add(scrollPane, BorderLayout.CENTER);

	}// END: Constructor

	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	private class PartitionTableModelListener implements TableModelListener {
		public void tableChanged(TableModelEvent ev) {

			System.out.println("TODO");

		}
	}// END: InteractiveTableModelListener

}// END: class
