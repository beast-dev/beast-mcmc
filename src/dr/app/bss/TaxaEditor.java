package dr.app.bss;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;

public class TaxaEditor {

	private PartitionDataList dataList;
	private MainFrame frame;
	
	// Window
	private JDialog window;
	private Frame owner;

	// Menubar
	private JMenuBar menu;

	// Buttons with options
	private JButton load;
	private JButton save;
	private JButton done;

	// Strings for paths
	private File workingDirectory;

	// Data, model & stuff for JTable
	private JTable table;
	private TableModel tableModel;
	private String[] COLUMN_NAMES = { "Name", "Height"};

	public TaxaEditor(
			MainFrame frame, PartitionDataList dataList
			) {

		this.frame = frame;
		this.dataList = dataList;
		
		// Setup Main Menu buttons
		load = new JButton("Load", Utils.createImageIcon(Utils.TEXT_FILE_ICON));
		save = new JButton("Save", Utils.createImageIcon(Utils.SAVE_ICON));
		done = new JButton("Done", Utils.createImageIcon(Utils.CHECK_ICON));

		// Add Main Menu buttons listeners
		load.addActionListener(new ListenLoadTaxaFile());
		save.addActionListener(new ListenSaveTaxaFile());
		done.addActionListener(new ListenOk());

		// Setup menu
		menu = new JMenuBar();
		menu.setLayout(new BorderLayout());
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);
		buttonsHolder.add(load);
		buttonsHolder.add(save);
		buttonsHolder.add(done);
		menu.add(buttonsHolder, BorderLayout.WEST);

		// Setup table
		table = new JTable();
		//TODO: override with own
		tableModel = new DefaultTableModel();
		table.setModel(tableModel);
		table.setSurrendersFocusOnKeystroke(true);


		JScrollPane scrollPane = new JScrollPane(table,
				ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS,
				ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		RowNumberTable rowNumberTable = new RowNumberTable(table);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER, rowNumberTable
				.getTableHeader());
		
		// Setup window
		owner = Utils.getActiveFrame();
		window = new JDialog(owner, "Setup location coordinates...");
		window.getContentPane().add(menu, BorderLayout.NORTH);
		window.getContentPane().add(scrollPane);
		window.pack();
		window.setLocationRelativeTo(owner);

	}// END: LocationCoordinatesEditor()

	private class ListenLoadTaxaFile implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			System.out.println("TODO");

		}// END: actionPerformed
	}// END: ListenOpenLocations

	private class ListenSaveTaxaFile implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

				System.out.println("TODO");

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

	private class ListenOk implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			System.out.println("TODO");
			window.setVisible(false);

			
		}// END: actionPerformed
	}// END: ListenSaveTaxaFile

	public void launch(File workingDirectory) {

		this.workingDirectory = workingDirectory;

		// Display Frame
		window.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		window.setSize(new Dimension(300, 300));
		window.setMinimumSize(new Dimension(100, 100));
		window.setResizable(true);
		window.setVisible(true);
	}// END: launch

}// END: class
