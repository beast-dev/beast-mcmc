package dr.app.bss;

import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import dr.app.gui.components.RealNumberField;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class BranchSubstitutionModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;
	
	// Settings
	private OptionsPanel optionPanel;
	private DisabledItemsComboBox substitutionCombo;
	private RealNumberField[] substitutionParameterFields;

	//Buttons
	private JButton done;
	private JButton cancel;
	
	// Window
	private JDialog window;
	private Frame owner;
	
	public BranchSubstitutionModelEditor(PartitionDataList dataList, int row) {
		
		this.dataList = dataList;
		this.row = row;
		
		substitutionParameterFields = new RealNumberField[PartitionData.substitutionParameterNames.length];
		window = new JDialog(owner, "Setup substitution model for partition " + (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		
		substitutionCombo = new DisabledItemsComboBox();
		substitutionCombo.setOpaque(false);

		int indexOf = 0;
		for (String substitutionModel : PartitionData.substitutionModels) {

			if (PartitionData.substitutionCompatibleDataTypes[indexOf] == dataList.get(row).dataTypeIndex) {

				substitutionCombo.addItem(substitutionModel, false);

			} else {

				substitutionCombo.addItem(substitutionModel, true);

			}

			indexOf++;
		}// END: fill loop

		substitutionCombo.addItemListener(new ListenSubstitutionCombo());
		
		for (int i = 0; i < PartitionData.substitutionParameterNames.length; i++) {
			
			substitutionParameterFields[i] = new RealNumberField(0.0, Double.MAX_VALUE);
			substitutionParameterFields[i].setColumns(8);
			substitutionParameterFields[i].setValue(this.dataList.get(row).substitutionParameterValues[i]);
			
		}// END: fill loop

		setSubstitutionArguments();

		// Buttons
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);
		
		cancel = new JButton("Cancel", Utils.createImageIcon(Utils.CLOSE_ICON));
		cancel.addActionListener(new ListenCancel());
		buttonsHolder.add(cancel);
		
		done = new JButton("Done", Utils.createImageIcon(Utils.CHECK_ICON));
		done.addActionListener(new ListenOk());
		buttonsHolder.add(done);
		
		// Window
		owner = Utils.getActiveFrame();
		window.setLocationRelativeTo(owner);
		window.getContentPane().setLayout(new BorderLayout());
		window.getContentPane().add(optionPanel, BorderLayout.CENTER);
		window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
		window.pack();
		
		//return to the previously chosen index on start
		substitutionCombo.setSelectedIndex(dataList.get(row).substitutionModelIndex);
		
	}//END: Constructor
	
	private void setSubstitutionArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Substitution model:"), substitutionCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = substitutionCombo.getSelectedIndex();

		for (int i = 0; i < PartitionData.substitutionParameterIndices[index].length; i++) {

			int k = PartitionData.substitutionParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(substitutionParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(PartitionData.substitutionParameterNames[k] + ":", panel);

		}// END: indices loop

		window.validate();
		window.repaint();
	}// END: setSubstitutionArguments
	
	private class ListenSubstitutionCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setSubstitutionArguments();

		}// END: actionPerformed
	}// END: ListenSubstitutionCombo
	
	public void collectSettings() {

		dataList.get(row).substitutionModelIndex = substitutionCombo.getSelectedIndex();
		for (int i = 0; i < PartitionData.substitutionParameterNames.length; i++) {

			dataList.get(row).substitutionParameterValues[i] = substitutionParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	private class ListenOk implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			collectSettings();
			
		}// END: actionPerformed
	}// END: ListenOk
	
	private class ListenCancel implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			
		}// END: actionPerformed
	}// END: ListenCancel
	
	public void showWindow() {
		window.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		window.setSize(new Dimension(450, 400));
		window.setMinimumSize(new Dimension(100, 100));
		window.setResizable(true);
		window.setModal(true);
		window.setVisible(true);
	}//END: showWindow
	
	public void launch() {

		if (SwingUtilities.isEventDispatchThread()) {
			showWindow();
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showWindow();
				}
			});
		}// END: edt check

	}// END: launch
	
}// END: class
