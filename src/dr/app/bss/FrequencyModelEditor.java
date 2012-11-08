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
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.virion.jam.components.RealNumberField;

public class FrequencyModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;
	
	// Settings
	private OptionsPanel optionPanel;
	private JScrollPane scrollPane;
	private JComboBox frequencyCombo;
	private RealNumberField[] frequencyParameterFields = new RealNumberField[PartitionData.frequencyParameterNames.length];
	
	//Buttons
	private JButton done;
	private JButton cancel;
	
	// Window
	private JDialog window;
	private Frame owner;
	
	public FrequencyModelEditor(PartitionDataList dataList, int row) {
		
		this.dataList = dataList;
		this.row = row;

		window = new JDialog(owner, "Setup branch substitution model for partition " + (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		
		scrollPane = new JScrollPane();
        optionPanel.setOpaque(false);
        scrollPane = new JScrollPane(optionPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.getViewport().setOpaque(false);

		frequencyCombo = new JComboBox();

		for (String frequencyModel : PartitionData.frequencyModels) {
			frequencyCombo.addItem(frequencyModel);
		}// END: fill loop

		frequencyCombo.addItemListener(new ListenFrequencyCombo());

		for (int i = 0; i < PartitionData.frequencyParameterNames.length; i++) {
			frequencyParameterFields[i] = new RealNumberField();
			frequencyParameterFields[i].setColumns(8);
			frequencyParameterFields[i].setValue(dataList.get(0).frequencyParameterValues[i]);
		}// END: fill loop

		setFrequencyArguments();
		
		// Buttons
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);
		
		cancel = new JButton("Cancel", BeagleSequenceSimulatorApp.closeIcon);
		cancel.addActionListener(new ListenCancel());
		buttonsHolder.add(cancel);
		
		done = new JButton("Done", BeagleSequenceSimulatorApp.doneIcon);
		done.addActionListener(new ListenOk());
		buttonsHolder.add(done);
		
		// Window
		owner = Utils.getActiveFrame();
		window.setLocationRelativeTo(owner);
		window.getContentPane().setLayout(new BorderLayout());
		window.getContentPane().add(scrollPane, BorderLayout.CENTER);
		window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
		window.pack();
	}//END: Constructor

	private void setFrequencyArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Frequency model:"), frequencyCombo);
		
		optionPanel.addSeparator();
		optionPanel.addComponentWithLabel("Set parameter values:", new JLabel());

		int index = frequencyCombo.getSelectedIndex();

		for (int i = 0; i < dataList.get(0).frequencyParameterIndices[index].length; i++) {

			int k = dataList.get(0).frequencyParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(frequencyParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					PartitionData.frequencyParameterNames[k] + ":",
					panel);

		}// END: indices loop
		
		window.validate();
		window.repaint();
	}// END: setFrequencyArguments

	private class ListenFrequencyCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setFrequencyArguments();

		}// END: actionPerformed
	}// END: ListenClockCombo

	public void collectSettings() {

		dataList.get(row).frequencyModel = frequencyCombo.getSelectedIndex();
		for (int i = 0; i < PartitionData.frequencyParameterNames.length; i++) {

			dataList.get(row).frequencyParameterValues[i] = frequencyParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	private class ListenOk implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			collectSettings();
			
		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates
	
	private class ListenCancel implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			
		}// END: actionPerformed
	}// END: ListenCancel
	
	public void launch() {
		window.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		window.setSize(new Dimension(450, 400));
		window.setMinimumSize(new Dimension(100, 100));
		window.setResizable(true);
		window.setVisible(true);
	}//END: launch
	
}//END: class
