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
import javax.swing.SwingConstants;

import org.virion.jam.components.RealNumberField;

public class ClockRateModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;

	// Settings
	private OptionsPanel optionPanel;
	private JComboBox clockCombo;
	private RealNumberField[] clockParameterFields = new RealNumberField[PartitionData.clockParameterNames.length];

	// Buttons
	private JButton done;
	private JButton cancel;

	// Window
	private JDialog window;
	private Frame owner;

	public ClockRateModelEditor(PartitionDataList dataList, int row) {

		this.dataList = dataList;
		this.row = row;

		window = new JDialog(owner, "Setup clock rate model for partition "
				+ (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		clockCombo = new JComboBox();
		clockCombo.setOpaque(false);

		for (String clockModel : PartitionData.clockModels) {
			clockCombo.addItem(clockModel);
		}// END: fill loop

		clockCombo.addItemListener(new ListenClockCombo());

		for (int i = 0; i < PartitionData.clockParameterNames.length; i++) {
			clockParameterFields[i] = new RealNumberField();
			clockParameterFields[i].setColumns(8);
			clockParameterFields[i]
					.setValue(dataList.get(row).clockParameterValues[i]);
		}// END: fill loop

		setClockArguments();

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
		window.getContentPane().add(optionPanel, BorderLayout.CENTER);
		window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
		window.pack();
	}// END: Constructor

	private void setClockArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Clock rate model:"), clockCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = clockCombo.getSelectedIndex();

		for (int i = 0; i < dataList.get(row).clockParameterIndices[index].length; i++) {

			int k = dataList.get(row).clockParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(clockParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					PartitionData.clockParameterNames[k] + ":", panel);

		}// END: indices loop

		window.validate();
		window.repaint();
	}// END: setClockArguments

	private class ListenClockCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setClockArguments();

		}// END: actionPerformed
	}// END: ListenClockCombo

	public void collectSettings() {

		dataList.get(row).clockModel = clockCombo.getSelectedIndex();
		for (int i = 0; i < PartitionData.clockParameterNames.length; i++) {

			dataList.get(row).clockParameterValues[i] = clockParameterFields[i]
					.getValue();

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
	}// END: launch

}// END: class
