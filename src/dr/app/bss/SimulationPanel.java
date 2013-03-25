package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import dr.app.gui.components.WholeNumberField;

@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

	private MainFrame frame;
	private PartitionDataList dataList;
	private OptionsPanel optionPanel;

	private WholeNumberField simulationsNumberField;

	// Buttons
	private JButton simulate;
	private JButton generateXML;

	public SimulationPanel(final MainFrame frame,
			final PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		//TODO: spinner?
		simulationsNumberField = new WholeNumberField(1, Integer.MAX_VALUE);
		simulationsNumberField.setColumns(8);
		simulationsNumberField.setValue(dataList.simulationsCount);
		optionPanel.addComponentWithLabel("Number of simulations:", simulationsNumberField);
		
		// Buttons holder
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);

		// simulate button
		simulate = new JButton("Simulate",
				BeagleSequenceSimulatorApp.biohazardIcon);
		simulate.addActionListener(new ListenSimulate());
		buttonsHolder.add(simulate);

		generateXML = new JButton("Generate XML",
				BeagleSequenceSimulatorApp.hammerIcon);
		generateXML.addActionListener(new ListenGenerateXML());
		buttonsHolder.add(generateXML);

		setOpaque(false);
		setLayout(new BorderLayout());
		add(optionPanel, BorderLayout.NORTH);
		add(buttonsHolder, BorderLayout.SOUTH);

	}// END: SimulationPanel

	public final void collectSettings() {
		
		dataList.simulationsCount = simulationsNumberField.getValue();
		
	}// END: collectSettings

	private class ListenSimulate implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frame.doExport();

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

	private class ListenGenerateXML implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frame.doGenerateXML();

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

	public void setBusy() {
		simulate.setEnabled(false);
		generateXML.setEnabled(false);
	}// END: setBusy

	public void setIdle() {
		simulate.setEnabled(true);
		generateXML.setEnabled(true);
	}// END: setIdle

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
