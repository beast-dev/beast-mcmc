package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import dr.app.gui.components.WholeNumberField;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

	private MainFrame frame;
	private PartitionDataList dataList;
	private OptionsPanel optionPanel;

	private WholeNumberField simulationsNumberField;
	private WholeNumberField startingSeedNumberField;
	
	// Buttons
	private JButton simulate;
	private JButton generateXML;

	// Check boxes
	private JCheckBox setSeed;
	
	public SimulationPanel(final MainFrame frame,
			final PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		simulationsNumberField = new WholeNumberField(1, Integer.MAX_VALUE);
		simulationsNumberField.setColumns(10);
		simulationsNumberField.setValue(dataList.simulationsCount);
		optionPanel.addComponentWithLabel("Number of simulations:", simulationsNumberField);
		
		setSeed = new JCheckBox();
		setSeed.addItemListener(new CheckBoxListener());
		setSeed.setSelected(dataList.setSeed);
		optionPanel.addComponentWithLabel("Set seed:", setSeed);
		
		startingSeedNumberField = new WholeNumberField(1, Long.MAX_VALUE);
		startingSeedNumberField.setColumns(10);
		startingSeedNumberField.setValue(dataList.startingSeed);
		startingSeedNumberField.setEnabled(dataList.setSeed);
		optionPanel.addComponentWithLabel("Starting seed:", startingSeedNumberField);
		
		// Buttons holder
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);

		// simulate button
		simulate = new JButton("Simulate",
				Utils.createImageIcon(Utils.BIOHAZARD_ICON));
		simulate.addActionListener(new ListenSimulate());
		buttonsHolder.add(simulate);

		generateXML = new JButton("Generate XML",
				Utils.createImageIcon(Utils.HAMMER_ICON));
		generateXML.addActionListener(new ListenGenerateXML());
		buttonsHolder.add(generateXML);

		setOpaque(false);
		setLayout(new BorderLayout());
		add(optionPanel, BorderLayout.NORTH);
		add(buttonsHolder, BorderLayout.SOUTH);

	}// END: SimulationPanel

	public final void collectSettings() {

		dataList.simulationsCount = simulationsNumberField.getValue();
		if (dataList.setSeed) {
			dataList.startingSeed = startingSeedNumberField.getValue();
		}

	}// END: collectSettings

	private class CheckBoxListener implements ItemListener {
		public void itemStateChanged(ItemEvent e) {

			if (setSeed.isSelected()) {
				startingSeedNumberField.setEnabled(true);
				dataList.setSeed = true;
			} else {
				startingSeedNumberField.setEnabled(false);
				dataList.setSeed = false;
			}

		}
	}// END: CheckBoxListener
	
	private class ListenSimulate implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frame.doExport();

		}// END: actionPerformed
	}// END: ListenSimulate

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
