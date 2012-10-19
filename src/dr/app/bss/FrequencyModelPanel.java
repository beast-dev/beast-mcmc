package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.virion.jam.components.RealNumberField;

@SuppressWarnings("serial")
public class FrequencyModelPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame = null;
	private ArrayList<PartitionData> dataList = null;
	
	private OptionsPanel optionPanel;
	private JScrollPane scrollPane;
	private JComboBox frequencyCombo;
	private RealNumberField[] frequencyParameterFields = new RealNumberField[PartitionData.frequencyParameterNames.length];
	
	public FrequencyModelPanel(final BeagleSequenceSimulatorFrame frame,
			final ArrayList<PartitionData> dataList) {
		
		this.frame = frame;
		this.dataList = dataList;

		setOpaque(false);
		setLayout(new BorderLayout());

		scrollPane = new JScrollPane();
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
        optionPanel.setOpaque(false);
        scrollPane = new JScrollPane(optionPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.getViewport().setOpaque(false);

        add(scrollPane, BorderLayout.CENTER);
		
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
	}// END: setFrequencyArguments

	private class ListenFrequencyCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setFrequencyArguments();
			frame.fireModelChanged();

		}// END: actionPerformed
	}// END: ListenClockCombo

	public void collectSettings() {

		dataList.get(0).frequencyModel = frequencyCombo.getSelectedIndex();
		for (int i = 0; i < PartitionData.frequencyParameterNames.length; i++) {

			dataList.get(0).frequencyParameterValues[i] = frequencyParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	public JComponent getExportableComponent() {
		return this;
	}//END: getExportableComponent
	
}//END: class
