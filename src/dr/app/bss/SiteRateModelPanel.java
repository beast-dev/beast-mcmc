package dr.app.bss;

import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.virion.jam.components.RealNumberField;

@SuppressWarnings("serial")
public class SiteRateModelPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JComboBox siteCombo;
	private RealNumberField[] siteParameterFields = new RealNumberField[BeagleSequenceSimulatorData.siteParameterNames.length];

	public SiteRateModelPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {

		this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		siteCombo = new JComboBox();
		siteCombo.setOpaque(false);

		for (String siteModel : BeagleSequenceSimulatorData.siteModels) {
			siteCombo.addItem(siteModel);
		}// END: fill loop

		siteCombo.addItemListener(new ListenSiteCombo());

		for (int i = 0; i < BeagleSequenceSimulatorData.siteParameterNames.length; i++) {
			siteParameterFields[i] = new RealNumberField();
			siteParameterFields[i].setColumns(8);
			siteParameterFields[i].setValue(data.siteParameterValues[i]);
		}// END: fill loop

		setSiteArguments();

	}// END: Constructor

	private void setSiteArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Site Rate model:"), siteCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = siteCombo.getSelectedIndex();

		for (int i = 0; i < data.siteParameterIndices[index].length; i++) {

			int k = data.siteParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(siteParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					BeagleSequenceSimulatorData.siteParameterNames[k] + ":",
					panel);

		}// END: indices loop

//		validate();
//		repaint();
	}// END: setSiteArguments

	private class ListenSiteCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setSiteArguments();
			frame.fireModelChanged();

		}// END: actionPerformed
	}// END: ListenSiteCombo

	public void collectSettings() {

		data.siteModel = siteCombo.getSelectedIndex();
		for (int i = 0; i < BeagleSequenceSimulatorData.siteParameterNames.length; i++) {

			data.siteParameterValues[i] = siteParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
