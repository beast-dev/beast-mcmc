package dr.app.bss;

import java.awt.BorderLayout;

import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ScrollPaneConstants;

@SuppressWarnings("serial")
public class TerminalPanel extends JPanel {

//	private MainFrame frame;
//	private PartitionDataList dataList;
	
	private JTextArea textArea;

	public TerminalPanel(
//			final MainFrame frame,
//			final PartitionDataList dataList
			) {

//		this.frame = frame;
//		this.dataList = dataList;
		
		// Setup miscallenous
		setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));

		// Setup text area
		textArea = new JTextArea();
		textArea.setEditable(true);
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
		// textAreaTest.setContentType("text/html");
		// textAreaTest.setText("kutas");

//		ScrollPane scrollPane = new ScrollPane(ScrollPane.SCROLLBARS_ALWAYS);
//		scrollPane.add(textArea);
//		add(scrollPane, BorderLayout.CENTER);
		
		JScrollPane scrollPane = new JScrollPane(textArea,
				ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS,
				ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		add(scrollPane, BorderLayout.CENTER);

	}//END: Constructor

	public void setText(String text) {
		textArea.append(text);
	}

	public void clearTerminal() {
		textArea.setText("");
	}

	public void setTextAndClar(String text) {
		setText(text);
		clearTerminal();
	}

}// END: class
