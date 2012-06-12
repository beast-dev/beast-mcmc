package dr.app.bss;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.Box;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

public class BeagleSequenceSimulatorApp {

	public static final String VERSION = "0.0.1";
	public static final String DATE_STRING = "2012";
	private static final String BEAGLE_SEQUENCE_SIMULATOR = "Beagle Sequence Simulator";
	
	// Dimension
	private Dimension dimension;

	// Frame
	private JFrame frame;

	// Menubar
	private JMenuBar mainMenu;

	// Menus with items
	private FileMenu fileMenu;
	private EditMenu editMenu;
	private HelpMenu helpMenu;
	private JTextArea textArea;
	
	public BeagleSequenceSimulatorApp() throws ClassNotFoundException,
			InstantiationException, IllegalAccessException,
			UnsupportedLookAndFeelException {

		boolean lafLoaded = false;

		// Setup Look & Feel
		if (System.getProperty("os.name").toLowerCase().startsWith("mac os x")) {

			// Mac stuff
			System.setProperty("apple.awt.showGrowBox", "true");
			System.setProperty("apple.awt.brushMetalLook", "true");
			System.setProperty("apple.laf.useScreenMenuBar", "true");

			System.setProperty("apple.awt.graphics.UseQuartz", "true");
			System.setProperty("apple.awt.antialiasing", "true");
			System.setProperty("apple.awt.rendering", "VALUE_RENDER_QUALITY");

			System.setProperty("apple.laf.useScreenMenuBar", "true");
			System.setProperty("apple.awt.draggableWindowBackground", "true");
			System.setProperty("apple.awt.showGrowBox", "true");

			UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
			UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));

			try {

				// UIManager.setLookAndFeel(UIManager
				// .getSystemLookAndFeelClassName());

				UIManager.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
				lafLoaded = true;

			} catch (Exception e) {
				//
			}

		} else {

			try {

				// UIManager.setLookAndFeel(UIManager
				// .getSystemLookAndFeelClassName());
				
				UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
				lafLoaded = true;

			} catch (Exception e) {
				//
			}

		}

		if (!lafLoaded) {

			try {

				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

			} catch (Exception e1) {

				e1.printStackTrace();
				System.out.println("Specified l&f not found. Loading system default l&f");

			}
		}

		dimension = Toolkit.getDefaultToolkit().getScreenSize();
		Toolkit.getDefaultToolkit().setDynamicLayout(true);

		// Setup Main Frame
		frame = new JFrame(BEAGLE_SEQUENCE_SIMULATOR);
		frame.setLayout(new BorderLayout());
		frame.addWindowListener(new ListenCloseWdw());

		textArea = new JTextArea(frame.getWidth(), frame.getHeight());
		textArea.setEditable(true);
		JScrollPane scrollingText = new JScrollPane(textArea);
	
		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BorderLayout());
		mainPanel.add(scrollingText, BorderLayout.CENTER);
		frame.setLocationRelativeTo(null);
		frame.setContentPane(mainPanel);
		mainPanel.setOpaque(true);
		
		// Setup Main Menu items
		fileMenu = new FileMenu();
		editMenu = new EditMenu();
		helpMenu = new HelpMenu();

		// Setup Main Menu
		mainMenu = new JMenuBar();
//		mainMenu.setLayout(new BoxLayout(mainMenu, BoxLayout.PAGE_AXIS));
		mainMenu.add(fileMenu);
		mainMenu.add(editMenu);
		mainMenu.add(helpMenu);
		
		// Setup frame
		frame.setJMenuBar(mainMenu);
//		frame.add(mainMenu, BorderLayout.NORTH);
		frame.getContentPane().add(Box.createVerticalStrut(15), BorderLayout.SOUTH);
		frame.pack();

	}//END: Constructor

	private class ListenCloseWdw extends WindowAdapter {
		public void windowClosing(WindowEvent ev) {
			System.exit(0);
		}
	}

	public void launchFrame() {

		// Display Frame
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(new Dimension(dimension.width/4, dimension.height/4));
		frame.setMinimumSize(new Dimension(260, 100));
		frame.setResizable(true);
		frame.setVisible(true);
	}

	public static void main(String args[]) {

		// Start application's GUI from Event Dispatching Thread
		SwingUtilities.invokeLater(new Runnable() {

			public void run() {

				BeagleSequenceSimulatorApp gui;

				try {

					gui = new BeagleSequenceSimulatorApp();
					gui.launchFrame();

				} catch (UnsupportedClassVersionError e) {

					System.err.println("Your Java Runtime Environment is too old. Please update");
					e.printStackTrace();

				} catch (ClassNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (InstantiationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (UnsupportedLookAndFeelException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		});
	}// END: main

}// END: TestlabOutbreakApp

//dr.app.beauti.generator
//dr.app.coalgen
