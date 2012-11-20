package dr.app.bss;

import jam.framework.SingleDocApplication;

import java.awt.Font;
import java.awt.Image;
import java.awt.Toolkit;
import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

public class BeagleSequenceSimulatorApp {

	public final static boolean DEBUG = true;
	
	// Share those if neccessary
	private Image beagleSequenceSimulatorImage;
	private static final String BEAGLE_SEQUENCE_SIMULATOR = "Beagle Sequence Simulator";

	// Icons
	public static ImageIcon doneIcon;
	public static ImageIcon errorIcon;
	public static ImageIcon hammerIcon;
	public static ImageIcon closeIcon;
	public static ImageIcon nuclearIcon;
	
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

			UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN,
					13));
			UIManager.put("SmallSystemFont", new Font("Lucida Grande",
					Font.PLAIN, 11));

			try {

				// UIManager.setLookAndFeel(UIManager
				// .getSystemLookAndFeelClassName());

				UIManager
						.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
				lafLoaded = true;

			} catch (Exception e) {
				//
			}

		} else {

			try {

//				UIManager.setLookAndFeel(UIManager
//						.getSystemLookAndFeelClassName());

				 UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
				lafLoaded = true;

			} catch (Exception e) {
				//
			}

		}

		if (!lafLoaded) {

			try {

				UIManager.setLookAndFeel(UIManager
						.getSystemLookAndFeelClassName());
				System.out
						.println("Specified l&f not found. Loading system default l&f");

			} catch (Exception e) {

				e.printStackTrace();

			}
		}

		SingleDocApplication app = new SingleDocApplication(
				new BeagleSequenceSimulatorMenuFactory(), //
				BEAGLE_SEQUENCE_SIMULATOR, //
				PartitionDataList.VERSION.concat(" ").concat(
						PartitionDataList.DATE_STRING), //
				null //
		);

		beagleSequenceSimulatorImage = CreateImage("icons/bss.png");
		BeagleSequenceSimulatorFrame frame = new BeagleSequenceSimulatorFrame(
				BEAGLE_SEQUENCE_SIMULATOR);
		frame.setIconImage(beagleSequenceSimulatorImage);
		app.setDocumentFrame(frame);

		// Setup icons
		doneIcon = CreateImageIcon("icons/check.png");
		errorIcon = CreateImageIcon("icons/error.png");
		hammerIcon = CreateImageIcon("icons/hammer.png");
		closeIcon = CreateImageIcon("icons/close.png");
		nuclearIcon = CreateImageIcon("icons/nuclear.png");
		
	}// END: Constructor

	public static void main(String args[]) {

		if (args.length > 1) {

			System.out.println("Command-line interface not yet implemented");

		} else {

			try {

				Thread.setDefaultUncaughtExceptionHandler(new ExceptionHandler());
				new BeagleSequenceSimulatorApp();
				// Cool sub-title for a paper
				System.out.println("Do the evolution baby!");

			} catch (UnsupportedClassVersionError e) {

				System.err
						.println("Your Java Runtime Environment is too old. Please update");
				e.printStackTrace();

			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			} catch (InstantiationException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			} catch (UnsupportedLookAndFeelException e) {
				e.printStackTrace();
			}// END: try catch block

		}// END: command line check

	}// END: main

	private Image CreateImage(String path) {
		URL imgURL = this.getClass().getResource(path);
		Toolkit kit = Toolkit.getDefaultToolkit();
		Image img = kit.createImage(imgURL);

		if (img != null) {
			return img;
		} else {
			System.out.println("Couldn't find file: " + path + "\n");
			return null;
		}

	}// END: CreateImage

	private ImageIcon CreateImageIcon(String path) {

		ImageIcon icon = null;

		try {

			URL imgURL = getClass().getResource(path);
			if (imgURL != null) {
				icon = new ImageIcon(imgURL);
			} else {
				System.err.println("Couldn't find file: " + path + "\n");
			}

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

		return icon;
	}// END: CreateImageIcon

}// END: BeagleSequenceSimulatorApp

