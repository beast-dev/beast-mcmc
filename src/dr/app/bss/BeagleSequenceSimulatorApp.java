package dr.app.bss;

import jam.framework.SingleDocApplication;

import java.awt.Font;
import java.awt.Image;
import java.awt.Toolkit;
import java.net.URL;
import java.util.Locale;

import javax.swing.ImageIcon;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class BeagleSequenceSimulatorApp {

	public static final boolean DEBUG = false;
	public static final boolean VERBOSE = true;
	
	// Share those if neccessary
	public static final String SHORT_NAME = "\u03C0BUSS";
	public static final String LONG_NAME = "Parallel BEAST/BEAGLE Utility for Sequence Simulation";
	public static final String VERSION = "1.3.7";
	public static final String DATE = "2014";

	// Icons
	private Image beagleSequenceSimulatorImage;
	
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

				// UIManager.setLookAndFeel(UIManager
				// .getSystemLookAndFeelClassName());

				UIManager
						.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
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

		ImageIcon bssIcon = Utils.createImageIcon(Utils.BSS_ICON);
		SingleDocApplication app = new SingleDocApplication(new MenuBarFactory(), //
				LONG_NAME, //
				VERSION.concat(" ").concat(
				DATE), //
				bssIcon //
		);

		beagleSequenceSimulatorImage = CreateImage(Utils.BSS_ICON);
		MainFrame frame = new MainFrame(SHORT_NAME);
		frame.setIconImage(beagleSequenceSimulatorImage);
		app.setDocumentFrame(frame);

	}// END: Constructor

	public static void main(String args[]) {

		 Locale.setDefault(Locale.US);
		
		if (args.length > 0) {
			
				BeagleSequenceSimulatorConsoleApp app = new BeagleSequenceSimulatorConsoleApp();
				app.simulate(args);
				
		} else {

			try {

				Thread.setDefaultUncaughtExceptionHandler(new ExceptionHandler());
				new BeagleSequenceSimulatorApp();

				// sub-title for a software note
				System.out.println("Do the evolution baby!");

			} catch (UnsupportedClassVersionError e) {

				System.out
						.println("Your Java Runtime Environment is too old. Please update");

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

}// END: class

