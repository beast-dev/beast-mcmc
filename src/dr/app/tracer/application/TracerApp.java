package dr.app.tracer.application;

import org.virion.jam.framework.*;

import javax.swing.*;
import java.awt.*;

public class TracerApp {
	public TracerApp() {
		java.net.URL url = TracerApp.class.getResource("images/Tracer.png");
		Icon icon = null;

		if (url != null) {
			icon = new ImageIcon(url);
		}

		String nameString = "Tracer";
		String aboutString = "MCMC Trace Analysis Tool\nVersion 1.4\n \nCopyright 2003-2006 Andrew Rambaut & Alexei Drummond\nUniversity of Oxford\nAll Rights Reserved.";

		MultiDocApplication app = new MultiDocApplication(nameString, aboutString, icon);

		Application.getMenuBarFactory().registerMenuFactory(new AnalysisMenuFactory());

		app.setDocumentFrameFactory(new DocumentFrameFactory() {
			public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
				return new TracerFrame("Tracer");
			}
		});


		app.initialize();

        app.doNew();
	}

	// Main entry point
	static public void main(String[] args) {
		System.setProperty("com.apple.macos.useScreenMenuBar","true");
		System.setProperty("apple.laf.useScreenMenuBar","true");
		System.setProperty("apple.awt.showGrowBox","true");
		System.setProperty("apple.awt.antialiasing","on");
		System.setProperty("apple.awt.textantialiasing","on");
		System.setProperty("apple.awt.rendering","VALUE_RENDER_SPEED");

		// set the Quaqua Look and Feel in the UIManager
		try {
            //System.setProperty("Quaqua.Debug.showClipBounds","true");
            //System.setProperty("Quaqua.Debug.showVisualBounds","true");
			UIManager.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
			// set UI manager properties here that affect Quaqua
			UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
			UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));

        } catch (Exception e) {
			try {

				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			} catch (Exception e1) {
				e1.printStackTrace();
			}
		}

		new TracerApp();
	}

}