package dr.app.tracer.application;

import dr.inference.trace.LogFileTraces;
import org.virion.jam.framework.*;
import org.virion.jam.mac.Utils;

import javax.swing.*;
import java.awt.*;
import java.io.File;

public class TracerApp extends MultiDocApplication {

	public TracerApp(String nameString, String aboutString, Icon icon,
	                 String websiteURLString, String helpURLString) {
		super(new TracerMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
	}

	// Main entry point
	static public void main(String[] args) {
		boolean lafLoaded = false;

        if (Utils.isMacOSX()) {
            System.setProperty("apple.awt.graphics.UseQuartz", "true");
            System.setProperty("apple.awt.antialiasing","true");
            System.setProperty("apple.awt.rendering","VALUE_RENDER_QUALITY");

            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.draggableWindowBackground","true");
            System.setProperty("apple.awt.showGrowBox","true");

            // set the Quaqua Look and Feel in the UIManager
//            try {
//                LookAndFeel lafClass;
//
//                if (!Utils.getMacOSXVersion().startsWith("10.5")) {
////                    lafClass = ch.randelshofer.quaqua.subset.Quaqua14ColorChooserLAF.class.newInstance();
////                } else {
//                    lafClass = QuaquaLookAndFeel.class.newInstance();
//                    UIManager.setLookAndFeel(lafClass);
//                    lafLoaded = true;
//                }
//
//
//            } catch (Exception e) {
//            }

            UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
            UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
        }

        if (!lafLoaded) {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception e1) {
                e1.printStackTrace();
            }
        }

		try {
			java.net.URL url = TracerApp.class.getResource("images/Tracer.png");
			Icon icon = null;

			if (url != null) {
				icon = new ImageIcon(url);
			}

			final String nameString = "Tracer";
			final String versionString = "v1.4.1";
			String aboutString = "<html><center><p>MCMC Trace Analysis Tool<br>" +
					"Version " + versionString + ", 2003-2008</p>" +
					"<p>by<br>" +
					"Andrew Rambaut and Alexei J. Drummond</p>" +
					"<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
					"<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
					"<p>Department of Computer Science, University of Auckland<br>" +
					"<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
					"<p>Available from the BEAST site:<br>" +
					"<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
					"<p>Source code distributed under the GNU LGPL:<br>" +
					"<a href=\"http://code.google.com/p/beast-mcmc/\">http://code.google.com/p/beast-mcmc/</a></p>" +
					"<p>Thanks for contributions to: Joseph Heled, Oliver Pybus, Benjamin Redelings & Marc Suchard</p>" +
					"</center></html>";

			String websiteURLString = "http://beast.bio.ed.ac.uk/";
			String helpURLString = "http://beast.bio.ed.ac.uk/Tracer";

			TracerApp app = new TracerApp(nameString, aboutString, icon, websiteURLString, helpURLString);
			app.setDocumentFrameFactory(new DocumentFrameFactory() {
				public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
					return new TracerFrame("Tracer");
				}
			});
			app.initialize();

			app.doNew();

			if (args.length > 0) {
				TracerFrame frame = (TracerFrame) app.getDefaultFrame();
				for (String fileName : args) {

					File file = new File(fileName);
					LogFileTraces[] traces = new LogFileTraces[] { new LogFileTraces(fileName, file) };

					frame.processTraces(traces);
				}
			}
		} catch (Exception e) {
			JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
					"Please report this to the authors",
					JOptionPane.ERROR_MESSAGE);
			e.printStackTrace();
		}
	}

}