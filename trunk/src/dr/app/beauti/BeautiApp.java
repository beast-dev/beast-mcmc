/*
 * BeautiApp.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package dr.app.beauti;

import dr.app.beast.BeastVersion;
import dr.util.Version;
import org.virion.jam.framework.*;

import javax.swing.*;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: BeautiApp.java,v 1.18 2006/09/09 16:07:05 rambaut Exp $
 */
public class BeautiApp extends MultiDocApplication {
	private final static Version version = new BeastVersion();

    public BeautiApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new BeautiMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    /**
     * In a departure from the standard UI, there is no "Open" command for this application
     * Instead, the user can create a New window, Import a NEXUS file and Apply a Template file.
     * None of these operations result in a file being associated with the DocumentFrame. All
     * these actions are located in the BeautiFrame class. This overriden method should never
     * be called and throw a RuntimeException if it is.
     * @return the action
     */
    public Action getOpenAction() {
        throw new UnsupportedOperationException("getOpenAction is not supported");
    }

    // Main entry point
    static public void main(String[] args) {


        if (args.length > 1) {

            if (args.length != 3) {
                System.err.println("Usage: beauti <input_file> <template_file> <output_file>");
                return;
            }

            String inputFileName = args[0];
            String templateFileName = args[1];
            String outputFileName = args[2];

            CommandLineBeauti beauti = new CommandLineBeauti(inputFileName, templateFileName, outputFileName);

        } else {

	        if (args.length == 1 && args[0].equalsIgnoreCase("-developer")) {
		        developer = true;
	        }

            System.setProperty("com.apple.macos.useScreenMenuBar","true");
            System.setProperty("apple.laf.useScreenMenuBar","true");

            try {

                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

                java.net.URL url = BeautiApp.class.getResource("images/beauti.png");
                Icon icon = null;

                if (url != null) {
                    icon = new ImageIcon(url);
                }

                final String nameString = "BEAUti";
                final String versionString = version.getVersionString();
                String aboutString = "<html><center><p>Bayesian Evolutionary Analysis Utility<br>" +
                        "Version " + versionString + ", " + version.getDateString() + "</p>" +
                        "<p>by<br>" +
                        "Andrew Rambaut and Alexei J. Drummond</p>" +
                        "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                        "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                        "<p>Department of Computer Science, University of Auckland<br>" +
                        "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                        "<p>Part of the BEAST package:<br>" +
                        "<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
                        "</center></html>";

                String websiteURLString = "http://beast.bio.ed.ac.uk/";
                String helpURLString = "http://beast.bio.ed.ac.uk/BEAUti/";

                BeautiApp app = new BeautiApp(nameString, aboutString, icon,
                        websiteURLString, helpURLString);
                app.setDocumentFrameFactory(new DocumentFrameFactory() {
                    public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                        return new BeautiFrame(nameString);
                    }
                });
                app.initialize();
                app.doNew();
            } catch (Exception e) {
                JOptionPane.showMessageDialog(new JFrame(), "Fatal exception: " + e,
                        "Please report this to the authors",
                        JOptionPane.ERROR_MESSAGE);
                e.printStackTrace();
            }
        }
    }

	public static boolean developer = false;
}
