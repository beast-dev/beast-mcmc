package dr.app.SnAPhyl;

import dr.app.beast.BeastVersion;
import dr.app.beauti.util.CommandLineBeauti;
import dr.app.util.OSType;
import dr.util.Version;
import org.virion.jam.framework.*;

import javax.swing.*;
import java.awt.*;


public class SnAPhylApp extends MultiDocApplication {
    private final static Version version = new BeastVersion();

    public SnAPhylApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new BeautiMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
    }

    /**
     * In a departure from the standard UI, there is no "Open" command for this application
     * Instead, the user can create a New window, Import a NEXUS file and Apply a Template file.
     * None of these operations result in a file being associated with the DocumentFrame. All
     * these actions are located in the BeautiFrame class. This overriden method should never
     * be called and throw a RuntimeException if it is.
     *
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

            new CommandLineBeauti(inputFileName, templateFileName, outputFileName);

        } else {

            if (args.length == 1 && args[0].equalsIgnoreCase("-advanced")) {
                advanced = true;
            }

            boolean lafLoaded = false;

            if (OSType.isMac()) {
                System.setProperty("apple.awt.graphics.UseQuartz", "true");
                System.setProperty("apple.awt.antialiasing","true");
                System.setProperty("apple.awt.rendering","VALUE_RENDER_QUALITY");

                System.setProperty("apple.laf.useScreenMenuBar","true");
                System.setProperty("apple.awt.draggableWindowBackground","true");
                System.setProperty("apple.awt.showGrowBox","true");

                // set the Quaqua Look and Feel in the UIManager
                try {
                    UIManager.setLookAndFeel(
                            "ch.randelshofer.quaqua.QuaquaLookAndFeel"
                    );
                    lafLoaded = true;


                } catch (Exception e) {

                }

                UIManager.put("SystemFont", new Font("Lucida Grande", Font.PLAIN, 13));
                UIManager.put("SmallSystemFont", new Font("Lucida Grande", Font.PLAIN, 11));
            }

            try {

                if (!lafLoaded) {
                    UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
                }

                java.net.URL url = SnAPhylApp.class.getResource("images/beauti.png");
                Icon icon = null;

                if (url != null) {
                    icon = new ImageIcon(url);
                }

                final String nameString = "SnAPhyl";
                final String versionString = version.getVersionString();
                String aboutString = "<html><div style=\"font-family:sans-serif;\"><center>" +
                        "<div style=\"font-size:12;\"><p>Bayesian Evolutionary Analysis Utility<br>" +
                        "Version " + versionString + ", " + version.getDateString() + "</p>" +
                        "<p>by Alexei J. Drummond, Andrew Rambaut and Walter Xie</p></div>" +
                                "<hr><div style=\"font-size:10;\">Part of the BEAST package:" +
                        version.getHTMLCredits() +
                        "</div></center></div></html>";

                String websiteURLString = "http://beast.bio.ed.ac.uk/BEAUti";
                String helpURLString = "http://beast.bio.ed.ac.uk/BEAUti";

                System.setProperty("BEAST & BEAUTi Version", version.getVersion());

                SnAPhylApp app = new SnAPhylApp(nameString, aboutString, icon,
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

    public static boolean advanced = false;
}
