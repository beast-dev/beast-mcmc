package dr.app.phylogeography.spread;

import jam.framework.*;
import jam.framework.MenuBarFactory;

import javax.swing.*;
import java.awt.*;

import jam.mac.Utils;

/**
 * @author Andrew Rambaut
 * @version $Id: BeautiApp.java,v 1.18 2006/09/09 16:07:05 rambaut Exp $
 */
public class SpreadApp extends MultiDocApplication {

    public SpreadApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new SpreadMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);
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

//            if (args.length != 3) {
//                System.err.println("Usage: phylogeography <input_file> <template_file> <output_file>");
//                return;
//            }
//
//            String inputFileName = args[0];
//            String templateFileName = args[1];
//            String outputFileName = args[2];
//
//            new CommandLineBeauti(inputFileName, templateFileName, outputFileName);

        } else {

            boolean lafLoaded = false;

            if (Utils.isMacOSX()) {
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

                java.net.URL url = SpreadApp.class.getResource("images/phylogeography.png");
                Icon icon = null;

                if (url != null) {
                    icon = new ImageIcon(url);
                }

                final String nameString = "S.P.R.E.A.D.";
                final String versionString = "v1.0";
                String aboutString = "<html><div style=\"font-family:sans-serif;\"><center>" +
                        "<div style=\"font-size:12;\"><p>Spatial Phylogenetic Reconstruction<br>" +
                        "Version " + versionString + ", 2009</p>" +
                        "<p>by Philippe Lemey, Marc Suchard & Andrew Rambaut</p></div>" +
                        "</div></html>";

                String websiteURLString = "http://beast.bio.ed.ac.uk/";
                String helpURLString = "http://beast.bio.ed.ac.uk/phylogeography/";

                SpreadApp app = new SpreadApp(nameString, aboutString, icon,
                        websiteURLString, helpURLString);
                app.setDocumentFrameFactory(new DocumentFrameFactory() {
                    public DocumentFrame createDocumentFrame(final Application application, final MenuBarFactory menuBarFactory) {
                        return new SpreadFrame(nameString);
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
