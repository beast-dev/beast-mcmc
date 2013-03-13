package dr.app.mapper.application;

import dr.app.mapper.application.menus.MapperMenuBarFactory;
import dr.app.util.OSType;
import dr.util.Version;
import jam.framework.*;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.Locale;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MapperApp extends MultiDocApplication {
    private final static Version version = new Version() {
        private static final String VERSION = "1.0pre";

        public String getVersion() {
            return VERSION;
        }

        public String getVersionString() {
            return "v" + VERSION;
        }

        public String getDateString() {
            return "2013";
        }

        public String getBuildString() {
            return "Build r3656";
        }

        public String[] getCredits() {
            return new String[0];
        }

        public String getHTMLCredits() {
            return "<p>by<br>" +
                    "Andrew Rambaut, Trevor Bedford & Marc A. Suchard</p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a>, <a href=\"mailto:t.bedford@ed.ac.uk\">t.bedford@ed.ac.uk</a></p>" +

                    "<p>Departments of Biomathematics, Biostatistics and Human Genetics, UCLA<br>" +
                    "<a href=\"mailto:msuchard@ucla.edu\">msuchard@ucla.edu</a></p>" +

                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>";
        }

    };


    public MapperApp(String nameString, String aboutString, Icon icon,
                     String websiteURLString, String helpURLString) {
        super(new MapperMenuBarFactory(), nameString, aboutString, icon, websiteURLString, helpURLString);

        addPreferencesSection(new GeneralPreferencesSection());
    }

    // Main entry point
    static public void main(String[] args) {
        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

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
                //
            }

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
            java.net.URL url = MapperApp.class.getResource("images/Mapper.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            final String nameString = "Mapper";
            final String versionString = version.getVersionString();
            String aboutString = "<html><center><p>Antigenic Mapping Investigation Tool<br>" +
                    "Version " + versionString + ", " + version.getDateString() + "</p>" +
                    version.getHTMLCredits() +
                    "</center></html>";

            String websiteURLString = "http://tree.bio.ed.ac.uk/";
            String helpURLString = "http://tree.bio.ed.ac.uk/software/mapper";

            MapperApp app = new MapperApp(nameString, aboutString, icon, websiteURLString, helpURLString);
            app.setDocumentFrameFactory(new DocumentFrameFactory() {
                public DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory) {
                    return new MapperFrame("Mapper");
                }
            });
            app.initialize();

            app.doNew();

            if (args.length > 0) {
                MapperFrame frame = (MapperFrame) app.getDefaultFrame();
                for (String fileName : args) {

                    File file = new File(fileName);
//                    LogFileTraces[] traces = { new LogFileTraces(fileName, file) };

//                    frame.processTraces(traces);
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