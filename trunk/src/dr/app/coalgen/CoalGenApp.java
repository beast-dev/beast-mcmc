/** 
 * CoalGenApp.java
 *
 * Title:			CoalGen
 * Description:		An application for analysing MCMC tree trace files.
 * @author			Andrew Rambaut	
 * @author			Alexei Drummond	
 * @version			$Id: CoalGenApp.java,v 1.3 2004/10/01 22:40:02 alexei Exp $
 */

package dr.app.coalgen;

import org.virion.jam.framework.SingleDocApplication;

import javax.swing.*;

import dr.app.beauti.BeautiMenuBarFactory;

public class CoalGenApp  {
    public CoalGenApp() {


        try {

            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            java.net.URL url = CoalGenApp.class.getResource("/images/CoalGen.png");
            Icon icon = null;

            if (url != null) {
                icon = new ImageIcon(url);
            }

            String nameString = "CoalGen";
            String aboutString = "Coalescent Tree Simulator\nVersion 1.0\n \nCopyright 2004 Andrew Rambaut and Alexei Drummond\nAll Rights Reserved.";


            SingleDocApplication app = new SingleDocApplication(new CoalGenMenuFactory(), nameString, aboutString, icon);


            CoalGenFrame frame = new CoalGenFrame(app, nameString);
            app.setDocumentFrame(frame);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    // Main entry point
    static public void main(String[] args) {
        System.setProperty("com.apple.macos.useScreenMenuBar","true");
        System.setProperty("apple.laf.useScreenMenuBar","true");

        new CoalGenApp();
    }

}