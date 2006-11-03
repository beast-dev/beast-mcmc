/** 
 * CoalGenMenuFactory.java
 *
 * Title:			CoalGen
 * Description:		An application for analysing MCMC tree trace files.
 * @author			Andrew Rambaut	
 * @author			Alexei Drummond	
 * @version			$Id: CoalGenMenuFactory.java,v 1.3 2004/10/01 22:40:02 alexei Exp $
 */

package dr.app.coalgen;

import org.virion.jam.framework.Application;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.framework.SingleDocMenuBarFactory;

import javax.swing.*;

public class CoalGenMenuFactory extends SingleDocMenuBarFactory {
	
	public void populateMenuBar(JMenuBar menuBar, DocumentFrame documentFrame, 
								Application application) {
								
		super.populateMenuBar(menuBar, documentFrame);
		
		documentFrame.getSaveAction().setEnabled(false);
		documentFrame.getSaveAsAction().setEnabled(false);

		documentFrame.getCutAction().setEnabled(false);
		documentFrame.getCopyAction().setEnabled(true);
		documentFrame.getPasteAction().setEnabled(false);
		documentFrame.getDeleteAction().setEnabled(false);
		documentFrame.getSelectAllAction().setEnabled(false);
		documentFrame.getFindAction().setEnabled(false);

		documentFrame.getZoomWindowAction().setEnabled(false);
	}

}