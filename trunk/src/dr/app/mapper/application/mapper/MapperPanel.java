package dr.app.mapper.application.mapper;

import dr.app.mapper.application.MapperFrame;
import dr.inference.trace.TraceList;
import jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MapperPanel extends JPanel implements Exportable {

    public MapperPanel(MapperFrame mapperFrame) {
    }

    /**
     * This function takes a multiple statistics in a single log files
     */
    public void setTraces(TraceList[] traceLists, java.util.List<String> traces) {
    }

    public void doCopy() {

        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(getExportText());

        clipboard.setContents(selection, selection);

    }

    public String getExportText() {
        return "";
    }


   @Override
    public JComponent getExportableComponent() {
        return this;
    }
}
