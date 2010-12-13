package dr.app.treespace;

import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TimelinePanel extends JPanel {

    private static final long serialVersionUID = -3710586474593827540L;

    private final TreeSpaceFrame frame;

    private final TreeSpaceDocument document;

    public TimelinePanel(final TreeSpaceFrame parent, final TreeSpaceDocument document) {

        this.frame = parent;
        this.document = document;

        setLayout(new BorderLayout());

        OptionsPanel optionsPanel = new OptionsPanel(12, 24);

        setOpaque(false);
        optionsPanel.setOpaque(false);

        optionsPanel.addSeparator();


        add(optionsPanel, BorderLayout.CENTER);

    }

}