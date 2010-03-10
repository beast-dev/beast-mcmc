package dr.app.phylogeography.spread;

import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TimelinePanel extends JPanel {

    private static final long serialVersionUID = -3710586474593827540L;

    private final SpreadFrame frame;

    private final SpreadDocument document;

    public TimelinePanel(final SpreadFrame parent, final SpreadDocument document) {

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