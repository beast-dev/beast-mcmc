package dr.app.phylogeography.spread;

import dr.app.phylogeography.generator.Generator;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.List;
import java.io.File;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class OutputPanel extends JPanel {

    private static final long serialVersionUID = -3710586474593827540L;

    private final SpreadFrame frame;

    private final JComboBox formatCombo;

    private final SpreadDocument document;

    public OutputPanel(final SpreadFrame parent, final SpreadDocument document, final List<Generator> generators) {
        this.frame = parent;
        this.document = document;

        setLayout(new BorderLayout());

        OptionsPanel optionsPanel = new OptionsPanel(12, 24);

        setOpaque(false);
        optionsPanel.setOpaque(false);

        formatCombo = new JComboBox(generators.toArray());
        optionsPanel.addComponentWithLabel("Output format:", formatCombo);

        optionsPanel.addSeparator();

        add(optionsPanel, BorderLayout.CENTER);

    }

}