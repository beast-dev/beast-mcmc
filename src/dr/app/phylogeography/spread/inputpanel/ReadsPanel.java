package pyromania.sourcespanel;

import pyromania.annotationspanel.AnnotationsPanel;
import pyromania.database.Pyromania;
import pyromania.PyromaniaFrame;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class ReadsPanel extends JPanel {

    private PyromaniaFrame frame = null;

    private final Pyromania document;

    public ReadsPanel(final PyromaniaFrame parent, final Pyromania document) {
        super(new BorderLayout());

        setBackground(Color.white);

        this.frame = parent;
        this.document = document;

        AnnotationsPanel panel1 = new AnnotationsPanel(parent, document);
        AnnotationsPanel panel2 = new AnnotationsPanel(parent, document);

        JPanel panel3 = new JPanel();
        panel3.setSize(800,400);
        panel3.setBackground(Color.white);

        JSplitPane splitPane1 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, panel2);
        splitPane1.setContinuousLayout(true);
        splitPane1.setResizeWeight(0.5);

        splitPane1.setDividerSize(5);
        splitPane1.putClientProperty("Quaqua.SplitPane.style", "bar");

        JSplitPane splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, splitPane1, panel3);
        splitPane2.setContinuousLayout(true);

        splitPane2.setDividerSize(5);
        splitPane2.putClientProperty("Quaqua.SplitPane.style", "bar");

        add(splitPane2, BorderLayout.CENTER);
    }



}