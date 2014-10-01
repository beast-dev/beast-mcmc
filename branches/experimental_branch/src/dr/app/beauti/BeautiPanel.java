package dr.app.beauti;

import jam.framework.Exportable;

import javax.swing.*;

import dr.app.beauti.options.BeautiOptions;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public abstract class BeautiPanel extends JPanel implements Exportable {

    public abstract void setOptions(BeautiOptions options);

    public abstract void getOptions(BeautiOptions options);
}
