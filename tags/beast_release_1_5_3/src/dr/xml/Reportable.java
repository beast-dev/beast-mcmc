package dr.xml;

/**
 * This is an alternative interface for items added to reports. By default 'toString' will be called
 * but for some objects this causes issues with debugging because the debugger also calls toString which
 * may cause a whole lot of computation which messes up the debugging. If a class implements this interface
 * then this will be called instead.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface Reportable {
    String getReport();
}
