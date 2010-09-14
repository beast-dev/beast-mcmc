package dr.util;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface TaskListener {
    /**
     * Provides the proportion of the task complete
     * @param progress in range (0, 1)
     */
    void progress(double progress);
}
