package dr.evomodel.coalescent;

/**
 * An interface to be used with DemographicLogger. Should provide
 * the time intervals between change points and the popsizes that
 * define the intervals. For 'stepwise' their is one popsize per
 * interval. For 'linear' and 'exponential' change there is one
 * fewer interval than popsizes.
 * 
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface DemographicReconstructor {

    enum ChangeType {
        STEPWISE,
        LINEAR,
        EXPONENTIAL
    };

    String getTitle();

    ChangeType getChangeType();

    double[] getIntervals();
    double[] getPopSizes();
    
}
