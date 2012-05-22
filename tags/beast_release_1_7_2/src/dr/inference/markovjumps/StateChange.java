package dr.inference.markovjumps;

/**
 * A class to represent a state-change in a continuous-time Markov chain.
 * This class is equivalent to ColourChange in terms of parameters, but labels events
 * more generically as either forward or backwards (instead of just backwards) in time.
 * <p/>
 * TODO Remove dr.evolution.tree.ColourChange (wrong package)
 *
 * @author Alexei Drummond
 * @author Marc A. Suchard
 *         <p/>
 *         This work is supported by NSF grant 0856099
 */

public class StateChange {

    // time of the change in colour
    private double time;

    // the state associated with this time
    private int state;

    public StateChange(StateChange change) {
        this(change.time, change.state);
    }

    public StateChange(double time, int state) {
        this.time = time;
        this.state = state;
    }

    /**
     * @return the time of this change in colour
     */
    public final double getTime() {
        return time;
    }

    /**
     * @return the state associated with the time
     */
    public final int getState() {
        return state;
    }

    /**
     * @param state sets the state of this event
     */
    public final void setState(int state) {
        this.state = state;
    }

    public final void setTime(double time) {
        this.time = time;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("{").append(getTime()).append(",").append(getState()).append("}");
        return sb.toString();
    }
}
