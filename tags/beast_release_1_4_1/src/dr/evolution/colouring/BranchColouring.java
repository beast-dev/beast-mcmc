package dr.evolution.colouring;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface BranchColouring {
    int getParentColour();

    int getChildColour();

    int getNumEvents();

    int getForwardColourBelow( int i );

    double getForwardTime( int i );

    int getBackwardColourAbove( int i );

    double getBackwardTime( int i );

    int getNextForwardEvent( double height );

    double getTimeInColour(int colour, double parentHeight, double childHeight);

    List getColourChanges();
}
