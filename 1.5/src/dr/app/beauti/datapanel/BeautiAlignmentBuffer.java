package dr.app.beauti.datapanel;


import java.util.List;

import dr.app.beauti.alignmentviewer.AlignmentBuffer;
import dr.app.beauti.alignmentviewer.AlignmentBufferListener;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;

/**
 * @author Andrew Rambaut
 * @version $Id: SimpleAlignmentBuffer.java,v 1.2 2005/12/11 22:41:25 rambaut Exp $
 */
public class BeautiAlignmentBuffer implements AlignmentBuffer {

    public BeautiAlignmentBuffer(Alignment alignment) {
        this.alignment = alignment;
        DataType type = alignment.getDataType();

        stateTable = new String[type.getAmbiguousStateCount()];
        for (int i = 0; i < stateTable.length; i++) {
            stateTable[i] = Character.toString(type.getChar(i));
        }

        gapState = (byte)type.getGapState();
    }

    public int getSequenceCount() {
        return alignment.getSequenceCount();
    }

    public int getSiteCount() {
        return alignment.getSiteCount();
    }

    public String getTaxonLabel(int i) {
        return alignment.getTaxonId(i);
    }

    public String[] getStateTable() {
        return stateTable;
    }

    public void getStates(int sequenceIndex, int fromSite, int toSite, byte[] states) {
        Sequence sequence = alignment.getSequence(sequenceIndex);
        int j = 0;
        for (int i = fromSite; i <= toSite; i++) {
            states[j] = (byte)sequence.getState(i);
            j++;
        }
    }

    public void addAlignmentBufferListener(AlignmentBufferListener listener) {
    }

    private final Alignment alignment;
    private final String[] stateTable;
    private final byte gapState;

}
