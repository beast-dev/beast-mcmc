package dr.evomodel.bigfasttree;

import java.util.BitSet;

/**
 * An interface for a clade reference for use in a clade model
 *
 * @author  JT McCrone
 */
public interface CladeRef {
    int getNumber();
    BitSet getBits();
}
