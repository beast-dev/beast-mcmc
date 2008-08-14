/**
 * 
 */
package dr.evolution.tree;

import java.util.BitSet;

/**
 * @author Sebastian Hoehna
 *
 */
public class Clade implements Comparable<Clade> {

    public Clade(final BitSet bits, final double height) {
        this.bits = bits;
        this.height = height;
        size = bits.cardinality();
    }

    public BitSet getBits() {
        return bits;
    }

    public double getHeight() {
        return height;
    }

    public int getSize() {
        return size;
    }

    public int compareTo(Clade clade) {
		int setBitIndexI = 0;
		int setBitIndexJ = 0;
		do {
			setBitIndexI = bits.nextSetBit(setBitIndexI + 1);
			setBitIndexJ = clade.getBits().nextSetBit(setBitIndexJ + 1);
		} while (setBitIndexI == setBitIndexJ && setBitIndexI != -1);

		return (setBitIndexI < setBitIndexJ ? -1
				: (setBitIndexI > setBitIndexJ ? 1 : 0));
	}

    private final BitSet bits;
    private final double height;
    private final int size;
}
