/*
 * CladeKey.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.app.tools.treeannotator;

import jebl.math.Random;

import java.util.Arrays;
import java.util.BitSet;

public class CladeKey {
    /*
     * BitSets are packed into arrays of "words."  Currently a word is
     * a long, which consists of 64 bits, requiring 6 address bits.
     * The choice of word size is determined purely by performance concerns.
     */
    private static final int ADDRESS_BITS_PER_WORD = 6;

    private final int maxIndex;

    /**
     * The internal field corresponding to the serialField "bits".
     */
    private long[] words;

    /**
     * The number of words in the logical size of this BitSet.
     */
    private transient int wordsInUse = 0;
    private int hashCache = 0;

    /**
     * Given a bit index, return word index containing it.
     */
    private static int wordIndex(int bitIndex) {
        return bitIndex >> ADDRESS_BITS_PER_WORD;
    }

    /**
     * Creates a new bit set to accommodate bit indices up to maxIndex. All bits are initially {@code false}.
     */
    public CladeKey(int maxIndex) {
        this.maxIndex = maxIndex;
        initWords(maxIndex);
    }

    /**
     * Sets the field wordsInUse to the logical size in words of the bit set.
     * WARNING:This method assumes that the number of words actually in use is
     * less than or equal to the current value of wordsInUse!
     */
    private void recalculateWordsInUse() {
        // Traverse the bitset until a used word is found
        int i;
        for (i = wordsInUse-1; i >= 0; i--) {
            if (words[i] != 0)
                break;
        }
        wordsInUse = i+1; // The new logical size
    }


    private void initWords(int nbits) {
        words = new long[wordIndex(nbits) + 1];
    }

    public void set(int bitIndex) {
        assert bitIndex >= 0 && bitIndex <= maxIndex;

        int wordIndex = wordIndex(bitIndex);

        wordsInUse = Math.max(wordsInUse, wordIndex + 1);

        words[wordIndex] |= (1L << bitIndex);

        hashCache = 0;
    }

    public void setTo(CladeKey key) {
        assert this != key;

        wordsInUse = key.wordsInUse;
        System.arraycopy(key.words, 0,
                words, 0,
                wordsInUse);
        hashCache = 0;
    }

    public void clear() {
        while (wordsInUse > 0) {
            words[--wordsInUse] = 0;
        }
        hashCache = 0;
    }

    /**
     * Returns the number of bits set to {@code true} in this {@code BitSet}.
     *
     * @return the number of bits set to {@code true} in this {@code BitSet}
     * @since  1.4
     */
    public int cardinality() {
        int sum = 0;
        for (int i = 0; i < wordsInUse; i++) {
            sum += Long.bitCount(words[i]);
        }
        return sum;
    }

    public void and(CladeKey key) {
        assert this != key;

        while (wordsInUse > key.wordsInUse) {
            words[--wordsInUse] = 0;
        }

        // Perform logical AND on words in common
        for (int i = 0; i < wordsInUse; i++) {
            words[i] &= key.words[i];
        }

        recalculateWordsInUse();
        hashCache = 0;
    }

    public void and(CladeKey key1, CladeKey key2) {
        assert this != key1 && this != key2 && key1 != key2;

        if (key1.wordsInUse >= key2.wordsInUse) {
            wordsInUse = key1.wordsInUse;
            words = Arrays.copyOf(key1.words, key1.wordsInUse);
            for (int i = 0; i < key2.wordsInUse; i++) {
                words[i] &= key2.words[i];
            }
            for (int i = key2.wordsInUse; i < wordsInUse; i++) {
                words[i] = 0;
            }
        } else {
            wordsInUse = key2.wordsInUse;
            words = Arrays.copyOf(key2.words, key2.wordsInUse);
            for (int i = 0; i < key1.wordsInUse; i++) {
                words[i] &= key1.words[i];
            }
            for (int i = key1.wordsInUse; i < wordsInUse; i++) {
                words[i] = 0;
            }
        }
        hashCache = 0;
    }


    public void or(CladeKey key) {
        assert this != key;
        assert getMaxIndex() >= key.getMaxIndex();
        
        int wordsInCommon = Math.min(wordsInUse, key.wordsInUse);
        wordsInUse = Math.max(wordsInUse, key.wordsInUse);

        // Perform logical OR on words in common
        for (int i = 0; i < wordsInCommon; i++) {
            words[i] |= key.words[i];
        }

        // Copy any remaining words
        if (wordsInCommon < key.wordsInUse) {
            System.arraycopy(key.words, wordsInCommon,
                    words, wordsInCommon,
                    wordsInUse - wordsInCommon);
        }
        hashCache = 0;
    }

    public void or(CladeKey key1, CladeKey key2) {
        assert this != key1 && this != key2 && key1 != key2;

        if (key1.wordsInUse >= key2.wordsInUse) {
            wordsInUse = key1.wordsInUse;
            System.arraycopy(key1.words, 0,
                    words, 0,
                    wordsInUse);

//            words = Arrays.copyOf(key1.words, key1.wordsInUse);
            for (int i = 0; i < key2.wordsInUse; i++) {
                words[i] |= key2.words[i];
            }
        } else {
            wordsInUse = key2.wordsInUse;
            System.arraycopy(key2.words, 0,
                    words, 0,
                    wordsInUse);
//            words = Arrays.copyOf(key2.words, key2.wordsInUse);
            for (int i = 0; i < key1.wordsInUse; i++) {
                words[i] |= key1.words[i];
            }
        }
        hashCache = 0;
    }


    /**
     * Returns the hash code value for this bit set. The hash code depends
     * only on which bits are set within this {@code BitSet}.
     *
     * <p>The hash code is defined to be the result of the following
     * calculation:
     *  <pre> {@code
     * public int hashCode() {
     *     long h = 1234;
     *     long[] words = toLongArray();
     *     for (int i = words.length; --i >= 0; )
     *         h ^= words[i] * (i + 1);
     *     return (int)((h >> 32) ^ h);
     * }}</pre>
     * Note that the hash code changes if the set of bits is altered.
     *
     * @return the hash code value for this bit set
     */
    public int hashCode() {
        if (hashCache == 0) {
            hashCache = computeHashCode();
        }
        return hashCache;
    }

    private int computeHashCode() {
        long h = 1234;
        for (int i = wordsInUse; --i >= 0; ) {
            h ^= words[i] * (i + 1);
        }
        return (int)((h >> 32) ^ h);
    }

    public boolean equals(Object obj) {
        CladeKey key = (CladeKey)obj;

        if (wordsInUse != key.wordsInUse)
            return false;

        // Check words in use by both BitSets
        for (int i = 0; i < wordsInUse; i++)
            if (words[i] != key.words[i])
                return false;

        return true;
    }

    public int getMaxIndex() {
        return maxIndex;
    }

    public static void main(String[] argv) {
        int size = 10000;
        int count = 100000000;

        BitSet bitset1 = new BitSet();
        BitSet bitset2 = new BitSet();
        for (int i = 0; i < 1000; i++) {
            bitset1.set(Random.nextInt(size));
            bitset2.set(Random.nextInt(size));
        }


        BitSet bitset = new BitSet();

        long startTime = System.currentTimeMillis();
        for (int i = 0; i < count; i++) {
            bitset.clear();
            bitset.or(bitset1);
            bitset.or(bitset2);
        }
        long time = System.currentTimeMillis() - startTime;
        System.out.println("Bitset clear/or/or: " + time);

        CladeKey cladeKey1 = new CladeKey(size);
        CladeKey cladeKey2 = new CladeKey(size);
        for (int i = 0; i < 1000; i++) {
            cladeKey1.set(Random.nextInt(size));
            cladeKey2.set(Random.nextInt(size));
        }

        CladeKey cladeKey = new CladeKey(size);

        startTime = System.currentTimeMillis();
        for (int i = 0; i < count; i++) {
//            cladeKey.setTo(cladeKey1);
            cladeKey.clear();
            cladeKey.or(cladeKey1);
            cladeKey.or(cladeKey2);
        }
         time = System.currentTimeMillis() - startTime;
        System.out.println("CladeKey clear/or/or: " + time);

        startTime = System.currentTimeMillis();
        for (int i = 0; i < count; i++) {
            cladeKey.setTo(cladeKey1);
            cladeKey.or(cladeKey2);
        }
        time = System.currentTimeMillis() - startTime;
        System.out.println("CladeKey setTo/or: " + time);

        // count = 10^9
        // size = 1024
//        Bitset clear/or/or: 12977
//        CladeKey clear/or/or: 12798
//        CladeKey setTo/or: 8594

        // size = 1600
//        Bitset clear/or/or: 15151
//        CladeKey clear/or/or: 15373
//        CladeKey setTo/or: 10699

        // size = 10000
//        Bitset clear/or/or: 59138
//        CladeKey clear/or/or: 59450
//        CladeKey setTo/or: 36380


    }

    private static final int BITS_PER_WORD = 1 << ADDRESS_BITS_PER_WORD;
    private static final long WORD_MASK = 0xffffffffffffffffL;

    public int nextSetBit(int fromIndex) {
        if (fromIndex < 0)
            throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);

        int u = wordIndex(fromIndex);
        if (u >= wordsInUse)
            return -1;

        long word = words[u] & (WORD_MASK << fromIndex);

        while (true) {
            if (word != 0)
                return (u * BITS_PER_WORD) + Long.numberOfTrailingZeros(word);
            if (++u == wordsInUse)
                return -1;
            word = words[u];
        }
    }

    public int nextClearBit(int fromIndex) {
        // Neither spec nor implementation handle bitsets of maximal length.
        // See 4816253.
        if (fromIndex < 0)
            throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);

        int u = wordIndex(fromIndex);
        if (u >= wordsInUse)
            return fromIndex;

        long word = ~words[u] & (WORD_MASK << fromIndex);

        while (true) {
            if (word != 0)
                return (u * BITS_PER_WORD) + Long.numberOfTrailingZeros(word);
            if (++u == wordsInUse)
                return wordsInUse * BITS_PER_WORD;
            word = ~words[u];
        }
    }

    public String toString() {
        final int MAX_INITIAL_CAPACITY = Integer.MAX_VALUE - 8;
        int numBits = (wordsInUse > 128) ?
                cardinality() : wordsInUse * BITS_PER_WORD;
        // Avoid overflow in the case of a humongous numBits
        int initialCapacity = (numBits <= (MAX_INITIAL_CAPACITY - 2) / 6) ?
                6 * numBits + 2 : MAX_INITIAL_CAPACITY;
        StringBuilder b = new StringBuilder(initialCapacity);
        b.append('{');

        int i = nextSetBit(0);
        if (i != -1) {
            b.append(i);
            while (true) {
                if (++i < 0) break;
                if ((i = nextSetBit(i)) < 0) break;
                int endOfRun = nextClearBit(i);
                do { b.append(", ").append(i); }
                while (++i != endOfRun);
            }
        }

        b.append('}');
        return b.toString();
    }
}
