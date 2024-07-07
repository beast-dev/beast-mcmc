package dr.app.tools.treeannotator;

import dr.evolution.util.Taxon;
import dr.util.Pair;

import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
class BiClade implements Clade {

    /**
     * Clade for a tip
     * @param index number of the tip
     */
    public BiClade(int index) {
        this.index = index;

        count = 0;
        credibility = 1.0;
        size = 1;

        key = index;
    }

    public BiClade(Clade child1, Clade child2) {
        assert child1 instanceof BiClade;
        assert child2 instanceof BiClade;

        count = 0;
        credibility = 0.0;
        BiClade left = (BiClade)child1;
        BiClade right = (BiClade)child2;
        if (right.index < left.index) {
            right = (BiClade)child1;
            left = (BiClade)child2;
        }
        size = left.size + right.size;
        index = left.index;
        addSubClades(left, right);

        key = BiClade.getKey(left, right);
    }

    public void addSubClades(Clade child1, Clade child2) {
        // arrange with the lowest index on the left
        BiClade left = (BiClade)child1;
        BiClade right = (BiClade)child2;
        if (right.index < left.index) {
            right = (BiClade)child1;
            left = (BiClade)child2;
        }
        assert left.size + right.size == size;
        assert left.index == index;
        subClades.add(new Pair<>(left, right));
    }

    @Override
    public int getCount() {
        return count;
    }

    @Override
    public void setCount(int count) {
        this.count = count;
    }

    @Override
    public double getCredibility() {
        return credibility;
    }

    @Override
    public void setCredibility(double credibility) {
        this.credibility = credibility;
    }

    @Override
    public void addAttributeValues(Object[] values) {
        if (attributeValues == null) {
            attributeValues = new ArrayList<>();
        }
        attributeValues.add(values);
    }

    @Override
    public List<Object[]> getAttributeValues() {
        return attributeValues;
    }

    @Override
    public int getSize() {
        return size;
    }

    @Override
    public Taxon getTaxon() {
        return null;
    }

    @Override
    public Clade getBestLeft() {
        return bestLeft;
    }

    @Override
    public Clade getBestRight() {
        return bestRight;
    }

    public Set<Pair<BiClade, BiClade>> getSubClades() {
        return subClades;
    }

    @Override
    public Object getKey() {
        return key;
    }

    public static Object getKey(BiClade child1, BiClade child2) {
        BitSet bits = new BitSet();
        if (child1.size == 1) {
            bits.set(child1.index);
        } else {
            bits.or((BitSet) child1.getKey());
        }
        if (child2.size == 1) {
            bits.set(child2.index);
        } else {
            bits.or((BitSet) child2.getKey());
        }
        return bits;
    }

//    public boolean equals(Object o) {
//        if (this == o) return true;
//        if (o == null || getClass() != o.getClass()) return false;
//
//        if (((BiClade) o).size != size) return false;
//
//        return !(bits != null ? !Arrays.equals(bits, ((BiClade) o).bits) : ((BiClade) o).bits != null);
//
//    }
//
//    public int hashCode() {
//        return left.hashCode() ^ right.hashCode();
//    }

    public String toString() {
        return "clade " + hashCode();
    }

    private int count;
    private double credibility;
    private final int size;
    private final int index;

    private final Object key;


    private final Set<Pair<BiClade, BiClade>> subClades = new HashSet<>();
     BiClade bestLeft = null;
     BiClade bestRight = null;


     double bestSubTreeCredibility;

    private List<Object[]> attributeValues = null;
}
