package dr.evomodel.treedatalikelihood.continuous;

import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

public class HashedMissingArray {

    final private double[] array;

    public HashedMissingArray(final double[] array) {
        this.array = array;
    }

    public double[] getArray() {
        return array;
    }

    public double get(int index) {
        return array[index];
    }

    public int getLength() {
        return array.length;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(array);
    }

    @Override
    public boolean equals(Object obj) {
        return obj instanceof HashedMissingArray && Arrays.equals(array,
                ((HashedMissingArray) obj).array);
    }

    public String toString() {
        return new Vector(array).toString();
    }
}
