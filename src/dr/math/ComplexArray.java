package dr.math;

import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */
public class ComplexArray {

    public final double[] real;
    public final double[] complex;
    public final int length;

    public ComplexArray(double[] real) {
        this(real, new double[real.length]);
    }

    public ComplexArray(double[] real, double[] complex) {
        this.real = real;
        this.complex = complex;
        this.length = real.length;
    }

    public void conjugate() {
        for (int i = 0; i < length; ++i) {
            complex[i] = -complex[i];
        }
    }

    public ComplexArray product(ComplexArray x) {
        double[] newReal = new double[length];
        double[] newComplex = new double[length];
        for (int i = 0; i < length; ++i) {
            // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
            final double a = real[i];
            final double b = complex[i];
            final double c = x.real[i];
            final double d = x.complex[i];
            newReal[i] = a * c - b * d;
            newComplex[i] = a * d + b * c;
        }
        return new ComplexArray(newReal, newComplex);
    }

    public String toString() {
        return "\nReal   : " + new Vector(real).toString() +
               "\nComplex: " + new Vector(complex).toString();
    }

    public static double[] interleave(double[] real, double[] complex) {
        final int length = real.length;
        double[] out = new double[length * 2];
        for (int i = 0; i < length; i++) {
            out[2 * i] = real[i];
            out[2 * i + 1] = complex[i];
        }
        return out;
    }
}
