/*
 * PrecisionType.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treedatalikelihood.continuous.cdi;

/**
 * @author Marc A. Suchard
 */
public enum PrecisionType {
    ELEMENTARY("elementary data vector", "elementary", 0) {
        @Override
        public void fillPrecisionInPartials(double[] partial, int offset, int index, double precision,
                                            int dimTrait) {
            partial[offset + dimTrait] = precision;
        }

        @Override
        public void fillEffDimInPartials(double[] partial, int offset, int effDim, int dimTrait) {
            //Do nothing
        }

        @Override
        public void copyObservation(double[] partial, int pOffset, double[] data, int dOffset, int dimTrait) {
            data[dOffset] = partial[pOffset];
        }

        @Override
        public int getPrecisionOffset(int dimTrait) {
            return dimTrait;
        }

        @Override
        public int getVarianceOffset(int dimTrait) {
            return -1;
        }

        @Override
        public int getEffectiveDimensionOffset(int dimTrait) {
            return -1;
        }

        @Override
        public double[] getScaledPrecision(double[] partial, int offset, double[] diffusionPrecision, int dimTrait) {
            double scalar = partial[offset + getPrecisionOffset(dimTrait)];
            return PrecisionType.scale(diffusionPrecision, scalar);
        }
    },

    SCALAR("proportional scaling per branch", "scalar", 0) {
        @Override
        public void fillPrecisionInPartials(double[] partial, int offset, int index, double precision,
                                            int dimTrait) {
            if (index == 0) {
                partial[offset + dimTrait] = precision;
            } else {
                if (partial[offset + dimTrait] != 0.0) {
                    partial[offset + dimTrait] = precision;
                }
            }
        }

        @Override
        public void fillEffDimInPartials(double[] partial, int offset, int effDim, int dimTrait) {
            //Do nothing
        }

        @Override
        public void copyObservation(double[] partial, int pOffset, double[] data, int dOffset, int dimTrait) {
            for (int i = 0; i < dimTrait; ++i) {
                data[dOffset + i] = Double.isInfinite(partial[pOffset + dimTrait]) ?
                        partial[pOffset + i] : Double.NaN;
            }
        }

        @Override
        public int getPrecisionOffset(int dimTrait) {
            return dimTrait;
        }

        @Override
        public int getVarianceOffset(int dimTrait) {
            return -1;
        }

        @Override
        public int getEffectiveDimensionOffset(int dimTrait) {
            return -1;
        }

        @Override
        public double[] getScaledPrecision(double[] partial, int offset, double[] diffusionPrecision, int dimTrait) {
            double scalar = partial[offset + getPrecisionOffset(dimTrait)];
            return PrecisionType.scale(diffusionPrecision, scalar);
        }
    },

    MIXED("mixed method", "mixed", 1) {
        @Override
        public void fillPrecisionInPartials(double[] partial, int offset, int index, double precision,
                                            int dimTrait) {
            partial[offset + dimTrait + index] = precision;
        }

        @Override
        public void fillEffDimInPartials(double[] partial, int offset, int effDim, int dimTrait) {
            //Do nothing
        }

        @Override
        public void copyObservation(double[] partial, int pOffset, double[] data, int dOffset, int dimTrait) {
            for (int i = 0; i < dimTrait; ++i) {
                data[dOffset + i] = Double.isInfinite(partial[pOffset + dimTrait + i]) ?
                        partial[pOffset + i] : Double.NaN;
            }
        }

        @Override
        public int getPrecisionOffset(int dimTrait) {
            return dimTrait;
        }

        @Override
        public int getVarianceOffset(int dimTrait) {
            return -1;
        }

        @Override
        public int getEffectiveDimensionOffset(int dimTrait) {
            return -1;
        }

        @Override
        public double[] getScaledPrecision(double[] partial, int offset, double[] diffusionPrecision, int dimTrait) {
            throw new RuntimeException("Not yet implemented");
        }
    },

    FULL("full precision matrix per branch", "full", 2) {
        // partial structure:
        //      [mean (p), precision (p^2), variance (p^2), fullPrecision (1), effective dimension (1), determinant (1), remainder (1)]

        @Override
        public void fillPrecisionInPartials(double[] partial, int offset, int index, double precision,
                                            int dimTrait) {
            final int offs = offset + dimTrait + index * dimTrait + index;
            partial[offs] = precision;
            partial[offs + dimTrait * dimTrait] = Double.isInfinite(precision) ? 0.0 : 1.0 / precision;
            partial[offset + dimTrait + 2 * dimTrait * dimTrait] = Double.POSITIVE_INFINITY;
        }

        @Override
        public void fillEffDimInPartials(double[] partial, int offset, int effDim, int dimTrait) {
            int effDimOffset = this.getEffectiveDimensionOffset(dimTrait);
            partial[offset + effDimOffset] = effDim;
        }

        @Override
        public void fillDeterminantInPartials(double[] partial, int offset, double determinant, int dimTrait) {
            int determinantOffset = getDeterminantOffset(dimTrait);
            partial[offset + determinantOffset] = determinant;
        }

        @Override
        public void fillNoDeterminantInPartials(double[] partial, int offset, int dimTrait) {
            fillDeterminantInPartials(partial, offset, Double.NaN, dimTrait); //TODO: is it bad to assume NaN is missing?
        }

        @Override
        public void fillRemainderInPartials(double[] partials, int offset, double remainder, int dimTrait) {
            int remOffset = getRemainderOffset(dimTrait);
            partials[offset + remOffset] = remainder;
        }

        @Override
        public void copyObservation(double[] partial, int pOffset, double[] data, int dOffset, int dimTrait) {
            for (int i = 0; i < dimTrait; ++i) {
                data[dOffset + i] = Double.isInfinite(partial[pOffset + dimTrait + i * dimTrait + i]) ?
                        partial[pOffset + i] : Double.NaN;
            }
        }

        @Override
        public int getPrecisionLength(int dimTrait) {
            return super.getMatrixLength(dimTrait);
        }

        @Override
        public int getVarianceLength(int dimTrait) {
            return super.getMatrixLength(dimTrait);
        }

        @Override
        public int getPrecisionOffset(int dimTrait) {
            return dimTrait;
        }

        @Override
        public int getEffectiveDimensionOffset(int dimTrait) {
            return dimTrait + dimTrait * dimTrait * 2 + 1;
        }

        @Override
        public int getDeterminantOffset(int dimTrait) {
            return dimTrait + dimTrait * dimTrait * 2 + 2;
        }

        @Override
        public int getRemainderOffset(int dimTrait) {
            return dimTrait + dimTrait * dimTrait * 2 + 3;
        }

        @Override
        public int getVarianceOffset(int dimTrait) {
            return dimTrait + dimTrait * dimTrait;
        }

        @Override
        public double[] getScaledPrecision(double[] partial, int offset, double[] diffusionPrecision, int dimTrait) {

            double[] precision = new double[dimTrait * dimTrait];
            System.arraycopy(partial, offset + getPrecisionOffset(dimTrait),
                    precision, 0,
                    dimTrait * dimTrait);

            return precision;
        }

        @Override
        public boolean hasEffectiveDimension() {
            return true;
        }

        @Override
        public boolean hasDeterminant() {
            return true;
        }

        @Override
        public boolean hasRemainder() {
            return true;
        }

        @Override
        public int getPartialsDimension(int dimTrait) {
            return dimTrait + 2 * getMatrixLength(dimTrait) + 4;
        }
    };

    private final int power;
    private final String name;
    private final String tag;

    PrecisionType(String name, String tag, int power) {
        this.name = name;
        this.tag = tag;
        this.power = power;
    }

    public String toString() {
        return name;
    }

    public int getPower() {
        return power;
    }

    protected int getMatrixLength(int dimTrait) {
        int length = 1;
        final int pow = getPower();
        for (int i = 0; i < pow; ++i) {
            length *= dimTrait;
        }
        return length;
    }

    public int getPrecisionLength(int dimTrait) {
        return getMatrixLength(dimTrait);
    }

    public int getVarianceLength(int dimTrait) {
        return 0;
    }

    public static double getObservedPrecisionValue(final boolean missing) {
        return missing ? 0.0 : Double.POSITIVE_INFINITY;
    }

    abstract public void fillPrecisionInPartials(double[] partial, int offset, int index, double precision,
                                                 int dimTrait);

    abstract public void fillEffDimInPartials(double[] partial, int offset, int effDim, int dimTrait);

    public void fillDeterminantInPartials(double[] partial, int offset, double det, int dimTrait) {
        // do nothing
    }

    public void fillNoDeterminantInPartials(double[] partial, int offset, int dimTrait) {
        // do nothing
    }

    public boolean isMissingDeterminantValue(double value) {
        return Double.isNaN(value);
    }

    public double getMissingDeterminantValue() {
        return Double.NaN;
    }

    abstract public void copyObservation(double[] partial, int pOffset, double[] data, int dOffset, int dimTrait);

    public int getMeanOffset(int dimTrait) {
        return 0;
    }

    @SuppressWarnings("unused")
    abstract public int getPrecisionOffset(int dimTrait);

    @SuppressWarnings("unused")
    abstract public int getVarianceOffset(int dimTrait);

    @SuppressWarnings("unused")
    abstract public int getEffectiveDimensionOffset(int dimTrait);

    public int getDeterminantOffset(int dimTrait) {
        return -1;
    }

    public int getRemainderOffset(int dimTrait) {
        return -1;
    }

    public void fillRemainderInPartials(double[] partials, int offset, double remainder, int dimTrait) {
        throw new RuntimeException("precision type " + tag + " does not store remainders");
    }

    abstract public double[] getScaledPrecision(double[] partial, int offset, double[] diffusionPrecision, int dimTrait);

    public int getPartialsDimension(int dimTrait) {
        return this.getMatrixLength(dimTrait) + dimTrait;
    }

    public boolean hasEffectiveDimension() {
        return false;
    }

    public boolean hasDeterminant() {
        return false;
    }

    public boolean hasRemainder() {
        return false;
    }

    private static double[] scale(double[] in, double scalar) {

        double[] out = new double[in.length];
        for (int i = 0; i < in.length; ++i) {
            out[i] = scalar * in[i];
        }
        return out;
    }

    public String getTag() {
        return tag;
    }
}
