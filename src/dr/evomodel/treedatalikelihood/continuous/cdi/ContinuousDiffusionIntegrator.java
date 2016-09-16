/*
 * ContinuousDiffusionIntegrator.java
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
public interface ContinuousDiffusionIntegrator {
    int OPERATION_TUPLE_SIZE = 5;
    int NONE = -1;

    void finalize() throws Throwable;

//    void setPartialMean(int bufferIndex, final double[] mean);

//    void getPartialMean(int bufferIndex, final double[] mean);

//    void setPartialPrecision(int bufferIndex, final double[] precision);

//    void getPartialPrecision(int bufferIndex, final double[] precision);

    void setPartial(int bufferIndex, final double[] partial);

    void getPartial(int bufferIndex, final double[] partial);

    void setDiffusionPrecision(int diffusionIndex, final double[] matrix);

    void updatePartials(final int[] operations, int operationCount);

    InstanceDetails getDetails();

    void updateDiffusionMatrices(int offsetIndex, final int[] probabilityIndices, final double[] edgeLengths,
                                 int updateCount);

//    abstract class AbstractBase implements ContinuousDiffusionIntegrator {
//
//        private InstanceDetails details = new InstanceDetails();
//
//        @Override
//        public void finalize() throws Throwable {
//            super.finalize();
//        }
//
//        @Override
//        public void setPartialMean(int bufferIndex, double[] mean) {
//
//        }
//
//        @Override
//        public void getPartialMean(int bufferIndex, double[] mean) {
//
//        }
//
//        @Override
//        public void setPartialPrecision(int bufferIndex, double[] precision) {
//
//        }
//
//        @Override
//        public void getPartialPrecision(int bufferIndex, double[] precision) {
//
//        }
//
//
//
//        @Override
//        public abstract void setPartial(int bufferIndex, double[] partial);
//
//        @Override
//        public abstract void getPartial(int bufferIndex, double[] partial);
//
//        @Override
//        public abstract void setDiffusionPrecision(int diffusionIndex, double[] matrix);
//
//        @Override
//        public abstract void updatePartials(int[] operations, int operationCount);
//
//    }

    class Basic implements ContinuousDiffusionIntegrator {

        private int instance = -1;
        private InstanceDetails details = new InstanceDetails();

        private final PrecisionType precisionType;
        private final int numTraits;
        private final int dimTrait;
        private final int bufferCount;
        private final int diffusionCount;

        private final int dimPartial;

        public Basic(
                final PrecisionType precisionType,
                final int numTraits,
                final int dimTrait,
                final int bufferCount,
                final int diffusionCount
        ) {
            assert(numTraits > 0);
            assert(dimTrait > 0);
            assert(bufferCount > 0);
            assert(diffusionCount > 0);

            this.precisionType = precisionType;
            this.numTraits = numTraits;
            this.dimTrait = dimTrait;
            this.bufferCount = bufferCount;
            this.diffusionCount = diffusionCount;
            this.dimPartial = precisionType.getPartialLength(numTraits, dimTrait);

            allocateStorage();
        }

        @Override
        public void finalize() throws Throwable {
            super.finalize();
        }

        @Override
        public void setPartial(int bufferIndex, double[] partial) {
            assert(partial.length == dimPartial);
            assert(partials != null);

            System.arraycopy(partial, 0, partials, dimPartial * bufferIndex, dimPartial);
        }

        @Override
        public void getPartial(int bufferIndex, double[] partial) {
            assert(partial.length == dimPartial);
            assert(partials != null);

            System.arraycopy(partials, dimPartial * bufferIndex, partial, 0, dimPartial);
        }

        @Override
        public void setDiffusionPrecision(int precisionIndex, final double[] matrix) {

        }

        @Override
        public void updatePartials(final int[] operations, int operationCount) {

            System.err.println("Operations:");
            int offset = 0;
            for (int i = 0; i < operationCount; ++i) {
                System.err.println("\t" + getOperationString(operations, offset));

                updatePartial(
                        operations[offset + 0],
                        operations[offset + 1],
                        operations[offset + 2],
                        operations[offset + 3],
                        operations[offset + 4]
                );

                offset += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
            }
            System.err.println("End");
            System.err.println("");
//            System.exit(-1);
        }

        @Override
        public void updateDiffusionMatrices(int offsetIndex, final int[] probabilityIndices,
                                            final double[] edgeLengths, int updateCount) {

            System.err.println("Matrices:");
            for (int i = 0; i < updateCount; ++i) {
                System.err.println("\t" + probabilityIndices[i] + " <- " + edgeLengths[i]);

                // TODO Currently only writtern for SCALAR model
                final int dimMatrix = 1;
                diffusions[dimMatrix * probabilityIndices[i]] = edgeLengths[i];
            }
        }

        @Override
        public InstanceDetails getDetails() {
            return details;
        }

        //

        private double[] partials;
        private double[] diffusions;

        private void updatePartial(
                final int kBuffer,
                final int iBuffer,
                final int iMatrix,
                final int jBuffer,
                final int jMatrix
        ) {

        }

        private void allocateStorage() {
            partials = new double[dimPartial * bufferCount];
            diffusions = new double[dimTrait * dimTrait * diffusionCount];
        }

        private String getOperationString(final int[] operations, final int offset) {
            StringBuilder sb = new StringBuilder("op:");
            for (int i = 0; i < ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE; ++i) {
                sb.append(" ").append(operations[offset + i]);
            }
            return sb.toString();
        }
    }
}
