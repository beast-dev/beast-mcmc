/*
 * CDIFlags.java
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

public enum CDIFlag {
    PRECISION_SINGLE(1L, "double precision computation"),
    PRECISION_DOUBLE(2L, "single precision computation"),
    COMPUTATION_SYNCH(4L, "synchronous computation (blocking"),
    COMPUTATION_ASYNCH(8L, "asynchronous computation (non-blocking)"),
//    EIGEN_REAL(16L, "real eigenvalue computation"),
//    EIGEN_COMPLEX(32L, "complex eigenvalue computation"),
//    SCALING_MANUAL(64L, "manual scaling"),
//    SCALING_AUTO(128L, "auto-scaling on"),
//    SCALING_ALWAYS(256L, "scale at every update"),
//    SCALING_DYNAMIC(524288L, "manual scaling with dynamic checking"),
//    SCALERS_RAW(512L, "save raw scalers"),
//    SCALERS_LOG(1024L, "save log scalers"),
    VECTOR_SSE(2048L, "SSE vector computation"),
    VECTOR_NONE(4096L, "no vector computation"),
    THREADING_TBB(8192L, "TBB threading"),
    THREADING_NONE(16384L, "no threading"),
    PROCESSOR_CPU(32768L, "use CPU as main processor"),
    PROCESSOR_GPU(65536L, "use GPU as main processor"),
    PROCESSOR_FPGA(131072L, "use FPGA as main processor"),
    PROCESSOR_CELL(262144L, "use CELL as main processor"),
    FRAMEWORK_CUDA(4194304L, "use CUDA implementation with GPU resources"),
    FRAMEWORK_OPENCL(8388608L, "use OpenCL implementation with CPU or GPU resources"),
    FRAMEWORK_CPU(134217728L, "use CPU implementation");

    private final long mask;
    private final String meaning;

    CDIFlag(long mask, String meaning) {
        this.mask = mask;
        this.meaning = meaning;
    }

    public long getMask() {
        return this.mask;
    }

    public String getMeaning() {
        return this.meaning;
    }

    public boolean isSet(long bits) {
        return (bits & this.mask) != 0L;
    }

    public static String toString(long bits) {
        StringBuilder sb = new StringBuilder();
        CDIFlag[] flags = values();

        for (CDIFlag flag : flags) {
            if (flag.isSet(bits)) {
                sb.append(" ").append(flag.name());
            }
        }

        return sb.toString();
    }
}
