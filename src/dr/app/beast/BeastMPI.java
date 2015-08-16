/*
 * BeastMPI.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beast;

/**
 * @author Marc A. Suchard
 */

public class BeastMPI {

    private static final String msg = "Unable to load mpiJava or MPJ";

//    private static Method getMPIMethod(String className, String methodName) {
//        Method method = null;
//        try {
//            Class clazz = Class.forName(className);
//            method = clazz.getMethod(methodName, Class.class);
//        } catch (Exception e) {
//            System.err.println(e.getMessage());
//            throw new RuntimeException("Libraries MPJ or mpiJava do not appear in the classpath");
//        }
//        return method;
//    }

    public static void Init(String[] args) {
        try {
            mpi.MPI.Init(args);
        } catch (NoClassDefFoundError e) {
            throw new RuntimeException(msg);
        }
    }

    public static void Finalize() {
        try {
            mpi.MPI.Finalize();
        } catch (NoClassDefFoundError e) {
            throw new RuntimeException(msg);
        }
    }

    public static class COMM_WORLD {
        public static int Rank() {
            int rtnValue = -1;
            try {
             rtnValue = mpi.MPI.COMM_WORLD.Rank();
            } catch (NoClassDefFoundError e) {
                throw new RuntimeException(msg);
            }
            return rtnValue;
        }
    }
}
