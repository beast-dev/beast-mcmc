/*
 * MPIServices.java
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

package dr.inference.parallel;

import mpi.MPI;
import mpi.Status;

/**
 * @author Marc A. Suchard
 */

public class MPIServices {

    public static void requestTermination(int slave) {
        int[] msg = new int[]{ServiceRequest.terminateProcess.getId()};
        MPI.COMM_WORLD.Send(msg,
                0, 1, MPI.INT, slave, ServiceRequest.MSG_REQUEST_TYPE);
    }

    public static void requestLikelihood(int slave) {
        int[] msg = new int[]{ServiceRequest.calculateLikeliood.getId()};
        MPI.COMM_WORLD.Send(msg, 0, 1, MPI.INT, slave, ServiceRequest.MSG_REQUEST_TYPE);
    }

    public static ServiceRequest getRequest(int master) {
        int[] msg = new int[1];
        Status status = MPI.COMM_WORLD.Recv(msg, 0, 1, MPI.INT,
                master, ServiceRequest.MSG_REQUEST_TYPE);
        // todo check status and throw exception if error
        return ServiceRequest.getByID(msg[0]);
    }

    public static void sendDouble(double value, int dest) {
        double[] msg = new double[]{value};
//		System.err.println("Sending double "+value+" to process "+dest);
        MPI.COMM_WORLD.Send(msg, 0, 1, MPI.DOUBLE, dest, ServiceRequest.MSG_REQUEST_TYPE);
    }

    public static void sendDoubleArray(double[] values, int dest) {
        MPI.COMM_WORLD.Send(values, 0, values.length, MPI.DOUBLE,
                dest, ServiceRequest.MSG_REQUEST_TYPE);

    }

    public static void sendDoubleArray(Double[] value, int dest) {
        double[] v = new double[value.length];
        for (int i = 0; i < v.length; i++) {
            v[i] = value[i];
        }
        sendDoubleArray(v, dest);
    }

    public static void sendIntArray(int[] values, int dest) {
        MPI.COMM_WORLD.Send(values, 0, values.length, MPI.INT,
                dest, ServiceRequest.MSG_REQUEST_TYPE);

    }

    public static void sendInt(int value, int dest) {
        int[] msg = new int[]{value};
        MPI.COMM_WORLD.Send(msg, 0, 1, MPI.INT, dest, ServiceRequest.MSG_REQUEST_TYPE);
    }

    public static int receiveInt(int source) {
        int[] msg = new int[1];
        Status status = MPI.COMM_WORLD.Recv(msg, 0, 1, MPI.INT,
                source, ServiceRequest.MSG_REQUEST_TYPE);
        return msg[0];
    }

    public static double receiveDouble(int source) {
        double[] msg = new double[1];
        Status status = MPI.COMM_WORLD.Recv(msg, 0, 1, MPI.DOUBLE,
                source, ServiceRequest.MSG_REQUEST_TYPE);
        return msg[0];
    }

    public static double[] receiveDoubleArray(int source, int length) {
        double[] msg = new double[length];
        Status status = MPI.COMM_WORLD.Recv(msg, 0, length, MPI.DOUBLE,
                source, ServiceRequest.MSG_REQUEST_TYPE);
        return msg;
    }

    public static int[] receiveIntArray(int source, int length) {
        int[] msg = new int[length];
        Status status = MPI.COMM_WORLD.Recv(msg, 0, length, MPI.INT,
                source, ServiceRequest.MSG_REQUEST_TYPE);
        return msg;
    }

    /*public static double receiveSumOfDoubleReduction(int master) {

     }*/

}
