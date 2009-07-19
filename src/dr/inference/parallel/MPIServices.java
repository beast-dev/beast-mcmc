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
