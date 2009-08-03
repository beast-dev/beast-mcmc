package dr.inference.parallel;

/**
 * @author Marc A. Suchard
 */

public enum ServiceRequest {

	calculateLikeliood, terminateProcess;

	public static final int CALCULATION_LIKELIHOOD_ID = 1;
	public static final int TERMINATE_ID = 2;
	public static final int NONE_ID = 0;


	public static final int MSG_REQUEST_TYPE = 10;

	public int getId() {
		switch (this) {
			case calculateLikeliood:
				return CALCULATION_LIKELIHOOD_ID;
			case terminateProcess:
				return TERMINATE_ID;
			default:
				return NONE_ID;
		}
	}

	public static ServiceRequest getByID(int id) {
		switch (id) {
			case CALCULATION_LIKELIHOOD_ID:
				return calculateLikeliood;
			case TERMINATE_ID:
				return terminateProcess;
			default:
				return null;
		}
	}

}
