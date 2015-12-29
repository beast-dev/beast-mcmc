/*
 * ServiceRequest.java
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
