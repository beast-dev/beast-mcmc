/*
 * OperatorFailedException.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.operators;

/**
 * This exception provides a fast-fail mechanism for operators.
 * It represents the situation where an operator realises rejection
 * is guaranteed to occur during acceptance phase.
 *
 * @author Alexei Drummond
 *
 * @version $Id: OperatorFailedException.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public class OperatorFailedException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3759460515533187879L;

	public OperatorFailedException(String message) {
		super(message);
	}
}

