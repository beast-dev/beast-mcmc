/*
 * SimpleMCMCOperator.java
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

public abstract class SimpleMCMCOperator implements MCMCOperator {

	public double getTargetAcceptanceProbability() { return 0.234; }
	public double getMinimumAcceptanceLevel() { return 0.10; }
	public double getMaximumAcceptanceLevel() { return 0.40; }
	public double getMinimumGoodAcceptanceLevel() { return 0.20; }
	public double getMaximumGoodAcceptanceLevel() { return 0.30; }

	public abstract String getOperatorName();

	/**
	 * @return the weight of this operator.
	 */
	public int getWeight() { return weight; }

	/** 
	 * Sets the weight of this operator.
	 */
	public void setWeight(int w) {
		if (w > 0) {
			weight = w;
		} else throw new IllegalArgumentException("Weight must be a positive integer."); 
	}

	public final void accept(double deviation) {
		lastDeviation = deviation;
			
		if (!operateAllowed) {
			operateAllowed = true;
			accepted += 1;
			sumDeviation += deviation;
		} else throw new RuntimeException("Accept/reject methods called twice without operate called in between!");
	}
	
	public void reject() {
		if (!operateAllowed) {
			operateAllowed = true;
			rejected += 1;
		} else throw new RuntimeException("Accept/reject methods called twice without operate called in between!");	
	}
	
	public final void reset() {
		operateAllowed = true;
		accepted = 0;
		rejected = 0;
        lastDeviation = 0.0;
        sumDeviation = 0.0;
	}
	
	public final int getAccepted() { return accepted; }

    public final void setAccepted(int accepted) {
        this.accepted = accepted;
    }

    public final int getRejected() { return rejected; }

    public final void setRejected(int rejected) {
        this.rejected = rejected;
    }

    public final double getMeanDeviation() { return sumDeviation / accepted; }
	
	public final double getDeviation() {
		return lastDeviation;
	}

    public final double getSumDeviation() {
        return sumDeviation;
    }

    public final void setDumDeviation(double sumDeviation) {
        this.sumDeviation = sumDeviation;
    }

	public final double operate() throws OperatorFailedException {
		if (operateAllowed) {
			operateAllowed = false;
			return doOperation();
		} else throw new RuntimeException("Operate called twice without accept/reject in between!");
	}
	
	public final double getAcceptanceProbability() {
		return (double)accepted / (double)(accepted + rejected);
	}
	
	/**
	 * Called by operate(), does the actual operation. 
	 * @return the hastings ratio
     * @throws OperatorFailedException
	 */
	public abstract double doOperation() throws OperatorFailedException;

	protected int weight = 1;
	private int accepted = 0;
	private int rejected = 0;
	private double sumDeviation = 0.0;
	private double lastDeviation = 0.0;
	private boolean operateAllowed = true;
}
