/*
 * CoalescentIntervalStatistic.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent;

import dr.inference.model.Statistic;


public class CoalescentIntervalStatistic extends Statistic.Abstract {

//	private coalescentInterval coalInt;

    private final CoalescentIntervalProvider coalescent;

    public CoalescentIntervalStatistic(CoalescentIntervalProvider coalescent) {
        this.coalescent = coalescent;
    }

//	public CoalescentIntervalStatistic(GMRFSkyrideLikelihood acl) {
//		this.coalInt = new GMRFStatistic(acl);
//	}
//
//	public CoalescentIntervalStatistic(CoalescentLikelihood coal) {
//		this.coalInt = new CoalescentStatistic(coal);
//	}

    public int getDimension() {
        return coalescent.getCoalescentIntervalDimension();
    }

    public double getStatisticValue(int i) {
        return coalescent.getCoalescentInterval(i);
    }

//	private class GMRFStatistic implements coalescentInterval {
//
//		private GMRFSkyrideLikelihood acl;
//		private int dimension;
//
//		private GMRFStatistic(GMRFSkyrideLikelihood acl) {
//			this.acl = acl;
//			this.dimension = acl.getCoalescentIntervalHeights().length;
//		}
//
//		public int getDimension() {
//			return dimension;
//		}
//
//		public double getStatisticValue(int dim) {
//			return acl.getCoalescentIntervalHeights()[dim];
//		}
//
//	}
//
//	private class CoalescentStatistic implements coalescentInterval {
//
//		private CoalescentLikelihood coal;
//		private int dimension;
//
//		private CoalescentStatistic(CoalescentLikelihood coal) {
//			this.coal = coal;
//			this.dimension = coal.getCoalescentIntervalHeights().length;
//		}
//
//		public int getDimension() {
//			return dimension;
//		}
//
//		public double getStatisticValue(int dim) {
//			return coal.getCoalescentIntervalHeights()[dim];
//		}
//
//	}
//
//	private interface coalescentInterval {
//
//		public int getDimension();
//
//		public double getStatisticValue(int dim);
//
//	}

}
