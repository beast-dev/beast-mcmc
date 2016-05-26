/*
 * LikelihoodCoreFactory.java
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

/**
 *
 */
package dr.evomodel.newtreelikelihood;

import dr.app.beagle.evomodel.treelikelihood.AbstractTreeLikelihood;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 *
 */
public class LikelihoodCoreFactory {

	public static final boolean useFloats = false;

	public static LikelihoodCore loadLikelihoodCore(int[] configuration, AbstractTreeLikelihood treeLikelihood) {

		if (coreRegistry == null) {  // Lazy loading
			coreRegistry = new ArrayList<LikelihoodCoreLoader>();  // List libraries in order of load-priority
            coreRegistry.add(new ThreadedGPULikelihoodCore.LikelihoodCoreLoader());
            coreRegistry.add(new GPULikelihoodCore.LikelihoodCoreLoader());
            coreRegistry.add(new NativeLikelihoodCore.LikelihoodCoreLoader());
		}

		for(LikelihoodCoreLoader loader: coreRegistry) {
            System.out.print("Attempting to load core: " + loader.getLibraryName());
            LikelihoodCore core = loader.createLikelihoodCore(configuration, treeLikelihood);
			if (core != null) {
                System.out.println(" - SUCCESS");

				return core;
            }
            System.out.println(" - FAILED");
		}

		// No libraries/processes available

		int stateCount = configuration[0];
		if (useFloats)
			return new FloatGeneralLikelihoodCore(stateCount);
		else
			return new GeneralLikelihoodCore(stateCount);
	}

	private static List<LikelihoodCoreLoader> coreRegistry;

	protected interface LikelihoodCoreLoader {
		public String getLibraryName();

		/**
		 * Actual factory
		 * @param configuration
		 * @param treeLikelihood  Should remove this after migration bug is solved for GPU
		 * @return
		 */
		public LikelihoodCore createLikelihoodCore(int[] configuration, AbstractTreeLikelihood treeLikelihood);
	}

}
