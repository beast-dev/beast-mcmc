/*
 * HastingsRatioContributor.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.model;

/**
 * A likelihood component that makes proposals of its own while it is evaluated (for example a data
 * likelihood that resamples latent state after a topology change), and so contributes to the Hastings
 * ratio of the MCMC step. The Markov chain adds the contribution to the operator's log Hastings ratio
 * after evaluating the posterior, before the acceptance test, so that it is not tempered.
 *
 * @author Andrew Rambaut
 */
public interface HastingsRatioContributor {

    /**
     * @return the log Hastings ratio (log of the reverse over the forward proposal density) of the proposals
     * made since the step began, and clears it (a second call in the same step returns 0)
     */
    double takeLogHastingsRatio();
}
