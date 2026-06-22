/*
 * ContinuousCanonicalPostOrderMessageProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

package dr.evomodel.treedatalikelihood.continuous.message;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

/**
 * Provider for post-order canonical Gaussian messages on a continuous-trait tree.
 *
 * <p>All methods write into caller-owned output states. Exact identity tip
 * observations are degenerate Gaussian factors, so branch-bottom queries at such
 * tips are not guaranteed to be representable as a finite
 * {@link CanonicalGaussianState}. Branch-top messages remain finite after
 * propagation through a non-zero Gaussian branch transition.</p>
 */
public interface ContinuousCanonicalPostOrderMessageProvider extends ContinuousMessageProvider {

    void ensurePostOrderMessages();

    void fillPostOrderNode(int nodeIndex, CanonicalGaussianState out);

    void fillPostOrderBranchBottom(int childNodeIndex, CanonicalGaussianState out);

    void fillPostOrderBranchTop(int childNodeIndex, CanonicalGaussianState out);
}
