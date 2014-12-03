/*
 * Prior.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.prior;

import dr.inference.model.Model;

import java.io.Serializable;

/**
 * This interface provides for general priors on models.
 *
 * @author Alexei Drummond
 * @version $Id: Prior.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public interface Prior extends Serializable {

    public static final class UniformPrior implements Prior {
        public double getLogPrior(Model m) {
            return 0.0;
        }

        public String getPriorName() {
            return "Uniform";
        }

        public String toString() {
            return "Uniform";
        }
    }

    public static final UniformPrior UNIFORM_PRIOR = new UniformPrior();


    /**
     * @param model the model under inspection.
     * @return the log prior of some aspect of the given model.
     */
    public double getLogPrior(Model model);

    /**
     * Returns the logical name of this prior. This name should be
     * the same as the string returned in the name attribute of the
     * XML prior.
     *
     * @return the logical name of this prior.
     */
    public String getPriorName();
}
