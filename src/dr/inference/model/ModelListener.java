/*
 * ModelListener.java
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

package dr.inference.model;

import java.io.Serializable;

/**
 * An interface that provides a listener on a model.
 *
 * @version $Id: ModelListener.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */

public interface ModelListener extends Serializable {

	/**
	 * The model has changed. The model firing the event can optionally
	 * supply a reference to an object and an index if appropriate. Use
	 * of this extra information will be contingent on recognising what
	 * model it was that fired the event.
	 */
	void modelChangedEvent(Model model, Object object, int index);

    /**
     * The model has been restored.
     * Required only for notification of non-models (say pure likelihoods) which depend on
     * models.
     * @param model
     */
    void modelRestored(Model model);
}
