/*
 * ObjectStore.java
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

package dr.xml;

import dr.util.Identifiable;

import java.util.Collection;
import java.util.Set;

public interface ObjectStore {

    /**
     * @return the object with unique id or throws an ObjectNotFoundException
     */
    Object getObjectById(Object uid) throws ObjectNotFoundException;

    /**
     * @return true if an object with the given id exists in this ObjectStore.
     */
    boolean hasObjectId(Object uid);

    /**
     * Adds an object using the id returned by getId().
     *
     * @param force true if object should be placed in store even if it will replace
     *              an existing object of the same id, false otherwise.
     */
    void addIdentifiableObject(Identifiable object, boolean force);

    /**
     * @return a set of the UIDs of the objects in this store.
     */
    public Set getIdSet();

    /**
     * @return a collection of the objects in this store.
     */
    public Collection getObjects();
}
