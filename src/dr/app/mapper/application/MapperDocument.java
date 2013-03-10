/*
 * DataSet.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.mapper.application;

import dr.evolution.util.Taxon;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * A class to store the current mapper document (i.e. all the data)
 *
 * @author Andrew Rambaut
 */
public final class MapperDocument {

    public MapperDocument() {
    }

    public void addTaxa(Collection<Taxon> newTaxa) {
        taxa.addAll(newTaxa);

        fireTaxaChanged();
    }

    public List<Taxon> getTaxa() {
        return taxa;
    }

    List<Taxon> taxa = new ArrayList<Taxon>();

    // Listeners and broadcasting
    public void addListener(Listener listener) {
        listeners.add(listener);
    }

    public void removeListener(Listener listener) {
        listeners.remove(listener);
    }

    public void fireTaxaChanged() {
        for (Listener listener : listeners) {
            listener.taxaChanged();
        }
    }

    private final List<Listener> listeners = new ArrayList<Listener>();

    public interface Listener {
        void taxaChanged();
    }

}
