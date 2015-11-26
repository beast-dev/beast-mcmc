/*
 * TreeSpaceDocument.java
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

package dr.app.treespace;

import dr.evolution.io.*;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.AbstractCladeImportanceDistribution;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class TreeSpaceDocument {
    private final List<InputFile> inputFiles = new ArrayList<InputFile>();
    private final CladeSystem cladeSystem = new CladeSystem();

    public void addTreeFile(InputFile inputFile) {
        inputFiles.add(inputFile);
        fireDataChanged();
    }

    public List<InputFile> getInputFiles() {
        return inputFiles;
    }

    public CladeSystem getCladeSystem() {
        return cladeSystem;
    }

    public void fireDataChanged() {
        for (Listener listener : listeners) {
            listener.dataChanged();
        }
    }

    public void fireSettingsChanged() {
        for (Listener listener : listeners) {
            listener.settingsChanged();
        }
    }

    public void addListener(Listener listener) {
        listeners.add(listener);
    }

    public interface Listener {
        void dataChanged();
        void settingsChanged();
    }

    private final Set<Listener> listeners = new HashSet<Listener>();
}
