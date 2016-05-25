/*
 * PartitionDataList.java
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

package dr.app.bss;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.util.Taxa;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionDataList extends ArrayList<PartitionData> implements Serializable {

    public int simulationsCount = 1;
    public boolean useParallel = false;
    public boolean outputAncestralSequences = false;
    public SimpleAlignment.OutputType outputFormat = SimpleAlignment.OutputType.FASTA;

    //List of all Taxa displayed in Taxa Panel
    public Taxa allTaxa = new Taxa();
    public LinkedList<TreesTableRecord> recordsList = new LinkedList<TreesTableRecord>();

    // do not serialize this two
    public transient boolean setSeed = false;
    public transient long startingSeed;

    public PartitionDataList() {
        super();
        startingSeed = System.currentTimeMillis();
    }// END: Constructor

}// END:class
