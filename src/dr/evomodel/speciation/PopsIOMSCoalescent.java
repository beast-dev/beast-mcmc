/*
 * PopsIOMSCoalescent.java
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

package dr.evomodel.speciation;

import dr.inference.model.Likelihood;

/**
 * @author  Graham  Jones
 * Date: 10/05/2012
 */

public class PopsIOMSCoalescent extends Likelihood.Abstract {
    private final PopsIOSpeciesTreeModel piostm;
    private final PopsIOSpeciesBindings piosb;


    public PopsIOMSCoalescent(PopsIOSpeciesBindings piosb, PopsIOSpeciesTreeModel piostm) {
        super(piostm);
        this.piostm = piostm;
        this.piosb = piosb;

        piostm.addModelRestoreListener(this);

        final PopsIOSpeciesBindings.GeneTreeInfo[] trees = piosb.getGeneTrees();
        for(PopsIOSpeciesBindings.GeneTreeInfo geneTree : trees) {
            geneTree.tree.addModelListener(this);
        }

    }


    @Override
    protected double calculateLogLikelihood() {
for (int i = 0; i < piosb.numberOfGeneTrees(); i++) {
            if (!piosb.geneTreeFitsInNetwork(i, piostm)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        // grjtodo-oneday JH has compatible flags for efficiency. I'm checking
        // every time.

        piosb.fillSpeciesTreeWithCoalescentInfo(piostm);
        return piostm.logLhoodAllGeneTreesInSpeciesTree();
    }


    @Override
    protected boolean getLikelihoodKnown() {
        return false;
    }
}
