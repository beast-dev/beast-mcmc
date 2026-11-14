/*
 * EfficientBirthDeathSpeciationLikelihood.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.birthdeath.BirthDeathModel;
import dr.evomodel.birthdeath.EfficientBirthDeathLikelihood;

import java.util.Set;

/**
 * Compatibility wrapper for legacy speciationLikelihood XML that uses birth-death models.
 */
public class EfficientBirthDeathSpeciationLikelihood extends SpeciationLikelihood {

    private final EfficientBirthDeathLikelihood birthDeathLikelihood;

    public EfficientBirthDeathSpeciationLikelihood(Tree tree, BirthDeathModel birthDeathModel, Set<Taxon> exclude, String id) {
        super(tree, birthDeathModel, exclude, id);
        birthDeathLikelihood = new EfficientBirthDeathLikelihood(tree, birthDeathModel, exclude, id);
        addModel(birthDeathLikelihood);
    }

    public EfficientBirthDeathLikelihood getBirthDeathLikelihood() {
        return birthDeathLikelihood;
    }

    @Override
    double calculateLogLikelihood() {
        return birthDeathLikelihood.getLogLikelihood();
    }
}
