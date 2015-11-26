/*
 * PopsIOSpeciesTreePriorParser.java
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

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.evomodel.speciation.PopsIOSpeciesTreePrior;
import dr.evomodel.speciation.SpeciationModel;
import dr.xml.*;

/**
 * @author Graham Jones
 * Date: 10/05/12
 */
public class PopsIOSpeciesTreePriorParser extends AbstractXMLObjectParser {
    public static final String POPSIO_SPECIES_TREE_PRIOR = "PopsIOSpeciesTreePrior";
    public static final String MODEL = "model";
    public static final String PIO_TREE = "pioTree";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        System.out.println("PopsIOSpeciesTreePriorParser");
        final XMLObject mxo = xo.getChild(MODEL);
        final SpeciationModel sppm = (SpeciationModel) mxo.getChild(SpeciationModel.class);
        final XMLObject mulsptxo = xo.getChild(PIO_TREE);
        final PopsIOSpeciesTreeModel piostm = (PopsIOSpeciesTreeModel) mulsptxo.getChild(PopsIOSpeciesTreeModel.class);
        return new PopsIOSpeciesTreePrior(sppm, piostm);
    }


    private  XMLSyntaxRule[] speciationModelSyntax() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciationModel.class)
        };

    }

    private  XMLSyntaxRule[] piostmSyntax() {
        return new XMLSyntaxRule[]{
                new ElementRule(PopsIOSpeciesTreeModel.class)
        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MODEL, speciationModelSyntax()),
                new ElementRule(PIO_TREE, piostmSyntax()),

        };
    }

    @Override
    public String getParserDescription() {
        return "Prior for a species tree.";
    }

    @Override
    public Class getReturnType() {
        return PopsIOSpeciesTreePrior.class;
    }

    public String getParserName() {
        return POPSIO_SPECIES_TREE_PRIOR;
    }
}
