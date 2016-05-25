/*
 * PopsIOSpeciesTreeModelParser.java
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

import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Graham Jones
 * Date: 20/09/2013
 */


public class PopsIOSpeciesTreeModelParser extends AbstractXMLObjectParser {
    public static final String PIO_SPECIES_TREE = "pioSpeciesTree";
    public static final String PIO_POP_PRIOR_SCALE = "pioPopPriorScale";
    public static final String PIO_POP_PRIOR_INVGAMMAS = "pioPopPriorInvGammas";
    public static final String PIO_POP_PRIOR_COMPONENT = "pioPopPriorComponent";
    public static final String WEIGHT = "weight";
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        System.out.println("PopsIOSpeciesTreeModelParser");
        PopsIOSpeciesBindings piosb = (PopsIOSpeciesBindings) xo.getChild(PopsIOSpeciesBindings.class);

        final Parameter popPriorScale = (Parameter) xo.getElementFirstChild(PIO_POP_PRIOR_SCALE);
        final XMLObject pioppxo = xo.getChild(PIO_POP_PRIOR_INVGAMMAS);
        final int nComponents = pioppxo.getChildCount();
        PopsIOSpeciesTreeModel.PriorComponent [] components = new PopsIOSpeciesTreeModel.PriorComponent[nComponents];

        for (int nc = 0; nc < nComponents; ++nc) {
            Object child = pioppxo.getChild(nc);
            assert ((XMLObject) child).getName().equals(PIO_POP_PRIOR_COMPONENT);
            double weight = ((XMLObject) child).getDoubleAttribute(WEIGHT);
            double alpha = ((XMLObject) child).getDoubleAttribute(ALPHA);
            double beta = ((XMLObject) child).getDoubleAttribute(BETA);
            components[nc] = new PopsIOSpeciesTreeModel.PriorComponent(weight, alpha, beta);
        }


        PopsIOSpeciesTreeModel piostm = new PopsIOSpeciesTreeModel(piosb, popPriorScale, components);
        return piostm;
    }


    private XMLSyntaxRule[] priorComponentSyntax() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newDoubleRule(ALPHA),
                AttributeRule.newDoubleRule(BETA)
        };
    }

    private XMLSyntaxRule[] popPriorSyntax() {
        return new XMLSyntaxRule[]{
             new ElementRule(PIO_POP_PRIOR_COMPONENT, priorComponentSyntax(), 1, Integer.MAX_VALUE)
        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
              new ElementRule(PopsIOSpeciesBindings.class),
              new ElementRule(PIO_POP_PRIOR_SCALE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
              new ElementRule(PIO_POP_PRIOR_INVGAMMAS, popPriorSyntax())
        };
    }

    @Override
    public String getParserDescription() {
        return "Species tree with population size parameters on branches integrated out analytically";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return PopsIOSpeciesTreeModel.class;
    }

    public String getParserName() {
        return PIO_SPECIES_TREE;
    }
}
