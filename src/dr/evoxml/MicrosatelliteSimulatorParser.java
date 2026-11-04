/*
 * MicrosatelliteSimulatorParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

package dr.evoxml;

import dr.xml.*;
import dr.evolution.tree.Tree;
import dr.evolution.alignment.Patterns;
import dr.evolution.util.Taxa;
import dr.evolution.datatype.Microsatellite;
import dr.oldevomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.app.seqgen.MicrosatelliteSimulator;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for microsatellite simulator
 *
 */
public class MicrosatelliteSimulatorParser extends AbstractXMLObjectParser {
    public static final String MICROSATELLITE_SIMULATOR = "microsatelliteSimulator";

    public String getParserName(){
        return MICROSATELLITE_SIMULATOR;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Microsatellite msatDataType = (Microsatellite)xo.getChild(Microsatellite.class);
        Taxa taxa = (Taxa)xo.getChild(Taxa.class);
        Tree tree = (Tree)xo.getChild(Tree.class);
        MicrosatelliteModel msatModel = (MicrosatelliteModel)xo.getChild(MicrosatelliteModel.class);

        BranchRateModel brModel = (BranchRateModel)xo.getChild(BranchRateModel.class);
        if(brModel == null){
            brModel = new DefaultBranchRateModel();
        }

        MicrosatelliteSimulator msatSim = new MicrosatelliteSimulator(
                msatDataType,
                taxa,
                tree,
                new GammaSiteModel(msatModel),
                brModel
        );

        Patterns patterns =  msatSim.simulateMsatPattern();

        String msatPatId = xo.getAttribute("id","simMsatPat");
        patterns.setId(msatPatId);
        MicrosatellitePatternParser.printDetails(patterns);
        MicrosatellitePatternParser.printMicrosatContent(patterns);

        return patterns;

    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Microsatellite.class),
                new ElementRule(Taxa.class),
                new ElementRule(Tree.class),
                new ElementRule(MicrosatelliteModel.class),
                new ElementRule(BranchRateModel.class, true)

        };
    }



    public String getParserDescription(){
        return "This parser facilliates simulation of microsatellites " +
                "given a tree and infinitesimal rate model";
    }

    public Class getReturnType(){
        return Patterns.class;
    }
}
