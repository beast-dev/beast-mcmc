/*
 * CalibratedSpeciationLikelihoodParser.java
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

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.speciation.CalibratedSpeciationLikelihood;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.tree.TMRCAStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.SkewTDistribution;
import dr.inference.distribution.SoftBoundUniformDistribution;
import dr.inference.distribution.TruncatedSoftBoundCauchyDistribution;
import dr.math.distributions.Distribution;
import dr.xml.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CalibratedSpeciationLikelihoodParser extends AbstractXMLObjectParser {

    public static final String CALIBRATED_SPECIATION_LIKELIHOOD = "calibratedSpeciationLikelihood";
    public static final String CALIBRATION = "calibration";
    public static final String MCMCTREE_NEWICK = "MCMCTreeNewick";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SpeciationLikelihood speciationLikelihood = (SpeciationLikelihood) xo.getChild(SpeciationLikelihood.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        List<CalibratedSpeciationLikelihood.CalibrationLikelihood> calibrationLikelihoods = new ArrayList<>();
        if (xo.hasChildNamed(CALIBRATION)) {
            for (XMLObject cxo : xo.getAllChildren(CALIBRATION)) {
                TMRCAStatistic tmrcaStatistic = (TMRCAStatistic) cxo.getChild(TMRCAStatistic.class);
                Distribution distribution = (Distribution) cxo.getChild(Distribution.class);
                calibrationLikelihoods.add(new CalibratedSpeciationLikelihood.CalibrationLikelihood(tmrcaStatistic, distribution));
            }
        } else {
            String labeledNewickString = (String) xo.getChild(MCMCTREE_NEWICK).getStringChild(0);
            calibrationLikelihoods = parseMCMCTreeString(labeledNewickString, tree);
        }

        return new CalibratedSpeciationLikelihood(CALIBRATED_SPECIATION_LIKELIHOOD, speciationLikelihood, tree, calibrationLikelihoods);
    }

    private List<CalibratedSpeciationLikelihood.CalibrationLikelihood> parseMCMCTreeString(String labeledNewickString, TreeModel tree) throws XMLParseException {
        List<CalibratedSpeciationLikelihood.CalibrationLikelihood> calibrationLikelihoods = new ArrayList<>();
        NewickImporter newickImporter = new NewickImporter(labeledNewickString);
        FlexibleTree flexibleTree;
        try {
            flexibleTree = (FlexibleTree) newickImporter.importTree(null);
        } catch (IOException ioe) {
            throw new XMLParseException("error parsing tree in newick element");
        } catch (NewickImporter.BranchMissingException bme) {
            throw new XMLParseException("branch missing in tree in newick element");
        } catch (Importer.ImportException ime) {
            throw new XMLParseException("error parsing tree in newick element - " + ime.getMessage());
        }

        for (int i = 0; i < flexibleTree.getInternalNodeCount(); i++) {
            FlexibleNode node = (FlexibleNode) flexibleTree.getNode(flexibleTree.getExternalNodeCount() + i);
//            TreeUtils.getDescendantLeaves(node)
            if (node.getAttribute("label") != null) {
                String distributionString = (String) node.getAttribute("label");
                Distribution distribution = parseDistributionString(distributionString);
                try {
                    TMRCAStatistic tmrcaStatistic = constructTMRCAStatistics(node, flexibleTree, tree);
                    calibrationLikelihoods.add(new CalibratedSpeciationLikelihood.CalibrationLikelihood(tmrcaStatistic, distribution));
                } catch (TreeUtils.MissingTaxonException e) {
                    throw new RuntimeException("Non-compatible taxa set.");
                }
            }
        }
        return calibrationLikelihoods;
    }

    private TMRCAStatistic constructTMRCAStatistics(FlexibleNode node, FlexibleTree flexibleTree, Tree tree) throws TreeUtils.MissingTaxonException {
        String name = "tmrca(node." + String.valueOf(node.getNumber()) + ")";
        Set<String> leaves = TreeUtils.getDescendantLeaves(flexibleTree, node);
        Collection<Taxon> taxa = new ArrayList<>();
        for (String taxonID : leaves) {
            taxa.add(tree.getTaxon(tree.getTaxonIndex(taxonID)));
        }
        return new TMRCAStatistic(name, tree, new Taxa(taxa), false, false);
    }

    private final Pattern distributionPattern = Pattern.compile("(\\D+)\\(([\\d.,-]+)\\)");

    private Distribution parseDistributionString(String distributionString) {
        Matcher MCMCTreeDistribution = distributionPattern.matcher(distributionString.replaceAll("\\s", ""));
        if (MCMCTreeDistribution.find()) {
            String name = MCMCTreeDistribution.group(1);
            String numbers = MCMCTreeDistribution.group(2);
            return constructDistribution(name, numbers);
        } else {
            throw new RuntimeException("Invalid distribution string.");
        }
    }

    private Distribution constructDistribution(String distributionName, String distributionParameters) {
        String[] parameters = distributionParameters.split(",");
        Distribution distribution;
        if (distributionName.equals("B")) {
            distribution = new SoftBoundUniformDistribution(Double.parseDouble(parameters[0]), Double.parseDouble(parameters[1]),
                    Double.parseDouble(parameters[2]), Double.parseDouble(parameters[3]));
        } else if (distributionName.equals("ST")) {
            distribution = new SkewTDistribution(Double.parseDouble(parameters[0]), Double.parseDouble(parameters[1]),
                    Double.parseDouble(parameters[2]), Double.parseDouble(parameters[3]));
        } else if (distributionName.equals("L")) {
            if (parameters.length == 3) {
                distribution = new TruncatedSoftBoundCauchyDistribution(Double.parseDouble(parameters[0]), Double.parseDouble(parameters[1]),
                        Double.parseDouble(parameters[2]));
            } else {
                distribution = new TruncatedSoftBoundCauchyDistribution(Double.parseDouble(parameters[0]), Double.parseDouble(parameters[1]),
                        Double.parseDouble(parameters[2]), Double.parseDouble(parameters[3]));
            }
        } else {
            throw new RuntimeException("Not yet implemented");
        }
        return distribution;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SpeciationLikelihood.class),
            new ElementRule(TreeModel.class),
            new XORRule(
                    new ElementRule(CALIBRATION,
                    new XMLSyntaxRule[] {
                            new ElementRule(TMRCAStatistic.class),
                            new ElementRule(Distribution.class)
                    }, 1, Integer.MAX_VALUE),
                    new ElementRule(MCMCTREE_NEWICK, String.class, "The MCMCTree style NEWICK format tree with calibrated internal nodes. Tip labels are taken to be Taxon IDs"))
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CalibratedSpeciationLikelihood.class;
    }

    @Override
    public String getParserName() {
        return CALIBRATED_SPECIATION_LIKELIHOOD;
    }
}
