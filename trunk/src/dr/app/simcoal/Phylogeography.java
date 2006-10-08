/*
 * Phylogeography.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.simcoal;

import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.util.Taxon;
import dr.evolution.continuous.Continuous;
import dr.evolution.continuous.ContinuousTraitLikelihood;
import dr.evolution.continuous.Contrastable;
import dr.evolution.continuous.DiffusionLikelihood;
import dr.stats.DiscreteStatistics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 */
public class Phylogeography {


    // load location and time information
    public static Map<String, SpaceTime> loadBiekFile(String fileName) throws IOException {

        SortedMap<String, SpaceTime> map = new TreeMap<String, SpaceTime>();


        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        reader.readLine();

        String line = reader.readLine();
        while (line != null) {
            // process a single line

            StringTokenizer tokens = new StringTokenizer(line, "\t");
            String id = tokens.nextToken();
            String fullName = tokens.nextToken();
            tokens.nextToken();
            String latStr = tokens.nextToken();
            String longStr = tokens.nextToken();

            double latitude = Double.parseDouble(latStr);
            double longitude = Double.parseDouble(longStr);
            double time = Double.parseDouble(fullName.substring(fullName.indexOf('_') + 1));

            SpaceTime spaceTime = new SpaceTime(latitude, longitude, time);

            map.put(id, spaceTime);

            line = reader.readLine();
        }

        reader.close();

        return map;
    }

    public static String createBEASTTaxaBlock(Map<String, SpaceTime> map) {

        StringBuilder builder = new StringBuilder();

        builder.append("<!-- ntax=" + map.size() + " -->\n");
        builder.append("<taxa id=\"taxa\">\n");
        for (String key : map.keySet()) {

            SpaceTime spaceTime = map.get(key);

            builder.append("  <taxon id=\"" + key + "\">\n");
            builder.append("    <date value=\"" + spaceTime.time + "\" direction=\"forwards\" units=\"years\"/>\n");
            builder.append("  </taxon>\n");

        }
        builder.append("</taxa>\n");

        return builder.toString();
    }

    public static String createBEASTTaxaBlock(Map<String, SpaceTime> map, PopAssigner assigner, int pop, String name) {

        StringBuilder builder = new StringBuilder();

        builder.append("<taxa id=\"" + name + "\">\n");
        for (String key : map.keySet()) {

            SpaceTime spaceTime = map.get(key);

            if (assigner.getPopulation(spaceTime) == pop) {

                builder.append("  <taxon idref=\"" + key + "\"/>\n");
            }
        }
        builder.append("</taxa>\n");

        return builder.toString();
    }


    public static double[] processSpaceTimeTree(Map<String, SpaceTime> map, Tree tree, PopAssigner assigner) {

        int n = tree.getExternalNodeCount();
        double[] params = null;

        List<Double> times = new ArrayList<Double>();
        List<Double> pop1 = new ArrayList<Double>();
        List<Double> pop2 = new ArrayList<Double>();

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getTaxon(i);
            String id = taxon.getId();

            SpaceTime spaceTime1 = map.get(id);
            double index1 = assigner.getPopulation(spaceTime1);

            for (int j = i+1; j < tree.getExternalNodeCount(); j++) {
                Taxon taxon2 = tree.getTaxon(j);
                String id2 = taxon2.getId();

                SpaceTime spaceTime2 = map.get(id2);
                double index2 = assigner.getPopulation(spaceTime2);

                double time = getPatristicDistance(tree, tree.getExternalNode(i), tree.getExternalNode(j));

                times.add(time);
                pop1.add(index1);
                pop2.add(index2);
            }
        }


        double[] timeArray = new double[times.size()];
        double[] pop1Array = new double[pop1.size()];
        double[] pop2Array = new double[pop2.size()];
        for (int i = 0; i < timeArray.length; i++) {
            timeArray[i] = times.get(i);
            pop1Array[i] = pop1.get(i);
            pop2Array[i] = pop2.get(i);
        }

        double rSquared = rSquared(timeArray, pop1Array);
        double rSquared2 = rSquared(pop1Array, pop2Array);

        //System.out.println("r-squared = " + rSquared + " " + rSquared2);

        Notohara1990Function lsf = new Notohara1990Function(timeArray,pop1Array, pop2Array);
        params = lsf.optimize();

        //System.out.println(cutoff + "\t" + params[0] + "\t" + params[1] + "\t" + rSquared);

        return params;


    }

    static double rSquared(double[] x, double[] y) {

        return DiscreteStatistics.covariance(x, y);

    }

    static double getPatristicDistance(Tree tree, NodeRef node1, NodeRef node2) {
        HashSet set = new HashSet();
        set.add(tree.getTaxon(node1.getNumber()).getId());
        set.add(tree.getTaxon(node2.getNumber()).getId());

        NodeRef ancestor = Tree.Utils.getCommonAncestorNode(tree, set);

        double distance = 0.0;

        NodeRef node = node1;
        while (node != ancestor && !tree.isRoot(node)) {
            distance += tree.getBranchLength(node);
            node = tree.getParent(node);
        }
        node = node2;
        while (node != ancestor && !tree.isRoot(node)) {
            distance += tree.getBranchLength(node);
            node = tree.getParent(node);
        }
        return distance;
    }

    public static Tree analyzeContinuousTraits(Map<String, SpaceTime> map, Tree inTree) {

        MutableTree tree = new FlexibleTree(inTree);

        for (String name: map.keySet()) {
            SpaceTime spaceTime = map.get(name);

            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                if (name.equals(tree.getTaxon(i).getId())) {
                    tree.setTaxonAttribute(i,"latitude",new Continuous(spaceTime.latitude));
                    tree.setTaxonAttribute(i,"longitude",new Continuous(spaceTime.longitude));
                    tree.setTaxonAttribute(i,"coord",spaceTime.coord);
                }
            }

        }

        ContinuousTraitLikelihood ctLikelihood = new ContinuousTraitLikelihood();

        Contrastable[] mles = new Contrastable[2];

        System.out.println("kappa\tlogL\tlatitude\tlongitude");
        for (double kappa = 0.01; kappa < 2.0; kappa *= 1.1) {
            double logL = ctLikelihood.calculateLikelihood(tree, new String[] {"latitude", "longitude"}, mles, kappa);
            System.out.println(kappa + "\t" + logL + "\t" + mles[0] + "\t" + mles[1]);
        }


        Contrastable[] mle = new Contrastable[1];
        System.out.println("logL (latitude) = " + ctLikelihood.calculateLikelihood(tree, new String[] {"latitude"}, mle, 1.0));
        System.out.println("mle(latitude) = " + mle[0]);
        System.out.println("logL (longitude) = " + ctLikelihood.calculateLikelihood(tree, new String[] {"longitude"}, mle, 1.0));
        System.out.println("mle(longitude) = " + mle[0]);
        System.out.println("logL (coord) = " + ctLikelihood.calculateLikelihood(tree, new String[] {"coord"}, mle, 1.0));
        System.out.println("mle(coord) = " + mle[0]);

        DiffusionLikelihood lkl = new DiffusionLikelihood(tree, "latitude");
        double[] D = new double[1];
        double logL2 = lkl.optimize(D);
        System.out.println("latitude D=" + D[0] + " logL= " + logL2);
        //logL2 = lkl.optimize(Dbias);
        //System.out.println("latitude D=" + Dbias[0] + " bias=" + Dbias[1] + " logL= " + logL2);

        lkl = new DiffusionLikelihood(tree, "longitude");
        logL2 = lkl.optimize(D);
        System.out.println("longitude D=" + D[0] + " logL= " + logL2);
        //logL2 = lkl.optimize(Dbias);
        //System.out.println("longitude D=" + Dbias[0] + " bias=" + Dbias[1] + " logL= " + logL2);

        lkl = new DiffusionLikelihood(tree, "coord");
        logL2 = lkl.optimize(D);
        System.out.println("coord D=" + D[0]+ " logL= " + logL2);
        //logL2 = lkl.optimize(Dbias);
        //System.out.println("coord D=" + Dbias[0] + " bias=" + Dbias[1] + " logL= " + logL2);

        return tree;
    }


    public static void main(String[] args) throws IOException, Importer.ImportException {

        Map<String, SpaceTime> map = loadBiekFile(args[0]);

        double latitudeCutoff = 49.5;
        PopAssigner latitudeBarrier = new LatitudeBasedPopAssigner(latitudeCutoff);

        System.out.println(createBEASTTaxaBlock(map));
        System.out.println(createBEASTTaxaBlock(map,latitudeBarrier,0,"before" + latitudeCutoff));
        System.out.println(createBEASTTaxaBlock(map,latitudeBarrier,1,"after" + latitudeCutoff));

        NexusImporter importer = new NexusImporter(new FileReader(args[1]));

        Tree tree = importer.importNextTree();


        analyzeContinuousTraits(map, tree);

        System.exit(0);

        FlexibleTree treeCopy = new FlexibleTree(tree);
        for (int i = 0; i < treeCopy.getExternalNodeCount(); i++) {
            String name = treeCopy.getTaxonId(i);
            int num = Integer.parseInt(name);
            if (num>42) num -= 1;
            treeCopy.setTaxonId(i, num+"");
        }
        System.out.println(Tree.Utils.newick(treeCopy));

        SpatialParsimony.analyzeParsimony(tree, map);


        double[] params = processSpaceTimeTree(map, tree, latitudeBarrier);

        System.out.print(latitudeCutoff + "\t");
        for (int i = 0; i < params.length; i++) {
            System.out.print(params[i] + "\t");
        }
        System.out.println();
    }

}
