/*
 * SpatialParsimony.java
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

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.stats.DiscreteStatistics;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: adru001
 * Date: Jun 13, 2006
 * Time: 9:40:49 AM
 * To change this template use File | Settings | File Templates.
 */
public class SpatialParsimony {


    public static void analyzeParsimony(Tree tree, Map<String, SpaceTime> map) {

        Random random = new Random();

        List<SpaceTime> spaceTimes = new ArrayList<SpaceTime>();
        Taxa taxa = new Taxa();

        for (String name : map.keySet()) {
            SpaceTime spaceTime = map.get(name);
            spaceTime.setName(name);
            spaceTimes.add(spaceTime);
            taxa.addTaxon(new Taxon(name));
        }

        Collections.sort(spaceTimes, SpaceTime.getLatitudeComparator());

        Set<String> leafSet = new HashSet<String>();
        for (int i = 0; i < spaceTimes.size() - 1; i++) {

            SpaceTime spaceTime = spaceTimes.get(i);

            //System.out.println("Latitude = " + spaceTime.latitude);

            leafSet.add(spaceTime.getName());

            // null distribution
            double[] steps = new double[1000];
            for (int j = 0; j < steps.length; j++) {
                Set<String> leafSet2 = new HashSet<String>();
                List<SpaceTime> spaceTimes2 = new ArrayList<SpaceTime>(spaceTimes);
                while (leafSet2.size() < leafSet.size()) {

                    int index = random.nextInt(spaceTimes2.size());
                    SpaceTime st = spaceTimes2.remove(index);

                    leafSet2.add(st.getName());

                }

                //System.out.println(leafSet2.size());

                steps[j] = Tree.Utils.getParsimonySteps(tree, leafSet2);

            }


            int parsimonySteps = Tree.Utils.getParsimonySteps(tree, leafSet);

            System.out.println(spaceTime.getName() + "\t" +
                    spaceTime.latitude + "\t" +
                    parsimonySteps + "\t" +
                    DiscreteStatistics.mean(steps) + "\t" +
                    DiscreteStatistics.quantile(0.025, steps) + "\t" +
                    DiscreteStatistics.quantile(0.975, steps) + "\t");
        }

    }


}
