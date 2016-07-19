/*
 * UncertainSiteList.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;

/**
 * Created by msuchard on 2016-05-19.
 */
public class UncertainSiteList extends SimpleSiteList implements Citable {

    public UncertainSiteList(DataType dataType, TaxonList taxonList) {
        super(dataType, taxonList);
    }

    public int addPattern(int[] pattern) {
        throw new IllegalArgumentException("Do not call directly");
    }

    public void addPattern(double[][] uncertainPattern) {
        uncertainSitePatterns.add(uncertainPattern);

        int[] map = new int[uncertainPattern.length];
        for (int i = 0; i < uncertainPattern.length; ++i) {
            map[i] = getMostProbable(uncertainPattern[i]);
        }
        super.addPattern(map);
    }

    private static int getMostProbable(double[] probabilities) {
        int map = 0;
        double maxProb = probabilities[0];

        for (int i = 1; i < probabilities.length; ++i) {
            if (probabilities[i] > maxProb) {
                maxProb = probabilities[i];
                map = i;
            }
        }
        return map;
    }

    public void fillPartials(final int sequenceIndex, final int site, double[] partials, final int offset) {
        double[][] sitePatterns = uncertainSitePatterns.get(site);
        System.arraycopy(sitePatterns[sequenceIndex], 0, partials, offset, getDataType().getStateCount());
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.DATA_MODELS;
    }

    @Override
    public String getDescription() {
        return "Uncertain site list";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(
                new Citation(
                        new Author[]{
                                new Author("MA", "Suchard"),
                                new Author("P", "Lemey"),
                                new Author("M", "Scotch"),

                        },
                        Citation.Status.IN_PREPARATION
                ));
    }

    private List<double[][]> uncertainSitePatterns = new ArrayList<double[][]>();
}
