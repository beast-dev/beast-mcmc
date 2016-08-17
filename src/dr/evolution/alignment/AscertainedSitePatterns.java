/*
 * AscertainedSitePatterns.java
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

package dr.evolution.alignment;

import dr.evolution.util.TaxonList;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * Package: AscertainedSitePatterns
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 10, 2008
 * Time: 12:50:36 PM
 */
public class AscertainedSitePatterns extends SitePatterns implements Citable {

    protected int[] includePatterns;
    protected int[] excludePatterns;
    protected int ascertainmentIncludeCount;
    protected int ascertainmentExcludeCount;

    /**
     * Constructor
     */
    public AscertainedSitePatterns(Alignment alignment) {
        super(alignment);
    }

    /**
     * Constructor
     */
    public AscertainedSitePatterns(Alignment alignment, TaxonList taxa) {
        super(alignment, taxa);
    }

    /**
     * Constructor
     */
    public AscertainedSitePatterns(Alignment alignment, int from, int to, int every) {
        super(alignment, from, to, every);
    }

    /**
     * Constructor
     */
    public AscertainedSitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every) {
        super(alignment, taxa, from, to, every);
    }

    /**
     * Constructor
     */
    public AscertainedSitePatterns(SiteList siteList) {
        super(siteList);
    }

    /**
     * Constructor
     */
    public AscertainedSitePatterns(SiteList siteList, int from, int to, int every) {
        super(siteList, from, to, every);
    }

    public AscertainedSitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every,
                                   int includeFrom, int includeTo,
                                   int excludeFrom, int excludeTo) {
        super(alignment, taxa, from, to, every);
        int[][] newPatterns = new int[patterns.length + (includeTo - includeFrom) + (excludeTo - excludeFrom)][];
        double[] newWeights = new double[patterns.length + (includeTo - includeFrom) + (excludeTo - excludeFrom)];
        for (int i = 0; i < patterns.length; ++i) {
            newPatterns[i] = patterns[i];
            newWeights[i] = weights[i];
        }
        patterns = newPatterns;
        weights = newWeights;

        if (includeTo - includeFrom >= 1)
            includePatterns(includeFrom, includeTo, every);
        if (excludeTo - excludeFrom >= 1)
            excludePatterns(excludeFrom, excludeTo, every);

    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.DATA_MODELS;
    }

    @Override
    public String getDescription() {
        return "Ascertained Site Patterns";
    }

    public List<Citation> getCitations() {
        return Collections.singletonList(
                CommonCitations.ALEKSEYENKO_2008
        );
    }

    public int getIncludePatternCount() {
        return ascertainmentIncludeCount;
    }

    public int[] getIncludePatternIndices() {
        return includePatterns;
    }

    protected void includePatterns(int includeFrom, int includeTo, int every) {
        if (includePatterns == null) {
            includePatterns = new int[includeTo - includeFrom];
        }
        for (int i = includeFrom; i < includeTo; i += every) {
            int[] pattern = siteList.getPattern(i);
            int index = addAscertainmentPattern(pattern);
            includePatterns[ascertainmentIncludeCount] = index;
            ascertainmentIncludeCount += 1;
        }
    }

    public int getExcludePatternCount() {
        return ascertainmentExcludeCount;
    }

    public int[] getExcludePatternIndices() {
        return excludePatterns;
    }

    protected void excludePatterns(int excludeFrom, int excludeTo, int every) {
        if (excludePatterns == null)
            excludePatterns = new int[excludeTo - excludeFrom];

        for (int i = excludeFrom; i < excludeTo; i += every) {
            int[] pattern = siteList.getPattern(i);
            int index = addAscertainmentPattern(pattern);
            weights[index] = 0.0; // Site is excluded, so set weight = 0
            excludePatterns[ascertainmentExcludeCount] = index;
            ascertainmentExcludeCount += 1;
        }

    }

    public double getAscertainmentCorrection(double[] patternLogProbs) {
        double excludeProb = 0, includeProb = 0, returnProb = 1.0;

        int[] includeIndices = getIncludePatternIndices();
        int[] excludeIndices = getExcludePatternIndices();
        for (int i = 0; i < getIncludePatternCount(); i++) {
            int index = includeIndices[i];
            includeProb += Math.exp(patternLogProbs[index]);
        }
        for (int j = 0; j < getExcludePatternCount(); j++) {
            int index = excludeIndices[j];
            excludeProb += Math.exp(patternLogProbs[index]);
        }
        if (includeProb == 0.0) {
            returnProb -= excludeProb;
        } else if (excludeProb == 0.0) {
            returnProb = includeProb;
        } else {
            returnProb = includeProb - excludeProb;
        }
        return Math.log(returnProb);
    }

    private int addAscertainmentPattern(int[] pattern) {
        for (int i = 0; i < patternCount; i++) {
            if (comparePatterns(patterns[i], pattern)) {
                return i;
            }
        }
        int index = patternCount;
        patterns[index] = pattern;
        weights[index] = 0.0;  /* do not affect weight */
        patternCount++;

        return index;
    }
}
