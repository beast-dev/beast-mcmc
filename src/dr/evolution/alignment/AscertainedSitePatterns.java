package dr.evolution.alignment;

import dr.evolution.util.TaxonList;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.List;
import java.util.ArrayList;

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

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                CommonCitations.ALEKSEYENKO_2008
        );
        return citations;
    }    

    public int getIncludePatternCount() {
        return ascertainmentIncludeCount;
    }

    public int[] getIncludePatternIndices() {
        return includePatterns;
    }

    protected void includePatterns(int includeFrom, int includeTo, int every) {
        if (includePatterns == null) {
            includePatterns = new int[includeTo - includeFrom + 1];
        }
        for (int i = includeFrom; i <= includeTo; i += every) {
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
            excludePatterns = new int[excludeTo - excludeFrom + 1];

        for (int i = excludeFrom; i <= excludeTo; i += every) {
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
