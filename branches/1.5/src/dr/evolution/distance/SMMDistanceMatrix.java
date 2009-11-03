package dr.evolution.distance;

import dr.evolution.alignment.PatternList;

/**
 * @author Chieh-Hsi Wu
 * Date: 31/07/2009
 * Time: 4:37:38 PM
 * This class is used to calculate the distance between different microsatellite alleles
 */
public class SMMDistanceMatrix extends DistanceMatrix{

   /**
    * constructor taking a pattern source
    * 
    * @param patterns   a pattern of a microsatellite locus
    */
    public SMMDistanceMatrix(PatternList patterns) {
        super(patterns);
    }

    protected double calculatePairwiseDistance(int taxon1, int taxon2) {

        int[] pattern = patterns.getPattern(0);
        int state1 = pattern[taxon1];
        
        int state2 = pattern[taxon2];
        double distance = 0.0;

        if (!dataType.isAmbiguousState(state1) && !dataType.isAmbiguousState(state2))
            distance = Math.abs(state1 - state2);

        return distance;
    }
}
