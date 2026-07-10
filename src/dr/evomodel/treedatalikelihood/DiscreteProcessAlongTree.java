package dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.PatternList;
import dr.evomodel.siteratemodel.SiteRateModel;

public interface DiscreteProcessAlongTree {

    EvolutionaryProcessDelegate getEvolutionaryProcessDelegate();

    SiteRateModel getSiteRateModel();

    PatternList getPatternList();
}
