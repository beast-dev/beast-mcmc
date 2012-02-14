package dr.app.tracer.application;

import javax.swing.*;

/**
 * @author rambaut
 *         Date: Feb 24, 2005
 *         Time: 5:13:13 PM
 */
public interface AnalysisMenuHandler {

	Action getDemographicAction();

	Action getBayesianSkylineAction();

    Action getExtendedBayesianSkylineAction();

    Action getGMRFSkyrideAction();

	Action getLineagesThroughTimeAction();

    Action getTraitThroughTimeAction();

    Action getCreateTemporalAnalysisAction();

	Action getAddDemographicAction();

	Action getAddBayesianSkylineAction();

    Action getAddExtendedBayesianSkylineAction();

	Action getAddTimeDensityAction();

	Action getBayesFactorsAction();

    Action getConditionalPosteriorDistAction();
}
