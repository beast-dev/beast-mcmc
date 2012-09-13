package dr.app.tracer.application;

import jam.framework.AbstractFrame;
import jam.framework.MenuFactory;

import javax.swing.*;

/**
 * @author rambaut
 *         Date: Feb 24, 2005
 *         Time: 5:12:11 PM
 */
public class AnalysisMenuFactory implements MenuFactory {

    public static final String DEMOGRAPHIC_RECONSTRUCTION = "Demographic Reconstruction...";
    public static final String BAYESIAN_SKYLINE_RECONSTRUCTION = "Bayesian Skyline Reconstruction...";
    public static final String EXTENDED_BAYESIAN_SKYLINE_RECONSTRUCTION = "Extended Bayesian Skyline Reconstruction...";
    public static final String GMRF_SKYRIDE_RECONSTRUCTION = "GMRF Skyride Reconstruction...";
    public static final String LINEAGES_THROUGH_TIME = "Lineages Through Time...";
    public static final String TRAIT_THROUGH_TIME = "Trait Through Time...";

    public static final String CREATE_TEMPORAL_ANALYSIS = "Create Temporal Analysis...";
    public static final String ADD_DEMOGRAPHIC_RECONSTRUCTION = "Add Demographic Reconstruction...";
    public static final String ADD_BAYESIAN_SKYLINE_RECONSTRUCTION = "Add Bayesian Skyline Reconstruction...";
    public static final String ADD_EXTENDED_BAYESIAN_SKYLINE_RECONSTRUCTION = "Add Extended Bayesian Skyline Reconstruction...";
    public static final String ADD_TIME_DENSITY = "Add Time Density...";

    public static final String MODEL_COMPARISON = "Model Comparison...";

    public static final String CONDITIONAL_POST_DIST = "Find Conditional Posterior Distributions...";

    public String getMenuName() {
        return "Analysis";
    }

    public void populateMenu(JMenu menu, AbstractFrame frame) {
        JMenuItem item;

        if (frame instanceof AnalysisMenuHandler) {
            item = new JMenuItem(((AnalysisMenuHandler) frame).getDemographicAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getBayesianSkylineAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getExtendedBayesianSkylineAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getGMRFSkyrideAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getLineagesThroughTimeAction());
            menu.add(item);
//
//            item = new JMenuItem(((AnalysisMenuHandler)frame).getTraitThroughTimeAction());
//            menu.add(item);
//
            menu.addSeparator();

            item = new JMenuItem(((AnalysisMenuHandler) frame).getCreateTemporalAnalysisAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getAddDemographicAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getAddBayesianSkylineAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getAddExtendedBayesianSkylineAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler) frame).getAddTimeDensityAction());
            menu.add(item);

            menu.addSeparator();

            item = new JMenuItem(((AnalysisMenuHandler) frame).getBayesFactorsAction());
            menu.add(item);

            menu.addSeparator();

            item = new JMenuItem(((AnalysisMenuHandler) frame).getConditionalPosteriorDistAction());
            menu.add(item);

        } else {
            item = new JMenuItem(DEMOGRAPHIC_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(BAYESIAN_SKYLINE_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(EXTENDED_BAYESIAN_SKYLINE_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(GMRF_SKYRIDE_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

//		    item = new JMenuItem(LINEAGES_THROUGH_TIME);
//		    item.setEnabled(false);
//		    menu.add(item);
//
//            item = new JMenuItem(TRAIT_THROUGH_TIME);
//            item.setEnabled(false);
//            menu.add(item);
//
            menu.addSeparator();

            item = new JMenuItem(CREATE_TEMPORAL_ANALYSIS);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(ADD_DEMOGRAPHIC_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(ADD_BAYESIAN_SKYLINE_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(ADD_EXTENDED_BAYESIAN_SKYLINE_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(ADD_TIME_DENSITY);
            item.setEnabled(false);
            menu.add(item);

            menu.addSeparator();

            item = new JMenuItem(MODEL_COMPARISON);
            item.setEnabled(false);
            menu.add(item);

            menu.addSeparator();

            item = new JMenuItem(CONDITIONAL_POST_DIST);
            item.setEnabled(false);
            menu.add(item);
        }

    }

    public int getPreferredAlignment() {
        return LEFT;
    }
}
