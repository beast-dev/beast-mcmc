package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ProductStatistic;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Each case will belong to an infectious category (and in future could have a latent category as well) which
 * corresponds to one of a list of probability distributions (most likely gamma or exponential) for the length of the
 * infectious period. The XML rules for the outbreak class ask for at least one ParametricDistributionModel.
 * Assignment of outbreak to distributions should be handled in whatever script or GUI writes the XML.
 *
 * Intended for situations where no data on infection times exists.
 *
 * User: Matthew Hall
 * Date: 20/08/2013
 */

public class JeffreysCategoryOutbreak extends AbstractOutbreak {

    public static final String JEFFREYS_CATEGORY_OUTBREAK = "jeffreysCategoryOutbreak";
    private double[][] distances;

    public JeffreysCategoryOutbreak(String name, Taxa taxa, boolean hasGeography, boolean hasLatentPeriods){
        super(name, taxa, hasLatentPeriods, hasGeography);
        cases = new ArrayList<AbstractCase>();
    }

    public JeffreysCategoryOutbreak(String name, Taxa taxa, boolean hasGeography){
        this(name, taxa, false, hasGeography);
    }

    public JeffreysCategoryOutbreak(String name, Taxa taxa, boolean hasGeography, boolean hasLatentPeriods,
                                    ArrayList<AbstractCase> cases){
        this(name, taxa, hasGeography, hasLatentPeriods);
        this.cases.addAll(cases);
    }

    public JeffreysCategoryOutbreak(Taxa taxa, boolean hasGeography, boolean hasLatentPeriods){
        this(JEFFREYS_CATEGORY_OUTBREAK, taxa, hasGeography, hasLatentPeriods);
    }

    public JeffreysCategoryOutbreak(Taxa taxa, boolean hasGeography, boolean hasLatentPeriods,
                                    ArrayList<AbstractCase> cases){
        this(JEFFREYS_CATEGORY_OUTBREAK, taxa, hasGeography, hasLatentPeriods, cases);
    }

    private void addCase(String caseID, Date examDate, Date cullDate, Parameter coords, Taxa associatedTaxa){
        JeffreysCategoryCase thisCase = new JeffreysCategoryCase(caseID, examDate, cullDate, coords, associatedTaxa);
        cases.add(thisCase);
        addModel(thisCase);
    }


    public double getDistance(AbstractCase a, AbstractCase b) {
        if(distances==null){
            throw new RuntimeException("Distance matrix has not been initialised");
        }
        return distances[getCaseIndex(a)][getCaseIndex(b)];
    }

    private void buildDistanceMatrix(){
        distances = new double[cases.size()][cases.size()];

        for(int i=0; i<cases.size(); i++){
            for(int j=0; j<cases.size(); j++){
                distances[i][j]=SpatialKernel.EuclideanDistance(getCase(i).getCoords(),getCase(j).getCoords());
            }
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void storeState() {
        //nothing to do
    }

    protected void restoreState() {
        //nothing to do
    }

    protected void acceptState() {
        //nothing to do
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected class JeffreysCategoryCase extends AbstractCase{

        public static final String JEFFREYS_CATEGORY_CASE = "jeffreysCategoryCase";
        private final Parameter coords;

        protected JeffreysCategoryCase(String name, String caseID, Date examDate, Date cullDate,
                                       Parameter coords, Taxa associatedTaxa){
            super(name);
            this.caseID = caseID;
            this.examDate = examDate;
            endOfInfectiousDate = cullDate;
            this.associatedTaxa = associatedTaxa;
            this.coords = coords;
        }


        private JeffreysCategoryCase(String caseID, Date examDate, Date cullDate, Parameter coords,
                                     Taxa associatedTaxa){
            this(JEFFREYS_CATEGORY_CASE, caseID, examDate, cullDate, coords, associatedTaxa);
        }

        public Date getLatestPossibleInfectionDate() {
            Double doubleDate = examDate.getTimeValue();
            return Date.createTimeSinceOrigin(doubleDate, Units.Type.DAYS, examDate.getOrigin());
        }

        public boolean culledYet(double time) {
            return time>endOfInfectiousDate.getTimeValue();
        }

        public boolean examinedYet(double time) {
            return time>examDate.getTimeValue();
        }

        public double[] getCoords() {
            return new double[]{coords.getParameterValue(0), coords.getParameterValue(1)};
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) {
            fireModelChanged();
        }

        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireModelChanged();
        }

        protected void storeState() {
            // nothing to do
        }

        protected void restoreState() {
            // nothing to do
        }

        protected void acceptState() {
            //nothing to do
        }

        public Date getExamDate() {
            return examDate;
        }


    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        //for the outbreak

        public static final String HAS_GEOGRAPHY = "hasGeography";
        public static final String HAS_LATENT_PERIODS = "hasLatentPeriods";

        //for the outbreak

        public static final String CASE_ID = "caseID";
        public static final String CULL_DAY = "cullDay";
        public static final String EXAMINATION_DAY = "examinationDay";
        public static final String COORDINATES = "spatialCoordinates";
        public static final String INFECTION_TIME_BRANCH_POSITION = "infectionTimeBranchPosition";
        public static final String INFECTIOUS_TIME_POSITION = "infectiousTimePosition";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final boolean hasGeography = xo.hasAttribute(HAS_GEOGRAPHY) && xo.getBooleanAttribute(HAS_GEOGRAPHY);
            final boolean hasLatentPeriods = xo.hasAttribute(HAS_LATENT_PERIODS)
                    && xo.getBooleanAttribute(HAS_LATENT_PERIODS);
            final Taxa taxa = (Taxa) xo.getChild(Taxa.class);
            JeffreysCategoryOutbreak cases = new JeffreysCategoryOutbreak(null, taxa, hasGeography, hasLatentPeriods);
            for(int i=0; i<xo.getChildCount(); i++){
                Object cxo = xo.getChild(i);
                if(cxo instanceof XMLObject && ((XMLObject)cxo).getName().equals(JeffreysCategoryCase.JEFFREYS_CATEGORY_CASE)){
                    parseCase((XMLObject)cxo, cases);
                }
            }
            cases.buildDistanceMatrix();
            return cases;
        }

        public void parseCase(XMLObject xo, JeffreysCategoryOutbreak outbreak)
                throws XMLParseException {
            String farmID = (String) xo.getAttribute(CASE_ID);
            final Date cullDate = (Date) xo.getElementFirstChild(CULL_DAY);
            final Date examDate = (Date) xo.getElementFirstChild(EXAMINATION_DAY);
            final Parameter coords = xo.hasChildNamed(COORDINATES) ?
                    (Parameter) xo.getElementFirstChild(COORDINATES) : null;
            Taxa taxa = new Taxa();
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getChild(i) instanceof Taxon){
                    taxa.addTaxon((Taxon)xo.getChild(i));
                }
            }
            outbreak.addCase(farmID, examDate, cullDate, coords, taxa);
        }

        public String getParserDescription(){
            return "Parses a set of 'category' farm outbreak and the information that they all share";
        }

        public Class getReturnType(){
            return SimpleOutbreak.class;
        }

        public String getParserName(){
            return JEFFREYS_CATEGORY_OUTBREAK;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] caseRules = {
                new StringAttributeRule(CASE_ID, "The unique identifier for this farm"),
                new ElementRule(CULL_DAY, Date.class, "The date this farm was culled", false),
                new ElementRule(EXAMINATION_DAY, Date.class, "The date this farm was examined", false),
                new ElementRule(Taxon.class, 0, Integer.MAX_VALUE),
                new ElementRule(INFECTION_TIME_BRANCH_POSITION, Parameter.class, "The exact position on the branch" +
                        " along which the infection of this case occurs that it actually does occur"),
                new ElementRule(INFECTIOUS_TIME_POSITION, Parameter.class, "Parameter taking a value between 0 and" +
                        "1, indicating when from infection (0) to first caused infection (or cull if the outbreak" +
                        "causes no infections) (1) the case became infectious", true),
                new ElementRule(COORDINATES, Parameter.class, "The spatial coordinates of this case", true)
        };

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(JeffreysCategoryCase.JEFFREYS_CATEGORY_CASE, caseRules, 1, Integer.MAX_VALUE),
                new ElementRule(Taxa.class),
                AttributeRule.newBooleanRule(HAS_GEOGRAPHY, true),
                AttributeRule.newBooleanRule(HAS_LATENT_PERIODS, true)
        };
    };


}
