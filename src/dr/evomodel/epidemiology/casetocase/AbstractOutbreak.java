package dr.evomodel.epidemiology.casetocase;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Abstract class for outbreaks. Implements PatternList for ease of compatability with AbstractTreeLikelihood, but
 * there is one and only one pattern.
 *
 * User: Matthew Hall
 * Date: 14/04/13
 */

public abstract class AbstractOutbreak extends AbstractModel implements PatternList {

    protected GeneralDataType caseDataType;
    protected TaxonList taxa;
    private final String CASE_NAME = "caseID";

    public AbstractOutbreak(String name, Taxa taxa){
        super(name);
        ArrayList<String> caseNames = new ArrayList<String>();
        for(int i=0; i<taxa.getTaxonCount(); i++){
            caseNames.add((String)taxa.getTaxonAttribute(i, CASE_NAME));
        }
        caseDataType = new GeneralDataType(caseNames);
    }

    protected ArrayList<AbstractCase> cases;

    public ArrayList<AbstractCase> getCases(){
        return new ArrayList<AbstractCase>(cases);
    }

    public int getCaseIndex(AbstractCase thisCase){
        return cases.indexOf(thisCase);
    }

    public int size(){
        return cases.size();
    }

    public AbstractCase getCase(int i){
        return cases.get(i);
    }

    public AbstractCase getCase(String name){
        for(AbstractCase thisCase: cases){
            if(thisCase.getName().equals(name)){
                return thisCase;
            }
        }
        return null;
    }

    public TaxonList getTaxa(){
        return taxa;
    }

    public abstract double getLogLikelihood();

    public abstract double probXInfectedByYAtTimeT(AbstractCase X, AbstractCase Y, double T);

    public abstract double logProbXInfectedByYAtTimeT(AbstractCase X, AbstractCase Y, double T);

    public abstract double probXInfectedByYBetweenTandU(AbstractCase X, AbstractCase Y, double T, double U);

    public abstract double logProbXInfectedByYBetweenTandU(AbstractCase X, AbstractCase Y, double T, double U);

    public abstract double probYInfectiousByTimeT(AbstractCase Y, double T);

    public abstract double logProbYInfectiousByTimeT(AbstractCase Y, double T);

    public abstract double probYInfectedAtTimeT(AbstractCase Y, double T);

    public abstract double logProbYInfectedAtTimeT(AbstractCase Y, double T);

    public abstract double probYInfectedBetweenTandU(AbstractCase Y, double T, double U);

    public abstract double logProbYInfectedBetweenTandU(AbstractCase Y, double T, double U);


    //************************************************************************
    // PatternList implementation
    //************************************************************************

    // not considering the possibility that we are simultaneously reconstructing more than one transmission tree!

    public int getPatternCount(){
        return 1;
    }

    public int getStateCount(){
        return size();
    }

    public int getPatternLength(){
        return taxa.getTaxonCount();
    }

    // with an exact correspondence between taxa and states, the following five methods are ill-fitting, but here if
    // needed.
    // @todo if these are never going to be used, get them to throw exceptions

    public int[] getPattern(int patternIndex){
        int[] out = new int[cases.size()];
        for(int i=0; i<cases.size(); i++){
            out[i] = i;
        }
        return out;
    }

    public int getPatternState(int taxonIndex, int patternIndex){
        return taxonIndex;
    }

    public double getPatternWeight(int patternIndex){
        return 1;
    }

    public double[] getPatternWeights(){
        return new double[]{1};
    }

    public double[] getStateFrequencies(){
        double[] out = new double[cases.size()];
        Arrays.fill(out, 1/cases.size());
        return out;
    }

    public DataType getDataType(){
        return caseDataType;
    }

    //************************************************************************
    // TaxonList implementation
    //************************************************************************

    public int getTaxonCount(){
        return taxa.getTaxonCount();
    }

    public Taxon getTaxon(int taxonIndex){
        return taxa.getTaxon(taxonIndex);
    }

    public String getTaxonId(int taxonIndex){
        return taxa.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(String id){
        return taxa.getTaxonIndex(id);
    }

    public int getTaxonIndex(Taxon taxon){
        return taxa.getTaxonIndex(taxon);
    }

    public List<Taxon> asList(){
        return taxa.asList();
    }

    public Object getTaxonAttribute(int taxonIndex, String name){
        return taxa.getTaxonAttribute(taxonIndex, name);
    }




}
