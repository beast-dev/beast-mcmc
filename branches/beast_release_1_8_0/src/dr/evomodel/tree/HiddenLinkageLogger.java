package dr.evomodel.tree;

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

/**
 * @author Aaron Darling
 */
public class HiddenLinkageLogger extends MCLogger {

	HiddenLinkageModel hlm;
	public HiddenLinkageLogger(HiddenLinkageModel hlm, LogFormatter formatter, int logEvery) 
	{
		super(formatter, logEvery, false);
		this.hlm = hlm;
	}

    public void startLogging() {
    	String[] labels = new String[1 + hlm.getData().getReadsTaxa().getTaxonCount()];
    	labels[0] = "iter";
    	for(int i=1; i<labels.length; i++){
    		labels[i] = hlm.getData().getReadsTaxa().getTaxonId(i-1);
    	}
    	this.logLabels(labels);
    }
    public void log(long state) {
    	if(state % logEvery != 0)
    		return;
    	String[] values = new String[1 + hlm.getData().getReadsTaxa().getTaxonCount()];
    	values[0] = new Long(state).toString();
    	for(int i=1; i<values.length; i++){
    		values[i] = new Integer(hlm.getLinkageGroupId(hlm.getData().getReadsTaxa().getTaxon(i-1))).toString();
    	}
    	this.logValues(values);
    }
    
}
