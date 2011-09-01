package dr.inference.loggers;


import java.io.PrintWriter;

import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;

/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


public class AlloppDBUGTUNELogger implements Logger {

	private PrintWriter printWriter;
	private AlloppSpeciesNetworkModel apspnet;
	private AlloppSpeciesBindings apspb;
	
	public AlloppDBUGTUNELogger(PrintWriter pw, AlloppSpeciesNetworkModel apspnet,
			                    AlloppSpeciesBindings apspb) {
         printWriter = pw;
         this.apspnet = apspnet;
         this.apspb = apspb;
    }
	
	
	public void startLogging() {
		
	}
	
	
	public void log(int state) {
		if (state % 1000 == 0) {
			printWriter.println("state = " + state);
			printWriter.write(apspnet.mullabTreeAsText());
			int ngt = apspb.numberOfGeneTrees();
			for (int g = 0; g < ngt; g++) {
				printWriter.println("Gene tree " + g);
				printWriter.write(apspb.genetreeAsText(g));
				printWriter.println(apspb.seqassignsAsText(g));
			}
			printWriter.println();
		}
	}
	

	public void stopLogging() {
		printWriter.flush();
	}
	
}
