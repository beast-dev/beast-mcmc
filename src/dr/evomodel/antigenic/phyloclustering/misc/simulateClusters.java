package dr.evomodel.antigenic.phyloclustering.misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.LinkedList;
import java.util.Random;

import dr.evolution.tree.NodeRef;
import dr.evomodel.antigenic.phyloclustering.TreeClusteringVirusesPrior;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class simulateClusters  {
	
    public static final String CLASSNAME = "SimulateClustersAndHI";

    private TreeModel treeModel;
	private int[] clusterLabels;
//	private Parameter clusterLabels;
	private Parameter indicators;
	private MatrixParameter virusLocations;
    int []membershipToClusterLabelIndexes = null;
    private Parameter virusOffsets;
	private int []correspondingTreeIndexForVirus = null; //relates treeModels's indexing system to cluster label's indexing system of viruses. Gets assigned


	private double[][] distCluster;
	
	//maybe run BEAST to load in the phylogeny.. and then introduce the cut
	//do the cutting
	//and then exit.
	
	//cut the tree
	//simulate HI data


	public simulateClusters(TreeModel treeModel_in, Parameter indicators_in ,int nClusters, int seraPerCluster, MatrixParameter vLoc_in,Parameter virusOffsets_in,
			String HIFile, String distClusterFile, int seedNum, int minSize,
	        double meanHomologousTiter,
	        double unitPerTransition,
	        double measurementNoiseVar){	
		
			Random generator = new Random(seedNum);
			//Random generator = new Random();
			
		
			System.out.println("loading the constructor for SimulateCluster and HI");
			TreeModel treeModel = treeModel_in;

			
	        this.treeModel = treeModel_in;
			this.indicators = indicators_in;
			this.virusLocations = vLoc_in;
			this.virusOffsets = virusOffsets_in;


			System.out.println("sera per cluster = " + seraPerCluster);
			System.out.println("nClusters = " + nClusters);
			System.out.println("minimum cluster size = " + minSize);
			System.out.println("meanHomologousTiter = " + meanHomologousTiter);
			System.out.println("unit per transition = " + unitPerTransition);
			System.out.println("measurement noise variance = " + measurementNoiseVar);
			
			System.out.println("seedNum = " + seedNum);
			
	        //initialize excision points
	          //indicators.setDimension(treeModel.getNodeCount());
	          for(int i=0; i < treeModel.getNodeCount(); i++){
	        	  indicators.setParameterValue(i, 0);
	          }


	          

	          //need to know which are the viruses
	          
	      	int numdata = virusLocations.getColumnDimension();
	    	int numNodes = treeModel.getNodeCount();

	      	
//	      	correspondingTreeIndexForVirus = Tree_Clustering_Shared_Routines.setMembershipTreeToVirusIndexes(numdata, virusLocations, numNodes, treeModel);
	        setMembershipToClusterLabelIndexes(); 

   

	        distCluster = new double[nClusters][];
	        for(int i=0; i < nClusters; i++){
	        	distCluster[i] = new double[nClusters];
	        }
	        
	        int[] clusterSize = new int[1];
	        double[][] serumYear = new double[1][1];
	          
	        //this.max_dim = max_dim_in;
	        //this.virusLocations = virusLocations_in;
	        
	        indicators.setParameterValue( treeModel.getRoot().getNumber() ,1); 
	        int numAddedClusters = 1;
	        
	        while(numAddedClusters < nClusters){
			//Clustering routine:
	        int I_selected;
  			do{
  				//I_selected = (int) (Math.floor(Math.random()*numNodes));
  				I_selected = (int) (Math.floor(generator.nextDouble()*numNodes));
  			}while( (int) indicators.getParameterValue(I_selected) == 1 );

	        indicators.setParameterValue(I_selected, 1);
  			numAddedClusters++;


			
	        
	        //if each cluster has > 5 nodes, then this I_selected = 0
	        
	        //want to count the number of nodes in each cluster

	        setClusterLabelsUsingIndicatorsAndCalculateDiff();
	        
	        clusterSize = new int[numAddedClusters];
	        serumYear = new double[numAddedClusters][numdata];
	        for(int i=0; i < numdata; i++){
	        	serumYear[clusterLabels[i]][clusterSize[clusterLabels[i]]] = virusOffsets.getParameterValue(i);
	        	clusterSize[ clusterLabels[i]] ++;
	        	//System.out.println(clusterLabels[i] );
	        }
	        
	        
	        int reachMinClusterSize = 1;
	        for(int i=0; i < numAddedClusters; i++){
	        	if(clusterSize[i] < minSize){
	        		reachMinClusterSize = 0;
	        	}
	        }
	        
	        if(reachMinClusterSize==1){
		        for(int i=0; i < numAddedClusters; i++){
		        	System.out.print("cluster " + i +  " size = " + clusterSize[i] + "\t");
		        	for(int j=0; j < clusterSize[i]; j++){
		        		System.out.print(serumYear[i][j]+",");
		        	}
		        	System.out.println("");
		        }
		        System.out.println("===========================================================");
	        }
	        
	        if(reachMinClusterSize == 0){
	        	indicators.setParameterValue(I_selected, 0);
	        	numAddedClusters--;
	        }
	        
	        
	        }
	        
	        System.out.println("Cluster labels: ");
	        for(int i=0; i < numdata; i++){
	   		   Parameter v = virusLocations.getParameter(i);
	  		   String curName = v.getParameterName();
	        	System.out.println(  curName + "\t" + clusterLabels[i] );
	        }
	        
	        
	        /*
	        try {
	            //Whatever the file path is.
	            File statText3 = new File("/Users/charles/Documents/researchData/clustering/simulations/H3N2_clusterAssignment.txt");
	            FileOutputStream is3 = new FileOutputStream(statText3);
	            OutputStreamWriter osw3 = new OutputStreamWriter(is3);    
	            Writer w3 = new BufferedWriter(osw3);

		        for(int i=0; i < numdata; i++){
			   		   Parameter v = virusLocations.getParameter(i);
			  		   String curName = v.getParameterName();
			        	w3.write(  curName + "\t" + clusterLabels[i] +"\n");
			     }
	            
		         w3.close();
			    } catch (IOException e3) {
			            System.err.println("Problem writing to the file");
			    }
			    */
	        
	        System.out.println("==============================================================");
	        
	        System.out.println("Indicators that are on:");
	        for(int i=0; i < numNodes; i++){
	        	if((int)indicators.getParameterValue(i) ==1){
	        		System.out.println(i);
	        	}
	        }
	        
	        /*
	        try {
	            //Whatever the file path is.
	            File statText2 = new File("/Users/charles/Documents/researchData/clustering/simulations/H3N2_HI_onIndicators.txt");
	            FileOutputStream is2 = new FileOutputStream(statText2);
	            OutputStreamWriter osw2 = new OutputStreamWriter(is2);    
	            Writer w2 = new BufferedWriter(osw2);
	            
		        for(int i=0; i < numNodes; i++){
		        	if((int)indicators.getParameterValue(i) ==1){
		        		w2.write(i+"\n");
		        	}
		        }

		         w2.close();
			        } catch (IOException e2) {
			            System.err.println("Problem writing to the file");
			        }
			   */ 
	        
	        
			System.out.println("===============================================================");

	        try {
	        	System.out.println("distClusterFile = " + distClusterFile);
	            //Whatever the file path is.
	            File statText5 = new File(distClusterFile);
	            FileOutputStream is5 = new FileOutputStream(statText5);
	            OutputStreamWriter osw5 = new OutputStreamWriter(is5);    
	            Writer w5 = new BufferedWriter(osw5);
	   	        
				for(int i=0; i < nClusters; i++){
					for(int j=0; j < nClusters; j++){
						System.out.print(distCluster[i][j] + " ");
						w5.write(distCluster[i][j] + " ");
					}
					System.out.println("");
					w5.write("\n");
				}
		         w5.close();
			        } catch (IOException e5) {
			            System.err.println("Problem writing to the file");
			    
			        
			     }

			
			System.out.println("===============================================================");
			
	        //NEED
	        //a way of telling how many cluster differences between 2 clusters..
	        //this should be like a table.
			
			//When constructing clusters, also construct the distance table.
	        	//just inherits from the parent + 1.
	        //distClutser[i,j] = #
	        

	        
	        
	        
	        try {
	        	System.out.println("HIFile = " + HIFile);
	            //Whatever the file path is.
	            File statText = new File(HIFile);
	            FileOutputStream is = new FileOutputStream(statText);
	            OutputStreamWriter osw = new OutputStreamWriter(is);    
	            Writer w = new BufferedWriter(osw);
	            w.write("virusIsolate\tvirusStrain\tvirusYear\tserumIsolate\tserumStrain\tserumYear\ttiter\tsource\n");
	   
	        
	        System.out.println("virusIsolate    virusStrain     virusYear       serumIsolate    serumStrain     serumYear       titer   source");
	        for(int c=0; c< nClusters; c++){
	        	for(int s=0; s<seraPerCluster; s++ ){
	        		
	        		//sample serum year
	        		//int whichSample = (int) (Math.floor(Math.random()*clusterSize[c]));
	        		int whichSample = (int) (Math.floor(generator.nextDouble()*clusterSize[c]));
	        		double curSerumYear = serumYear[c][whichSample];
	        		
	        		for(int v=0; v< numdata; v++){
		        		double[] meanTiter = new double[1];
		        		meanTiter[0] = meanHomologousTiter - distCluster[c][clusterLabels[v]]*unitPerTransition;
		        		double[][] prec = new double[1][1];
		        		prec[0][0] = 1/measurementNoiseVar;
		        
		        		double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(meanTiter, prec);
		        		
		        		Parameter virus = virusLocations.getParameter(v);
		        		String curName = virus.getParameterName();
		        		//System.out.print("C" + clusterLabels[v] + "\t" + curName +"\t" + virusOffsets.getParameterValue(v) +"\t"  );
		        		//System.out.print("C" + c +"\t" + "c" + c + "s" + (s+1) +"\t" + curSerumYear +"\t");
		        		//System.out.println(Math.pow(2,values[0]) + "\tsimulation");
	
	//if(Math.random() < 0.1){
		        		//if(distCluster[c][clusterLabels[v]] <= 3){
		        		w.write("C" + clusterLabels[v] + "\t" + curName +"\t" + virusOffsets.getParameterValue(v) +"\t"  );
		        		w.write("c" + c + "s" + (s+1) +"\t"+ "c" + c + "s" + (s+1) +"\t" + curSerumYear +"\t");
		        		w.write(Math.pow(2,values[0]) + "\tsimulation\n");
	//	}       		
	        		}
	        	}
	        }

	         w.close();
		        } catch (IOException e) {
		            System.err.println("Problem writing to the file");
		    
		        
		     }
		    
	        
	        
	    	
	    	//Make XML
	    	//Run
	    	
			//System.exit(0);
	}

	
	
	

	private void setMembershipToClusterLabelIndexes(){

  	   //I suspect this is an expensive operation, so I don't want to do it many times,
  	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
		int numdata = virusLocations.getColumnDimension();
		int numNodes = treeModel.getNodeCount();
         membershipToClusterLabelIndexes = new int[numdata]; 
         clusterLabels = new int[numdata];
         for(int i=0; i < numdata; i++){
  		   Parameter v = virusLocations.getParameter(i);
  		   String curName = v.getParameterName();
  		  // System.out.println(curName);
  		   int isFound = 0;
      	   for(int j=0; j < numNodes; j++){
      		   String treeId = treeModel.getTaxonId(j);
      		   if(curName.equals(treeId) ){
      		//	   System.out.println("  isFound at j=" + j);
      			   membershipToClusterLabelIndexes[i] = j;
      			   isFound=1;
      			   break;
      		   }	   
      	   }
      	   if(isFound ==0){
      		   System.out.println("not found. Exit now.");
      		   System.exit(0);
      	   }     	   
         }
    }
	
	
	private void setClusterLabelsUsingIndicatorsAndCalculateDiff(){

        int []membership = determine_membership_v3(treeModel);
		int numdata = virusLocations.getColumnDimension();
        for(int i=0; i < numdata; i++){   
        	clusterLabels[i] = membership[membershipToClusterLabelIndexes[i]] ;
        }
	}
	
	

    //traverse down the tree, top down, do calculation
     int[] determine_membership_v3(TreeModel treeModel){
	    	
	    NodeRef root = treeModel.getRoot();
	
	    int numClusters = 1;
	    LinkedList<NodeRef> list = new LinkedList<NodeRef>();
	    list.addFirst(root);
	
	    int[] membership = new int[treeModel.getNodeCount()];
	    for(int i=0; i < treeModel.getNodeCount(); i++){
	    	membership[i] = -1;
	    }
	    membership[root.getNumber()] = 0; //root always given the first cluster
	          
	    while(!list.isEmpty()){
	    	//do things with the current object
	    	NodeRef curElement = list.pop();
	    	//String content = "node #" + curElement.getNumber() +", taxon=" + treeModel.getNodeTaxon(curElement) + " and parent is = " ;
	    	String content = "node #" + curElement.getNumber() +", taxon= " ;
	    	if(treeModel.getNodeTaxon(curElement)== null){
	    		content += "internal node\t";
	    	}
	    	else{
	    		content += treeModel.getNodeTaxon(curElement).getId() + "\t";
	    		//content += treeModel.getTaxonIndex(treeModel.getNodeTaxon(curElement)) + "\t";
	    	}
	    	
	       	if(treeModel.getParent(curElement)== null){
	    		//content += "no parent";
	    	}
	    	else{
	    		//content += "parent node#=" + treeModel.getParent(curElement).getNumber();
	    	}
	    	
	    	//cluster assignment:
	    	if(!treeModel.isRoot(curElement)){
	    		if( (int) indicators.getParameterValue(curElement.getNumber() ) == 1) {
	    			numClusters++ ;
	    			membership[ curElement.getNumber() ] = numClusters - 1; 
	    			
	    			int parentClusterLabel = membership[treeModel.getParent(curElement).getNumber()] ; //parent cluster label
	    			//assign distCluster
	    			for(int i=0; i < (numClusters-1); i++){
	    				distCluster[numClusters -1][i] = distCluster[ parentClusterLabel][i] +1;
	    				distCluster[i][numClusters -1] = distCluster[i][ parentClusterLabel] +1; 

	    			}
	      	 	}
	    		else{
	    			//inherit from parent's cluster assignment
	    			membership[curElement.getNumber()] = membership[treeModel.getParent(curElement).getNumber()]; 
	    		}        	
	    	}//is not Root
	    	content += " cluster = " + membership[curElement.getNumber()] ; 
	    	
	    //	System.out.println(content);
	
	    	
	        for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
	        	list.addFirst(treeModel.getChild(curElement,childNum));
	        }
	    }
	
	     return(membership);
    }

	



    
    
    //traverse down the tree, top down, do calculation
    int[] determine_from_membership_v2(TreeModel treeModel){
	    	//note: I set MAX_DIM as the most I would print, but in order to avoid bug, I 
    	//declare the number of nodes as the most active nodes I can have.
    	int[] fromMembership = new int[treeModel.getNodeCount()];
    	for(int i=0; i < treeModel.getNodeCount(); i++){
    		fromMembership[i ] = -1;
    	}
    	
	    NodeRef root = treeModel.getRoot();
	
	    int numClusters = 1;
	    LinkedList<NodeRef> list = new LinkedList<NodeRef>();
	    list.addFirst(root);
	
	    int[] membership = new int[treeModel.getNodeCount()];
	    for(int i=0; i < treeModel.getNodeCount(); i++){
	    	membership[i] = -1;
	    }
	    membership[root.getNumber()] = 0; //root always given the first cluster
	          
	    while(!list.isEmpty()){
	    	//do things with the current object
	    	NodeRef curElement = list.pop();
	    	//String content = "node #" + curElement.getNumber() +", taxon=" + treeModel.getNodeTaxon(curElement) + " and parent is = " ;
	    	String content = "node #" + curElement.getNumber() +", taxon= " ;
	    	if(treeModel.getNodeTaxon(curElement)== null){
	    		content += "internal node\t";
	    	}
	    	else{
	    		content += treeModel.getNodeTaxon(curElement).getId() + "\t";
	    		//content += treeModel.getTaxonIndex(treeModel.getNodeTaxon(curElement)) + "\t";
	    	}
	    	
	       	if(treeModel.getParent(curElement)== null){
	    		//content += "no parent";
	    	}
	    	else{
	    		//content += "parent node#=" + treeModel.getParent(curElement).getNumber();
	    	}
	    	
	    	//cluster assignment:
	    	if(!treeModel.isRoot(curElement)){
	    	 if( (int) indicators.getParameterValue(curElement.getNumber() ) == 1) {
	    //		 System.out.print("indicator # " + curElement.getNumber()  + " ");
	    		numClusters++ ;
	    		membership[ curElement.getNumber() ] = numClusters - 1; 
	    		fromMembership[numClusters -1] = membership[ treeModel.getParent(curElement).getNumber()];
	    //		System.out.println("    membership " + (numClusters-1) + " assigned from " + membership[ treeModel.getParent(curElement).getNumber()] );
	      	}
	    	else{
	    		//inherit from parent's cluster assignment
	    		membership[curElement.getNumber()] = membership[treeModel.getParent(curElement).getNumber()]; 
	    	 }
	    	        	
	    	}//is not Root
	    	content += " cluster = " + membership[curElement.getNumber()] ; 
	    	
	    //	System.out.println(content);
	
	    	
	        for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
	        	list.addFirst(treeModel.getChild(curElement,childNum));
	        }
	    }
	
	     return(fromMembership);
   }

	
    
    
	
	
	
    

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	   	
    	public final static String EXCISIONPOINTS = "excisionPoints";
    	public final static String CLUSTERLABELS = "clusterLabels";
    	public final static String CLUSTERLABELSTREENODE = "clusterLabelsTreeNode";

    	public final static String  MU = "mu";

    	public final static String OFFSETS = "offsets";
    	public final static String VIRUS_LOCATIONS = "virusLocations";
    	public final static String VIRUS_LOCATIONSTREENODE = "virusLocationsTreeNodes";
    	
    	public final static String INDICATORS = "indicators";

    	public final static String VIRUS_OFFSETS = "virusOffsets";
        boolean integrate = false;
        
        
     //   public final static String MUVARIANCE = "muVariance";
        public final static String MUPRECISION = "muPrecision";
        public final static String PROBACTIVENODE = "probActiveNode";
        
        public final static String INITIALNUMCLUSTERS = "numClusters";
        public final static String NUMSERA = "numSeraPerCluster";

        public final static String MUMEAN = "muMean";

        public final static String FILE_NAME = "HIFile";
        public final static String FILE_NAME2 = "distClusterFile";

        public final static String SEEDNUM = "seed";
        public final static String MINSIZECLUSTER = "minClusterSize";
        
        
        public final static String MEANHOMOLOGOUSTITER = "meanHomologousLog2Titer";
        public final static String UNITPERTRANSITION = "log2TiterDiffPerTransition";
        public final static String MEASUREMENTNOISEVAR = "measurementNoiseVariance";
        
        public String getParserName() {
            return CLASSNAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        	
        	
        		int initialK = 10;
            	if (xo.hasAttribute(INITIALNUMCLUSTERS)) {
            		initialK = xo.getIntegerAttribute(INITIALNUMCLUSTERS);
            	}

        		int numSera = 10;
            	if (xo.hasAttribute(NUMSERA)) {
            		numSera = xo.getIntegerAttribute(NUMSERA);
            	}
            	
        		int seedNum = (int) Math.floor(Math.random()*100000);
            	if (xo.hasAttribute(SEEDNUM)) {
            		seedNum = xo.getIntegerAttribute(SEEDNUM);
            	}
            	
            	int minSizeCluster = 5;
            	if (xo.hasAttribute(MINSIZECLUSTER)) {
            		minSizeCluster = xo.getIntegerAttribute(MINSIZECLUSTER);
            	}
            	

    	        double meanHomologousTiter = 10; //log2
    	        if(xo.hasAttribute(MEANHOMOLOGOUSTITER)){
    	        	meanHomologousTiter = xo.getDoubleAttribute(MEANHOMOLOGOUSTITER);
    	        }
    	        double unitPerTransition = 2;
    	        if(xo.hasAttribute(UNITPERTRANSITION)){
    	        	unitPerTransition = xo.getDoubleAttribute(UNITPERTRANSITION);
    	        }
    	        double measurementNoiseVar = 1;
            	if(xo.hasAttribute(MEASUREMENTNOISEVAR)){
            		measurementNoiseVar = xo.getDoubleAttribute(MEASUREMENTNOISEVAR);
            	}
            	
                TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
                XMLObject cxo  = xo.getChild(INDICATORS);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);
                cxo=xo.getChild(VIRUS_LOCATIONS);
                MatrixParameter virusLocations =(MatrixParameter) cxo.getChild(MatrixParameter.class);                
               
                cxo=xo.getChild(VIRUS_OFFSETS);
                Parameter virusOffsets =(Parameter) cxo.getChild(Parameter.class);      
                
                String fileName = xo.getStringAttribute(FILE_NAME);
                String fileName2 = xo.getStringAttribute(FILE_NAME2);
                
                

		        return new simulateClusters(treeModel, indicators, initialK, numSera, virusLocations, virusOffsets, fileName, fileName2, seedNum, minSizeCluster,
		        		meanHomologousTiter, unitPerTransition, measurementNoiseVar); 
            }
        
        	private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newStringRule(FILE_NAME2, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(INITIALNUMCLUSTERS, true, "the initial number of clusters"),
                AttributeRule.newIntegerRule(SEEDNUM, true, "the initial number of clusters"),
                AttributeRule.newIntegerRule(NUMSERA, true, "number of sera to simulate"),
                AttributeRule.newIntegerRule(MINSIZECLUSTER, true, "The minimum cluster size of a virus"),
                AttributeRule.newDoubleRule(MEANHOMOLOGOUSTITER, true, "the expected log2 titer of a homologous virus"),
                AttributeRule.newDoubleRule(UNITPERTRANSITION, true, "the expected decrease in log2 titer value per major transition"),
                AttributeRule.newDoubleRule(MEASUREMENTNOISEVAR, true, "the variance of the measurement noise in the log2 titer value"),
                new ElementRule(TreeModel.class),
                new ElementRule(INDICATORS, Parameter.class),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class),
                new ElementRule(VIRUS_OFFSETS, Parameter.class),                     
        };

            //************************************************************************
            // AbstractXMLObjectParser implementation
            //************************************************************************

            public String getParserDescription() {
                return "tree clustering viruses";
            }

            public Class getReturnType() {
                return TreeClusteringVirusesPrior.class;
            }

            public XMLSyntaxRule[] getSyntaxRules() {
                return rules;
            }
            
            

            
    };

    String Atribute = null;

	public Model getModel() {
		// TODO Auto-generated method stub
		return null;
	}





	public double getLogLikelihood() {
		// TODO Auto-generated method stub
		return 0;
	}





	public void makeDirty() {
		// TODO Auto-generated method stub
		
	}




	
}
