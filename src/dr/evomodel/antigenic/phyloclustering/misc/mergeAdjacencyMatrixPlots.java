package dr.evomodel.antigenic.phyloclustering.misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.ListIterator;

public class mergeAdjacencyMatrixPlots {

    public static void main(String[] args) {
    	
    	int NUM_BURNINs = 250;

    	// 	    	/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/

    	//int numdata = 363;

 		FileReader fileReader2;

 		try {	
 			
 	    	int numdata = 115; // H1N1

 			String input="/Users/charles/Documents/researchData/clustering/forManuscripts-moreReplicates/H1N1/mds0_8/H1N1_mds.clusterLabels.log";
 			String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript7-31-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_8-adjacencyMatrix.txt";
 	    	
 	    	
 		//	String input="/Users/charles/Documents/researchData/clustering/forManuscripts-moreReplicates/H1N1/mds0_7/H1N1_mds.clusterLabels.log";
 		//	String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript7-31-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_7-adjacencyMatrix.txt";
 	    	
 	    	
 		//	String input="/Users/charles/Documents/researchData/clustering/forManuscripts-moreReplicates/H1N1/mds0_6/H1N1_mds.clusterLabels.log";
 		//	String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript7-31-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_6-adjacencyMatrix.txt";

 	    	
 		//	String input="/Users/charles/Documents/researchData/clustering/forManuscripts-moreReplicates/H1N1/mds0_5/H1N1_mds.clusterLabels.log";
 		//	String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript7-31-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_5-adjacencyMatrix.txt";

 		//	String input="/Users/charles/Documents/researchData/clustering/forManuscripts-moreReplicates/H1N1/mds0_4/H1N1_mds.clusterLabels.log";
 		//	String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript7-31-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_4-adjacencyMatrix.txt";

 	    	
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_1/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_1-adjacencyMatrix.txt";
 			
 	    	//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_3/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_3-adjacencyMatrix.txt";

 	    	//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-adjacencyMatrix.txt";	    	

 	    	//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_1-sample2/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_1-sample2-adjacencyMatrix.txt";
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_1-sample3/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_1-sample3-adjacencyMatrix.txt";
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_1-sample4/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_1-sample4-adjacencyMatrix.txt";
 	    	//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_1-sample5/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_1-sample5-adjacencyMatrix.txt";
 	    	//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H1N1/C3b-mds0_1-sample6/H1N1_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H1N1/adjacencyMatrices/H1N1-C3b-mds0_1-sample6-adjacencyMatrix.txt";

 			
 	    //	int numdata = 402; //H3N2
			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_1/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_1-adjacencyMatrix.txt";
 			
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_05/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_05-adjacencyMatrix.txt";
			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_15/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_15-adjacencyMatrix.txt";
			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-adjacencyMatrix.txt";

 		
 			
 			
 	    	//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_05-sample2/H3N2_mds.clusterLabels.log";
 	    	//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_05-sample2-adjacencyMatrix.txt";
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_05-sample3/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_05-sample3-adjacencyMatrix.txt";
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_05-sample4/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_05-sample4-adjacencyMatrix.txt";
 			//String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_05-sample5/H3N2_mds.clusterLabels.log";
 			//String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_05-sample5-adjacencyMatrix.txt";
 		//	String input="/Users/charles/Documents/researchData/clustering/forManuscript/H3N2/C3b-mds0_05-sample6/H3N2_mds.clusterLabels.log";
 		//	String output="/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysisManuscript2-12-2015/H3N2/adjacencyMatrices/H3N2-C3b-mds0_05-sample6-adjacencyMatrix.txt";

 	    	
 			
	
 	    	
 	    	
 	    	fileReader2 = new FileReader(input);
 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter(output)); 


 			
  			 // fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialConditionWithInitialLocationDrift/H3N2ddCRP.log");
 			 //fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/testModel/longRun/H3N2ddCRP.log");
// 			fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/testModel/30Mrun/H3N2ddCRP-1to6M.log");
 		//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/varyParameters/fixK20/H3N2_mds.clusterLabels.log");
// 			  fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/varyParameters-fullData/default/H3N2_mds.clusterLabels.log");
 			
 			 
 			//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec18/H1N1/mds0_15/H1N1_mds.clusterLabels.log");
 			//BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/H1N1/summaryCooccurrence.txt")); 
 	    	//int numdata = 115; // H1N1

 			//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec19/H1N1/mds0_3/H1N1_mds.clusterLabels.log");
 			//BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/H1N1/summaryCooccurrence-mds0_3.txt")); 
 	    	//int numdata = 115; // H1N1
 	    	
 	    	
// 	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec18/Vic/mds0_15-prior2/Vic_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/Vic/summaryCooccurrence.txt")); 
// 	    	int numdata = 179; //Vic

// 	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec19/Vic/mds0_3-prior4/Vic_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/Vic/summaryCooccurrence-mds0_3-prior4.txt")); 
// 	    	int numdata = 179; //Vic

// 			fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec19/Vic/mds0_2-prior4/Vic_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/Vic/summaryCooccurrence-mds0_2-prior4.txt")); 
// 	    	int numdata = 179; //Vic

 			
 			
// 	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec18/Yam/mds0_15-prior2/Yam_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/Yam/summaryCooccurrence.txt")); 
// 	    	int numdata = 174; //Yam

// 	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec19/Yam/mds0_2-prior4/Yam_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/Yam/summaryCooccurrence-mds0_2-prior4.txt")); 
// 	    	int numdata = 174; //Yam


 	    	

// 	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec18/H3N2/mds0_15/H3N2_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/H3N2/summaryCooccurrence_0_15.txt")); 
// 	    	int numdata = 402; //H3N2
	
 			
//	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec18/H3N2/mds0_1/H3N2_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/H3N2/summaryCooccurrence_0_1.txt")); 
// 	    	int numdata = 402; //H3N2

//	    	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/tmpDec18/H3N2/mds0_05/H3N2_mds.clusterLabels.log");
// 			BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/H3N2/summaryCooccurrence_0_05.txt")); 
// 	    	int numdata = 402; //H3N2
			
 			
 	    	
 	    	//h3n2 full data doesn't fix mds
	    	//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/varyParameters-fullData/default/H3N2_mds.clusterLabels.log");
 			//BufferedWriter outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fluTypes/H3N2/old_summaryCooccurrence.txt")); 
 	    	//int numdata = 402; //H3N2
 			
 			 
// 			

 			  
 			  


 			 
		      //print the matrix into a file
		     
		    //outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialConditionWithInitialLocationDrift/summarizeClusters.txt"));
		    //outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/testModel/30Mrun/summarizeClusters1to6M.txt"));
		    //outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fixK20/summaryCooccurence.txt"));
		    //outputWriter = new BufferedWriter(new FileWriter("/Users/charles/Documents/research/antigenic/GenoPheno/driver/clustering/analysis/fullData-default/summaryCooccurence.txt"));

 			
 			 
 			 
		      BufferedReader bReader2 = new BufferedReader( fileReader2); 		      
 		     int[][] coOccur = new int[numdata][numdata];
 		     //assume initializes to 0
 		      
 		     String line2;
 		      
 		      //this routine may give false results if there are extra lines with spaces
 		     
 		     //ignore the first few lines of the file
 		      line2 = bReader2.readLine();
 		      line2 = bReader2.readLine();
 		      line2 = bReader2.readLine();

 		      for(int i=0; i < NUM_BURNINs; i++){
 		    	 line2 = bReader2.readLine();
 		      }

 		      
 		      System.out.println(line2);
 		      
 		      while((line2 = bReader2.readLine())   != null){
 		      System.out.println(line2); 		      
 		      
 		      String datavalue2[] = line2.split("\t");
 		          
 		     LinkedList<Integer>[] clusterLL = new LinkedList[numdata];
 		    for (int i = 0; i < clusterLL.length; i++){
 		        clusterLL[i] = new LinkedList<Integer>();
 		    }
 		     
 		      for (int i = 0; i < numdata; i++) { 
 		    	 //put data into cluster bin linked list	
 		    	  int num =  (int) Double.parseDouble(datavalue2[i+2]);
 		    	 // System.out.println(num);
 		    	 clusterLL[num ].add(i);
 		      }
 		      for (int i = 0; i < numdata; i++) {
 		    	  //System.out.print("Bin=" + i + " ");
 		    	 int clusterSize = clusterLL[i].size();
 		    	 if(clusterSize >0){
 	 		    	int[] clusterArray = new int[clusterSize]; 
 		    		ListIterator<Integer> iter = clusterLL[i].listIterator();
 		    		for(int j=0; j < clusterSize; j++){
 		    			clusterArray[j] = iter.next();
 		    			//System.out.print(clusterArray[j] + " ");
 		    		}
 		    	 //System.out.println("");
                 //for each cluster, increment each pair of points in the coocurrence matrix 		    	
 		    		for(int a=0; a< clusterSize; a++){
 		    		 for(int b=a; b< clusterSize; b++){
 		    			 coOccur[ clusterArray[a]][clusterArray[b]] ++;
 		    		 }
 		    		}
 		    	 }//clusterSize >0
                //          

                   
 		      }
 	 		 
 		      }//while
 		      bReader2.close();
 
 		      //write the lower triangle in the symmetric matrix
 		      for(int i=0; i< numdata; i++){
 		    	  for(int j=(i+1); j< numdata; j++){
 		 		     coOccur[j][i] = coOccur[i][j]; 		    		  
 		    	  }
 		      }

 		      

 		    for(int i=0; i < numdata; i++){
 		    	for(int j=0; j < numdata; j++){
 		    		outputWriter.write(coOccur[i][j] + " ");// Here I know i cant just write x[0] or anything. Do i need to loop in order to write the array?
 		    	}
 	 		    outputWriter.newLine();
 	 		    outputWriter.flush();  
 		    }

 		    outputWriter.close(); 
 		      
 		      
 		} catch (FileNotFoundException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}       
    	

    }
	
}

