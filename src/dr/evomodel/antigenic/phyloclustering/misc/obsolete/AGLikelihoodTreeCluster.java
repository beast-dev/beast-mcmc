package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.LogTricks;
import dr.math.distributions.NormalDistribution;
import dr.util.*;
import dr.xml.*;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Charles Cheung
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
/*
    Both virus locations and serum locations are shifted by the parameter locationDrift.
    A location is increased by locationDrift x offset.
    Offset is set to 0 for the earliest virus and increasing with difference in date from earliest virus.
    This makes the raw virusLocations and serumLocations parameters not directly interpretable.
*/
public class AGLikelihoodTreeCluster extends AbstractModelLikelihood implements Citable {
    private static final boolean CHECK_INFINITE = false;
    private static final boolean USE_THRESHOLDS = true;
    private static final boolean USE_INTERVALS = true;

    public final static String AG_LIKELIHOOD = "aglikelihoodtreecluster";

    // column indices in table
    private static final int VIRUS_ISOLATE = 0;
    private static final int VIRUS_STRAIN = 1;
    private static final int VIRUS_DATE = 2;
    private static final int SERUM_ISOLATE = 3;
    private static final int SERUM_STRAIN = 4;
    private static final int SERUM_DATE = 5;
    private static final int TITRE = 6;
    

    private  double oldLogLikelihood =0;

    public enum MeasurementType {
        INTERVAL,
        POINT,
        THRESHOLD,
        MISSING
    }

    public AGLikelihoodTreeCluster(
            int mdsDimension,
            Parameter mdsPrecisionParameter,
            Parameter locationDriftParameter,
            Parameter virusDriftParameter,
            Parameter serumDriftParameter,
            MatrixParameter virusLocationsParameter,
            MatrixParameter serumLocationsParameter,
            CompoundParameter tipTraitsParameter,
            Parameter virusOffsetsParameter,
            Parameter serumOffsetsParameter,
            Parameter serumPotenciesParameter,
            Parameter serumBreadthsParameter,
            Parameter virusAviditiesParameter,
            DataTable<String[]> dataTable,
            boolean mergeSerumIsolates,
            double intervalWidth,
            double driftInitialLocations, 
            boolean clusterMeans,
            Parameter clusterOffsetsParameter) {

        super(AG_LIKELIHOOD);

        this.intervalWidth = intervalWidth;
        boolean useIntervals = USE_INTERVALS && intervalWidth > 0.0;

        int thresholdCount = 0;
             
        double earliestDate = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dataTable.getRowCount(); i++) {

            String[] values = dataTable.getRow(i);

            String virusName = values[VIRUS_STRAIN];
            double virusDate = Double.parseDouble(values[VIRUS_DATE]);
            int virus = virusNames.indexOf(virusName);
            if (virus == -1) {
                virusNames.add(virusName);
                virusDates.add(virusDate);
                virus = virusNames.size() - 1;
            }

            String serumName = "";
            if (mergeSerumIsolates) {
                serumName = values[SERUM_STRAIN];
            } else {
                serumName = values[SERUM_ISOLATE];
            }
            double serumDate = Double.parseDouble(values[SERUM_DATE]);
            int serum = serumNames.indexOf(serumName);
            if (serum == -1) {
                serumNames.add(serumName);
                serumDates.add(serumDate);
                serum = serumNames.size() - 1;
            }
            
            boolean isThreshold = false;
            boolean isLowerThreshold = false;
            double rawTitre = Double.NaN;
            if (values[TITRE].length() > 0) {
                try {
                    rawTitre = Double.parseDouble(values[TITRE]);
                } catch (NumberFormatException nfe) {
                    // check if threshold below
                    if (values[TITRE].contains("<")) {
                        rawTitre = Double.parseDouble(values[TITRE].replace("<",""));
                        isThreshold = true;
                        isLowerThreshold = true;
                        thresholdCount++;
                    }
                    // check if threshold above
                    if (values[TITRE].contains(">")) {                	
                        rawTitre = Double.parseDouble(values[TITRE].replace(">",""));
                        isThreshold = true;
                        isLowerThreshold = false;
                        thresholdCount++;                    	
                        //throw new IllegalArgumentException("Error in measurement: unsupported greater than threshold at row " + (i+1));
                    }
                }
            }

            if (serumDate < earliestDate) {
                earliestDate = serumDate;
            }

            if (virusDate < earliestDate) {
                earliestDate = virusDate;
            }

            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(virus, serum, virusDate, serumDate, type, rawTitre, isLowerThreshold);

            if (USE_THRESHOLDS || !isThreshold) {
                measurements.add(measurement);
            }

        }

        double[] maxColumnTitres = new double[serumNames.size()];
        for (Measurement measurement : measurements) {
            double titre = measurement.log2Titre;
            if (Double.isNaN(titre)) {
                titre = measurement.log2Titre;
            }
            if (titre > maxColumnTitres[measurement.serum]) {
                maxColumnTitres[measurement.serum] = titre;
            }
        }

        this.mdsDimension = mdsDimension;

        this.mdsPrecisionParameter = mdsPrecisionParameter;
        addVariable(mdsPrecisionParameter);

        this.locationDriftParameter = locationDriftParameter;
        if (this.locationDriftParameter != null) {
            addVariable(locationDriftParameter);
        }

        this.virusDriftParameter = virusDriftParameter;
        if (this.virusDriftParameter != null) {
            addVariable(virusDriftParameter);
        }

        this.serumDriftParameter = serumDriftParameter;
        if (this.serumDriftParameter != null) {
            addVariable(serumDriftParameter);
        }

        this.virusLocationsParameter = virusLocationsParameter;
        if (this.virusLocationsParameter != null) {
            setupLocationsParameter(virusLocationsParameter, virusNames);
        }

        this.serumLocationsParameter = serumLocationsParameter;
        if (this.serumLocationsParameter != null) {
            setupLocationsParameter(serumLocationsParameter, serumNames);
        }

        this.tipTraitsParameter = tipTraitsParameter;
        if (tipTraitsParameter != null) {
            setupTipTraitsParameter(this.tipTraitsParameter, virusNames);
        }

        this.virusOffsetsParameter = virusOffsetsParameter;
        if (virusOffsetsParameter != null) {
            setupOffsetsParameter(virusOffsetsParameter, virusNames, virusDates, earliestDate);
        }

        this.serumOffsetsParameter = serumOffsetsParameter;
        if (serumOffsetsParameter != null) {
            setupOffsetsParameter(serumOffsetsParameter, serumNames, serumDates, earliestDate);
        }

        this.serumPotenciesParameter = setupSerumPotencies(serumPotenciesParameter, maxColumnTitres);
        this.serumBreadthsParameter = setupSerumBreadths(serumBreadthsParameter);
        this.virusAviditiesParameter = setupVirusAvidities(virusAviditiesParameter);

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicLikelihood:\n");
        sb.append("\t\t" + virusNames.size() + " viruses\n");
        sb.append("\t\t" + serumNames.size() + " sera\n");
        sb.append("\t\t" + measurements.size() + " assay measurements\n");
        if (USE_THRESHOLDS) {
            sb.append("\t\t" + thresholdCount + " thresholded measurements\n");
        }
        if (useIntervals) {
            sb.append("\n\t\tAssuming a log 2 measurement interval width of " + intervalWidth + "\n");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());

        virusLocationChanged = new boolean[this.virusLocationsParameter.getParameterCount()];
        serumLocationChanged = new boolean[this.serumLocationsParameter.getParameterCount()];
        virusEffectChanged = new boolean[virusNames.size()];
        serumEffectChanged = new boolean[serumNames.size()];
        logLikelihoods = new double[measurements.size()];
        storedLogLikelihoods = new double[measurements.size()];

       // driftInitialLocations = 1; //charles added - now specified in the xml
   //     setupInitialLocations(driftInitialLocations);
   //     loadInitialLocations(virusNames, serumNames);
        
        //System.out.println("Print now!");
		//      for (int i = 0; i < virusLocationsParameter.getParameterCount(); i++) {    	  
		 //   	 System.out.print(virusLocationsParameter.getParameter(i).getParameterValue(0) + " ");
		  //  	 System.out.print(virusLocationsParameter.getParameter(i).getParameterValue(1) + " ");  	  
		   //   }
		   //   System.out.println("");
     

		        if(clusterMeans){
		        	this.clusterMeans = clusterMeans;
		        	this.clusterOffsetsParameter = clusterOffsetsParameter;
		        	
		        	
		        	//if(clusterOffsetsParameter != null){
		        	//System.out.println("virusNames.size()="+ virusNames.size());
		        	//clusterOffsetsParameter.setDimension( virusNames.size());  
		        //    for (int i = 0; i < virusNames.size(); i++) {
		           // 	clusterOffsetsParameter.setId(virusNames.get(i));
		           // }
		            //addVariable(clusterOffsetsParameter);
		        	//}
		        	
		        	//stay null
		           if (clusterOffsetsParameter == null) {
		            //	clusterOffsetsParameter = new Parameter.Default("clusterOffsets");
		            } else {
		            	//clusterOffsetsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1000));
		                addVariable(clusterOffsetsParameter);
			            clusterOffsetsParameter.setDimension(virusNames.size());

		            }
		        	
		        	
		        	System.out.println(" clusterMeans = true");
		        	//System.exit(0);
		        }

		      
        makeDirty();
        
        loadInitialSerumLocations(serumNames);
        
        
        System.out.println("======================================");
    	double sum=0;
    	double s2 = 0;
		for(int i=0; i <  serumLocationsParameter.getParameterCount(); i++){
			sum+=serumLocationsParameter.getParameter(i).getParameterValue(0);
			s2+=serumLocationsParameter.getParameter(i).getParameterValue(1);
		}
		System.out.println("sum sera location dimension 1 = " + sum);
		System.out.println("sum sera location dimension 2 = " + s2);

        
        
        System.out.println("======================================");
       // System.exit(0);
        
        
    }

    
    
    //load initial serum location - load the last line
    private void loadInitialSerumLocations(List<String> serumNames) {

		FileReader fileReader2;
		try {
//			fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test23/run4/H3N2_mds.serumLocs.log");
			//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run64/H3N2_mds.serumLocs.log");
			fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run79/H3N2_mds.serumLocs.log");
		     /**
		       * Creating a buffered reader to read the file
		       */
		      BufferedReader bReader2 = new BufferedReader( fileReader2);

		      String line;
		      line = bReader2.readLine();
		      line = bReader2.readLine();
		      
		      line = bReader2.readLine();
		      System.out.println(line);
		      String namevalue[] = line.split("\t");

		      
		      //skip to the last line
		      String testLine;
		      while ((testLine = bReader2.readLine()) != null){
		    	  line = testLine;
		      }

		      System.out.println(line);
		      
		      String datavalue[] = line.split("\t");

		   //   double sumDim2=0;
		       //   System.out.println(serumLocationsParameter.getParameterCount());
		      for (int i = 0; i < serumLocationsParameter.getParameterCount(); i++) {
		    	  //int index = findStrain( namevalue[i*2+1], serumNames);  //don't enable this.. this will cause a bug because the serumNames are not unique.
		    	  //System.out.println("index=" + index);
		    	  
		    	  double dim1 = Double.parseDouble(datavalue[i*2+1])- serumDriftParameter.getParameterValue(0)*serumOffsetsParameter.getParameterValue(i);
		    	  double dim2 = Double.parseDouble(datavalue[i*2+2]);
		    	 // System.out.println(datavalue[i*2+1]);
		    	  serumLocationsParameter.getParameter(i).setParameterValue(0, dim1);
		    	  serumLocationsParameter.getParameter(i).setParameterValue(1, dim2);
		    //	  sumDim2+= dim2;
		    //	  System.out.print(dim2 + "\t");
		          //virusLocationsParameter.getParameter(i).setParameterValue(0, 1);
		   	  
		      }
		  //    System.out.println("\n sum dim2 = " + sumDim2);
//	    	  System.exit(0);

		      bReader2.close();
		
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}        
      
    	

    }
        
    
    
    private Parameter setupVirusAvidities(Parameter virusAviditiesParameter) {
        // If no row parameter is given, then we will only use the serum effects
        if (virusAviditiesParameter != null) {
            virusAviditiesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, Double.MIN_VALUE, 1));
            virusAviditiesParameter.setDimension(virusNames.size());
            addVariable(virusAviditiesParameter);
            String[] labelArray = new String[virusNames.size()];
            virusNames.toArray(labelArray);
            virusAviditiesParameter.setDimensionNames(labelArray);
            for (int i = 0; i < virusNames.size(); i++) {
                virusAviditiesParameter.setParameterValueQuietly(i, 0.0);
            }
        }
        return virusAviditiesParameter;
    }

    private Parameter setupSerumPotencies(Parameter serumPotenciesParameter, double[] maxColumnTitres) {
        // If no serum potencies parameter is given, make one to hold maximum values for scaling titres...
        if (serumPotenciesParameter == null) {
            serumPotenciesParameter = new Parameter.Default("serumPotencies");
        } else {
            serumPotenciesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            addVariable(serumPotenciesParameter);
        }

        serumPotenciesParameter.setDimension(serumNames.size());
        String[] labelArray = new String[serumNames.size()];
        serumNames.toArray(labelArray);
        serumPotenciesParameter.setDimensionNames(labelArray);
        for (int i = 0; i < maxColumnTitres.length; i++) {
            serumPotenciesParameter.setParameterValueQuietly(i, maxColumnTitres[i]);
        }

        return serumPotenciesParameter;
    }

    private Parameter setupSerumBreadths(Parameter serumBreadthsParameter) {
        // If no serum breadths parameter is given, then we will only use the serum potencies
        if (serumBreadthsParameter != null) {
            serumBreadthsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            serumBreadthsParameter.setDimension(serumNames.size());
            addVariable(serumBreadthsParameter);
            String[] labelArray = new String[serumNames.size()];
            serumNames.toArray(labelArray);
            serumBreadthsParameter.setDimensionNames(labelArray);
            for (int i = 0; i < serumNames.size(); i++) {
                serumBreadthsParameter.setParameterValueQuietly(i, 1.0);
            }
        }
        return serumBreadthsParameter;
    }

    protected void setupLocationsParameter(MatrixParameter locationsParameter, List<String> strains) {
        locationsParameter.setColumnDimension(mdsDimension);
        locationsParameter.setRowDimension(strains.size());
        for (int i = 0; i < strains.size(); i++) {
            locationsParameter.getParameter(i).setId(strains.get(i));
        }
        addVariable(locationsParameter);
    }

    private void setupOffsetsParameter(Parameter offsetsParameter, List<String> strainNames, List<Double> strainDates, double earliest) {

        offsetsParameter.setDimension(strainNames.size());
        String[] labelArray = new String[strainNames.size()];
        strainNames.toArray(labelArray);
        offsetsParameter.setDimensionNames(labelArray);
        for (int i = 0; i < strainNames.size(); i++) {
            Double offset = strainDates.get(i) - new Double(earliest);
            if (offset == null) {
                throw new IllegalArgumentException("Date missing for strain: " + strainNames.get(i));
            }
            offsetsParameter.setParameterValue(i, offset);
        }
        addVariable(offsetsParameter);
    }


    private void setupTipTraitsParameter(CompoundParameter tipTraitsParameter, List<String> strainNames) {
        tipIndices = new int[strainNames.size()];

        for (int i = 0; i < strainNames.size(); i++) {
            tipIndices[i] = -1;
        }

        for (int i = 0; i < tipTraitsParameter.getParameterCount(); i++) {
            Parameter tip = tipTraitsParameter.getParameter(i);
            String label = tip.getParameterName();
            int index = findStrain(label, strainNames);
            if (index != -1) {
                if (tipIndices[index] != -1) {
                    throw new IllegalArgumentException("Duplicated tip name: " + label);
                }

                tipIndices[index] = i;

                // rather than setting these here, we set them when the locations are set so the changes propagate
                // through to the diffusion model.
//                Parameter location = locationsParameter.getParameter(index);
//                for (int dim = 0; dim < mdsDimension; dim++) {
//                    tip.setParameterValue(dim, location.getParameterValue(dim));
//                }
            } else {
                // The tree may contain viruses not present in the assay data
                //       throw new IllegalArgumentException("Unmatched tip name in assay data: " + label);
            }
        }
        // we are only setting this parameter not listening to it:
//        addVariable(this.tipTraitsParameter);
    }

    private final int findStrain(String label, List<String> strainNames) {
        int index = 0;
        for (String strainName : strainNames) {
            if (label.startsWith(strainName)) {
                return index;
            }

            index ++;
        }
        return -1;
    }

    private void setupInitialLocations(double drift) {
    	//System.out.println("hihi");
        for (int i = 0; i < virusLocationsParameter.getParameterCount(); i++) {
            double offset = 0.0;
            if (virusOffsetsParameter != null) {
            	//System.out.print("virus Offset Parameter present"+ ": ");
            	//System.out.print( virusOffsetsParameter.getParameterValue(i) + " ");
            	//System.out.print(" drift= " + drift + " ");
                offset = drift * virusOffsetsParameter.getParameterValue(i);
            }
            else{
            	System.out.println("virus Offeset Parameter NOT present");
            }
            double r = MathUtils.nextGaussian() + offset;
            virusLocationsParameter.getParameter(i).setParameterValue(0, r);
           // System.out.println (  virusLocationsParameter.getParameter(i).getParameterValue(0));
            if (mdsDimension > 1) {
                for (int j = 1; j < mdsDimension; j++) {
                    r = MathUtils.nextGaussian();
                    virusLocationsParameter.getParameter(i).setParameterValue(j, r);
                }
            }
        }
        for (int i = 0; i < serumLocationsParameter.getParameterCount(); i++) {
            double offset = 0.0;
            if (serumOffsetsParameter != null) {
                offset = drift * serumOffsetsParameter.getParameterValue(i);
            }
            double r = MathUtils.nextGaussian() + offset;
            serumLocationsParameter.getParameter(i).setParameterValue(0, r);
            if (mdsDimension > 1) {
                for (int j = 1; j < mdsDimension; j++) {
                    r = MathUtils.nextGaussian();
                    serumLocationsParameter.getParameter(i).setParameterValue(j, r);
                }
            }
        }
    }

    
    
    //load initial
    private void loadInitialLocations(List<String> strainNames, List<String> serumNames) {

		FileReader fileReader;
		try {
			//fileReader = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialCondition/H3N2_mds.virusLocs.log");
			fileReader = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialConditionWithInitialLocationDrift/lastIteration/H3N2_mds.virusLocs.log");
		     /**
		       * Creating a buffered reader to read the file
		       */
		      BufferedReader bReader = new BufferedReader( fileReader);

		      String line;

		      
		      //this routine may give false results if there are extra lines with spaces
		      
		      line = bReader.readLine();
		      System.out.println(line);
		      String namevalue[] = line.split("\t");

		      
		      line = bReader.readLine();
		      System.out.println(line);
		      
		      String datavalue[] = line.split("\t");
		          
		      for (int i = 0; i < virusLocationsParameter.getParameterCount(); i++) {
		    	  
		    	  int index = findStrain( namevalue[i*2+1], strainNames);  //note. namevalue actually has the extra 1 or 2attached to it.. but it doesn't seem to matter
		    //	  System.out.println("name: " + virusLocationsParameter.getParameter(i).getParameterName() + " :" + index);
		    	 // System.out.println(datavalue[i*2+1]);
		    	  virusLocationsParameter.getParameter(index).setParameterValue(0, Double.parseDouble(datavalue[i*2+1]));
		    	  virusLocationsParameter.getParameter(index).setParameterValue(1, Double.parseDouble(datavalue[i*2+2]));
		          //virusLocationsParameter.getParameter(i).setParameterValue(0, 1);
			    	// System.out.print(virusLocationsParameter.getParameter(i).getParameterValue(0) + " ");
			    	// System.out.print(virusLocationsParameter.getParameter(i).getParameterValue(1) + " ");  	  

		      }
		      bReader.close();
		
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}          
 

		FileReader fileReader2;
		try {
			//fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialCondition/H3N2.serumLocs.log");
			fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialConditionWithInitialLocationDrift/lastIteration/H3N2.serumLocs.log");
			
		     /**
		       * Creating a buffered reader to read the file
		       */
		      BufferedReader bReader2 = new BufferedReader( fileReader2);

		      String line;
		      
		      line = bReader2.readLine();
		      System.out.println(line);
		      String namevalue[] = line.split("\t");

		      
		      line = bReader2.readLine();
		      System.out.println(line);
		      
		      String datavalue[] = line.split("\t");
		       //   System.out.println(serumLocationsParameter.getParameterCount());
		      for (int i = 0; i < serumLocationsParameter.getParameterCount(); i++) {
		    	  int index = findStrain( namevalue[i*2+1], serumNames);

		    	 // System.out.println(datavalue[i*2+1]);
		    	  serumLocationsParameter.getParameter(index).setParameterValue(0, Double.parseDouble(datavalue[i*2+1]));
		    	  serumLocationsParameter.getParameter(index).setParameterValue(1, Double.parseDouble(datavalue[i*2+2]));
		          //virusLocationsParameter.getParameter(i).setParameterValue(0, 1);
		   	  
		      }
		      bReader2.close();
		
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}        
      
    	

    }
    

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == virusLocationsParameter) {
            int loc = index / mdsDimension;
            virusLocationChanged[loc] = true;
            if (tipTraitsParameter != null && tipIndices[loc] != -1) {
                Parameter location = virusLocationsParameter.getParameter(loc);
                Parameter tip = tipTraitsParameter.getParameter(tipIndices[loc]);
                int dim = index % mdsDimension;
                tip.setParameterValue(dim, location.getParameterValue(dim));
            }
        } else if (variable == serumLocationsParameter) {
            int loc = index / mdsDimension;
            serumLocationChanged[loc] = true;
        } else if (variable == mdsPrecisionParameter) {
            setLocationChangedFlags(true);
        } else if (variable == locationDriftParameter) {
            setLocationChangedFlags(true);
        } else if (variable == virusDriftParameter) {
                setLocationChangedFlags(true);
        } else if (variable == serumDriftParameter) {
                setLocationChangedFlags(true);
        } else if (variable == serumPotenciesParameter) {
            serumEffectChanged[index] = true;
        } else if (variable == serumBreadthsParameter) {
            serumEffectChanged[index] = true;
        } else if (variable == virusAviditiesParameter) {
            virusEffectChanged[index] = true;
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        System.arraycopy(logLikelihoods, 0, storedLogLikelihoods, 0, logLikelihoods.length);
    }

    @Override
    protected void restoreState() {
        double[] tmp = logLikelihoods;
        logLikelihoods = storedLogLikelihoods;
        storedLogLikelihoods = tmp;

        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
 //uncommenting for testing only

    		
    		//System.exit(0);
      //  if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
     //   }
        
    	//System.out.println("		logL of AGLikelihoodTreeCluster =" + logLikelihood);
    	//System.out.println(" 	p(new)/p(old) = " + Math.exp(logLikelihood - oldLogLikelihood) );
         //oldLogLikelihood = logLikelihood;


// logLikelihood=0;       //for testing purpose only
//System.out.println("logLikelihood of AGLikelihoodCluster= " + logLikelihood);
        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    private double computeLogLikelihood() {
    	    	

        double precision = mdsPrecisionParameter.getParameterValue(0);
        double sd = 1.0 / Math.sqrt(precision);
        
        logLikelihood = 0.0;
        int i = 0;

        for (Measurement measurement : measurements) {

            if (virusLocationChanged[measurement.virus] || serumLocationChanged[measurement.serum] || virusEffectChanged[measurement.virus] || serumEffectChanged[measurement.serum]) {

                double expectation = calculateBaseline(measurement.virus, measurement.serum) - computeDistance(measurement.virus, measurement.serum);
                
                switch (measurement.type) {
                    case INTERVAL: {
                        double minTitre = measurement.log2Titre;
                        double maxTitre = measurement.log2Titre + intervalWidth;
                        logLikelihoods[i] = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
                    } break;
                    case POINT: {
                        logLikelihoods[i] = computeMeasurementLikelihood(measurement.log2Titre, expectation, sd);
                    } break;
                    case THRESHOLD: {
                    	if(measurement.isLowerThreshold){
                    		logLikelihoods[i] = computeMeasurementThresholdLikelihood(measurement.log2Titre, expectation, sd);
                    	}
                    	else{
                    		logLikelihoods[i] = computeMeasurementUpperThresholdLikelihood(measurement.log2Titre, expectation, sd);                  		
                    	}
                    } break;
                    case MISSING:
                        break;
                }
            }
            logLikelihood += logLikelihoods[i];
            i++;
        }
//System.out.println("\nlogLikelihood sum = " + logLikelihood);
        likelihoodKnown = true;

        setLocationChangedFlags(false);
        setSerumEffectChangedFlags(false);
        setVirusEffectChangedFlags(false);

        return logLikelihood;
    }

    private void setLocationChangedFlags(boolean flag) {
        for (int i = 0; i < virusLocationChanged.length; i++) {
            virusLocationChanged[i] = flag;
        }
        for (int i = 0; i < serumLocationChanged.length; i++) {
            serumLocationChanged[i] = flag;
        }
    }

    private void setSerumEffectChangedFlags(boolean flag) {
        for (int i = 0; i < serumEffectChanged.length; i++) {
            serumEffectChanged[i] = flag;
        }
    }

    private void setVirusEffectChangedFlags(boolean flag) {
        for (int i = 0; i < virusEffectChanged.length; i++) {
            virusEffectChanged[i] = flag;
        }
    }

    // offset virus and serum location when computing
    protected double computeDistance(int virus, int serum) {
        Parameter vLoc = virusLocationsParameter.getParameter(virus);
        Parameter sLoc = serumLocationsParameter.getParameter(serum);
        double sum = 0.0;

        // first dimension is shifted
        double vxOffset = 0.0;
        double sxOffset = 0.0;
        if(clusterMeans == true){      	
        	
        	if(virusDriftParameter!= null && virusOffsetsParameter != null && serumOffsetsParameter != null && clusterOffsetsParameter!=null){
//                vxOffset = virusDriftParameter.getParameterValue(0)* clusterOffsetsParameter.getParameterValue(virus);
        		sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
                //vxOffset = locationDriftParameter.getParameterValue(0)*  ;               
           //     System.out.println("clusterOffset =" + clusterOffsetsParameter.getParameterValue(virus));
                 	//System.out.println("offset = " + vxOffset);
                 
        	}
        	
        	//overwrite serum drift
	        if (serumDriftParameter != null && serumOffsetsParameter != null) {
	        //	System.out.println("hihi ya");
	            sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
	        }
	        
        }
        else{
	        if (locationDriftParameter != null && virusOffsetsParameter != null && serumOffsetsParameter != null) {
	            vxOffset = locationDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
                sxOffset = locationDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
	        }
	        if (virusDriftParameter != null && virusOffsetsParameter != null) {
	            vxOffset = virusDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
	        }
	        if (serumDriftParameter != null && serumOffsetsParameter != null) {
	            sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
	        }
        }

        double vxLoc = vLoc.getParameterValue(0) + vxOffset;
        double sxLoc = sLoc.getParameterValue(0) + sxOffset;

       // if(virus ==1){
        //	System.out.println("virus " + virus + " has vxLoc of " + vxLoc + " = " + vLoc.getParameterValue(0) + "+" + vxOffset);
        //}
        
        double difference = vxLoc - sxLoc;
        sum += difference * difference;

        // other dimensions are not
        for (int i = 1; i < mdsDimension; i++) {
            difference = vLoc.getParameterValue(i) - sLoc.getParameterValue(i);
            sum += difference * difference;
        }

        double dist = Math.sqrt(sum);

        if (serumBreadthsParameter != null) {
            double serumBreadth = serumBreadthsParameter.getParameterValue(serum);
            dist /= serumBreadth;
        }
        
        
        
        //if(serum ==0){
        	//System.out.println("The serum location is " + sxLoc +"," + sLoc.getParameterValue(1));
       // }
        
        return(dist);
    }
    
    
    
    
    // offset virus and serum location when computing
    protected double computeDistanceBasedOnArray(int virus, int serum, double[][] virusLocArray) {
       // Parameter vLoc = virusLocationsParameter.getParameter(virus);
        Parameter sLoc = serumLocationsParameter.getParameter(serum);
        double sum = 0.0;

        // first dimension is shifted
        double vxOffset = 0.0;
        double sxOffset = 0.0;
        if(clusterMeans == true){      	
        	
        	
        	if(virusDriftParameter!= null && virusOffsetsParameter != null && serumOffsetsParameter != null && clusterOffsetsParameter!=null){
                vxOffset = virusDriftParameter.getParameterValue(0)* clusterOffsetsParameter.getParameterValue(virus);
        		sxOffset = virusDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
                //vxOffset = locationDriftParameter.getParameterValue(0)*  ;               
           //     System.out.println("clusterOffset =" + clusterOffsetsParameter.getParameterValue(virus));
                 	//System.out.println("offset = " + vxOffset);
                 
        	}
        	
        	//overwrite serum drift
	        if (serumDriftParameter != null && serumOffsetsParameter != null) {
	        //	System.out.println("hihi ya");
	            sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
	        }
	        
        }
        else{
	        if (locationDriftParameter != null && virusOffsetsParameter != null && serumOffsetsParameter != null) {
	            vxOffset = locationDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
                sxOffset = locationDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
	        }
	        if (virusDriftParameter != null && virusOffsetsParameter != null) {
	            vxOffset = virusDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
	        }
	        if (serumDriftParameter != null && serumOffsetsParameter != null) {
	            sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
	        }
        }

        //double vxLoc = vLoc.getParameterValue(0) + vxOffset;
        double vxLoc = virusLocArray[virus][0] + vxOffset;
        double sxLoc = sLoc.getParameterValue(0) + sxOffset;

       // if(virus ==1){
        //	System.out.println("virus " + virus + " has vxLoc of " + vxLoc + " = " + vLoc.getParameterValue(0) + "+" + vxOffset);
        //}
        
        double difference = vxLoc - sxLoc;
        sum += difference * difference;

        // other dimensions are not
        for (int i = 1; i < mdsDimension; i++) {
           // difference = vLoc.getParameterValue(i) - sLoc.getParameterValue(i);
        	 difference = virusLocArray[virus][i] - sLoc.getParameterValue(i);
            sum += difference * difference;
        }

        double dist = Math.sqrt(sum);

        if (serumBreadthsParameter != null) {
            double serumBreadth = serumBreadthsParameter.getParameterValue(serum);
            dist /= serumBreadth;
        }

        return dist;
    }
    
    


    // Calculates the expected log2 titre when mapDistance = 0
    private double calculateBaseline(int virus, int serum) {
        double baseline = serumPotenciesParameter.getParameterValue(serum);
        if (virusAviditiesParameter != null) {
            baseline += virusAviditiesParameter.getParameterValue(virus);
        }
        return baseline;
    }

    private static double computeMeasurementLikelihood(double titre, double expectation, double sd) {

        double lnL = NormalDistribution.logPdf(titre, expectation, sd);

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite point measurement");
        }
        return lnL;
    }

    private static double computeMeasurementThresholdLikelihood(double titre, double expectation, double sd) {

        // real titre is somewhere between -infinity and measured 'titre'
        // want the lower tail of the normal CDF

        double lnL = NormalDistribution.cdf(titre, expectation, sd, true);          // returns logged CDF

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite threshold measurement");
        }
        return lnL;
    }


    private static double computeMeasurementUpperThresholdLikelihood(double titre, double expectation, double sd) {

        // real titre is somewhere between -infinity and measured 'titre'
        // want the lower tail of the normal CDF
    	double L = NormalDistribution.cdf(titre, expectation, sd, false);          // returns  CDF
    	double lnL = Math.log(1-L);  //get the upper tail probability, then log it

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite threshold measurement");
        }
        return lnL;
    }    
    
    
    private static double computeMeasurementIntervalLikelihood(double minTitre, double maxTitre, double expectation, double sd) {

        // real titre is somewhere between measured minTitre and maxTitre

        double cdf1 = NormalDistribution.cdf(maxTitre, expectation, sd, true);     // returns logged CDF
        double cdf2 = NormalDistribution.cdf(minTitre, expectation, sd, true);     // returns logged CDF
        double lnL = LogTricks.logDiff(cdf1, cdf2);

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            // this occurs when the interval is in the far tail of the distribution, cdf1 == cdf2
            // instead return logPDF of the point
            lnL = NormalDistribution.logPdf(minTitre, expectation, sd);
            if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
                throw new RuntimeException("infinite interval measurement");
            }
        }
        return lnL;
    }

    public void makeDirty() {
        likelihoodKnown = false;
        setLocationChangedFlags(true);
    }

    private class Measurement {
        private Measurement(final int virus, final int serum, final double virusDate, final double serumDate, final MeasurementType type, final double titre, final boolean isLowerThreshold) {
            this.virus = virus;
            this.serum = serum;
            this.virusDate = virusDate;
            this.serumDate = serumDate;
            this.type = type;
            this.titre = titre;
            this.log2Titre = Math.log(titre) / Math.log(2);
            this.isLowerThreshold = isLowerThreshold;
        }

        final int virus;
        final int serum;
        final double virusDate;
        final double serumDate;
        final MeasurementType type;
        final double titre;
        final double log2Titre;
        final boolean isLowerThreshold;

    };

    private final List<Measurement> measurements = new ArrayList<Measurement>();
    private final List<String> virusNames = new ArrayList<String>();
    private final List<String> serumNames = new ArrayList<String>();
    private final List<Double> virusDates = new ArrayList<Double>();
    private final List<Double> serumDates = new ArrayList<Double>();

    private final int mdsDimension;
    private final double intervalWidth;
    private final Parameter mdsPrecisionParameter;
    private final Parameter locationDriftParameter;
    private final Parameter virusDriftParameter;
    private final Parameter serumDriftParameter;

    private final MatrixParameter virusLocationsParameter;
    private final MatrixParameter serumLocationsParameter;

    private final Parameter virusOffsetsParameter;
    private final Parameter serumOffsetsParameter;

    private final CompoundParameter tipTraitsParameter;
    private int[] tipIndices;

    private final Parameter virusAviditiesParameter;
    private final Parameter serumPotenciesParameter;
    private final Parameter serumBreadthsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    private final boolean[] virusLocationChanged;
    private final boolean[] serumLocationChanged;
    private final boolean[] serumEffectChanged;
    private final boolean[] virusEffectChanged;
    private double[] logLikelihoods;
    private double[] storedLogLikelihoods;
    
    private boolean clusterMeans = false;
    private Parameter clusterOffsetsParameter;

// **************************************************************
// XMLObjectParser
// **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String FILE_NAME = "fileName";
        public final static String TIP_TRAIT = "tipTrait";
        public final static String VIRUS_LOCATIONS = "virusLocations";
        public final static String SERUM_LOCATIONS = "serumLocations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MERGE_SERUM_ISOLATES = "mergeSerumIsolates";
        public static final String DRIFT_INITIAL_LOCATIONS = "driftInitialLocations";
        public static final String INTERVAL_WIDTH = "intervalWidth";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String LOCATION_DRIFT = "locationDrift";
        public static final String VIRUS_DRIFT = "virusDrift";
        public static final String SERUM_DRIFT = "serumDrift";
        public static final String VIRUS_AVIDITIES = "virusAvidities";
        public static final String SERUM_POTENCIES = "serumPotencies";
        public static final String SERUM_BREADTHS = "serumBreadths";
        public final static String VIRUS_OFFSETS = "virusOffsets";
        public final static String SERUM_OFFSETS = "serumOffsets";
        public final static String CLUSTER_MEANS = "clusterMeans";
        public final static String CLUSTER_OFFSETS = "clusterOffsetsParameter";

        public String getParserName() {
            return AG_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<String[]> assayTable;
            try {
                assayTable = DataTable.Text.parse(new FileReader(fileName), true, false);
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }
            System.out.println("Loaded HI table file: " + fileName);

            boolean mergeSerumIsolates = xo.getAttribute(MERGE_SERUM_ISOLATES, false);
            
            boolean cluster_means = xo.getAttribute(CLUSTER_MEANS, false);

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);
            double intervalWidth = 0.0;
            if (xo.hasAttribute(INTERVAL_WIDTH)) {
                intervalWidth = xo.getDoubleAttribute(INTERVAL_WIDTH);
            }

            double driftInitialLocations = 0.0;
            if (xo.hasAttribute(DRIFT_INITIAL_LOCATIONS)) {
                driftInitialLocations = xo.getDoubleAttribute(DRIFT_INITIAL_LOCATIONS);
            }

            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter virusLocationsParameter = null;
            if (xo.hasChildNamed(VIRUS_LOCATIONS)) {
                virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
            }

            MatrixParameter serumLocationsParameter = null;
            if (xo.hasChildNamed(SERUM_LOCATIONS)) {
                serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);
            }

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter locationDrift = null;
            if (xo.hasChildNamed(LOCATION_DRIFT)) {
                locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
            }

            Parameter virusDrift = null;
            if (xo.hasChildNamed(VIRUS_DRIFT)) {
            	virusDrift = (Parameter) xo.getElementFirstChild(VIRUS_DRIFT);
            }

            Parameter serumDrift = null;
            if (xo.hasChildNamed(SERUM_DRIFT)) {
            	serumDrift = (Parameter) xo.getElementFirstChild(SERUM_DRIFT);
            }

            Parameter virusOffsetsParameter = null;
            if (xo.hasChildNamed(VIRUS_OFFSETS)) {
                virusOffsetsParameter = (Parameter) xo.getElementFirstChild(VIRUS_OFFSETS);
            }

            Parameter serumOffsetsParameter = null;
            if (xo.hasChildNamed(SERUM_OFFSETS)) {
                serumOffsetsParameter = (Parameter) xo.getElementFirstChild(SERUM_OFFSETS);
            }

            Parameter serumPotenciesParameter = null;
            if (xo.hasChildNamed(SERUM_POTENCIES)) {
                serumPotenciesParameter = (Parameter) xo.getElementFirstChild(SERUM_POTENCIES);
            }

            Parameter serumBreadthsParameter = null;
            if (xo.hasChildNamed(SERUM_BREADTHS)) {
                serumBreadthsParameter = (Parameter) xo.getElementFirstChild(SERUM_BREADTHS);
            }

            Parameter virusAviditiesParameter = null;
            if (xo.hasChildNamed(VIRUS_AVIDITIES)) {
                virusAviditiesParameter = (Parameter) xo.getElementFirstChild(VIRUS_AVIDITIES);
            }
            
            Parameter clusterOffsetsParameter = null;
            if (xo.hasChildNamed(CLUSTER_OFFSETS)) {
            	clusterOffsetsParameter = (Parameter) xo.getElementFirstChild(CLUSTER_OFFSETS);
            }


            AGLikelihoodTreeCluster AGL = new AGLikelihoodTreeCluster(
                    mdsDimension,
                    mdsPrecision,
                    locationDrift,
                    virusDrift,
                    serumDrift,
                    virusLocationsParameter,
                    serumLocationsParameter,
                    tipTraitParameter,
                    virusOffsetsParameter,
                    serumOffsetsParameter,
                    serumPotenciesParameter,
                    serumBreadthsParameter,
                    virusAviditiesParameter,
                    assayTable,
                    mergeSerumIsolates,
                    intervalWidth,
                    driftInitialLocations, 
                    cluster_means, 
                    clusterOffsetsParameter);
                        

            Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Utils.getCitationString(AGL));

            return AGL;
        }

//************************************************************************
// AbstractXMLObjectParser implementation
//************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of immunological assay data such as Hemagglutinin inhibition (HI) given vectors of coordinates" +
                    "for viruses and sera/antisera in some multidimensional 'antigenic' space.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                AttributeRule.newBooleanRule(MERGE_SERUM_ISOLATES, true, "Should multiple serum isolates from the same strain have their locations merged (defaults to false)"),
                AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
                AttributeRule.newDoubleRule(DRIFT_INITIAL_LOCATIONS, true, "The degree to drift initial virus and serum locations, defaults to 0.0"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "Optional parameter of tip locations from the tree", true),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "Parameter of locations of all virus"),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "Parameter of locations of all sera"),
                new ElementRule(VIRUS_OFFSETS, Parameter.class, "Optional parameter for virus dates to be stored", true),
                new ElementRule(SERUM_OFFSETS, Parameter.class, "Optional parameter for serum dates to be stored", true),
                new ElementRule(SERUM_POTENCIES, Parameter.class, "Optional parameter for serum potencies", true),
                new ElementRule(SERUM_BREADTHS, Parameter.class, "Optional parameter for serum breadths", true),
                new ElementRule(VIRUS_AVIDITIES, Parameter.class, "Optional parameter for virus avidities", true),
                new ElementRule(MDS_PRECISION, Parameter.class, "Parameter for precision of MDS embedding"),
                new ElementRule(LOCATION_DRIFT, Parameter.class, "Optional parameter for drifting locations with time", true),
                new ElementRule(VIRUS_DRIFT, Parameter.class, "Optional parameter for drifting only virus locations, overrides locationDrift", true),
                new ElementRule(SERUM_DRIFT, Parameter.class, "Optional parameter for drifting only serum locations, overrides locationDrift", true),
                AttributeRule.newBooleanRule(CLUSTER_MEANS, true, "Should we use cluster means to control the virus locations"),
               new ElementRule(CLUSTER_OFFSETS, Parameter.class, "Parameter of cluster offsets of all virus"),                
        };

        public Class getReturnType() {
            return AGLikelihoodTreeCluster.class;
        }
    };

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian Antigenic Cartography framework";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{
                                new Author("C", "Cheung"),
                                new Author("A", "Rambaut"),
                                new Author("P", "Lemey"),
                                new Author("MA", "Suchard"),
                                new Author("T", "Bedford")
                        },
                        Citation.Status.IN_PREPARATION
                ),
                CommonCitations.BEDFORD_2015_INTEGRATING);
    }
    
    public double getLogLikelihoodBasedOnPrecompute(int[] clusterLabel, int numClusters, int[] oldObservationCluster, double[] oldContribution, int[] newObservationCluster, double[] newContribution) {
	   	
        double precision = mdsPrecisionParameter.getParameterValue(0);
        double sd = 1.0 / Math.sqrt(precision);
    			
	     for (int i=0; i < measurements.size(); i++) {
	    	 newObservationCluster[i] = clusterLabel[measurements.get(i).virus];
	     }
	     
		double logL=0;
		for(int i=0; i< measurements.size(); i++ ){
			if(oldObservationCluster != null && newObservationCluster[i] == oldObservationCluster[i]){
				
				double newC = computeContribution(measurements.get(i), precision, sd);
				if(newC != oldContribution[i]){
					System.out.println("newObservationCluster[i]=" + newObservationCluster[i] + " and old =" + oldObservationCluster[i]);
					System.out.println("old contribution of i=" + i + " is " + oldContribution[i] + " but new is " + newC);
					
					System.out.println("They should be the same. Why are they different?");
					System.exit(0);
				}
				
			//	System.out.println("run!");
				newContribution[i] = oldContribution[i];				
			}
			else{
				//recompute.
				newContribution[i] = computeContribution(measurements.get(i), precision, sd);
			}
			logL += newContribution[i];
			
		}
		
		
			
		return(logL);
	}
    
    
    
    public double computeContribution(Measurement measurement, double precision, double sd){
    		double curLogL=0;
        

            double expectation = calculateBaseline(measurement.virus, measurement.serum) - computeDistance(measurement.virus, measurement.serum);
            switch (measurement.type) {
                case INTERVAL: {

                    double minTitre = measurement.log2Titre;
                    double maxTitre = measurement.log2Titre + intervalWidth;
                    curLogL = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
                } break;
                case POINT: {

                	curLogL = computeMeasurementLikelihood(measurement.log2Titre, expectation, sd);
                } break;
                case THRESHOLD: {

                	if(measurement.isLowerThreshold){
                		curLogL = computeMeasurementThresholdLikelihood(measurement.log2Titre, expectation, sd);
                	}
                	else{
                		curLogL = computeMeasurementUpperThresholdLikelihood(measurement.log2Titre, expectation, sd);                  		
                	}
                } break;
                case MISSING:
                    break;
            }
        

    return curLogL;
    }
    
    
    
    
    
    
    public double getClusterLogLikelihoodUpdate(int[] clusterLabel, int numClusters, int[] oldObservationCluster, double[] oldContribution, int[] newObservationCluster, double[] newContribution, double[] oldClusterSum, double[] newClusterSum) {

        double precision = mdsPrecisionParameter.getParameterValue(0);
        double sd = 1.0 / Math.sqrt(precision);
    	
       	
		List<Measurement>[] partition = new List[clusterLabel.length]; // the size may be an overkill, but it's ok
		for(int i=0; i < clusterLabel.length; i++){
			partition[i] = new LinkedList();
		}        
    
		int []needUpdateCluster = new int[numClusters];
		
		if(oldObservationCluster == null){
			for (int i=0; i < measurements.size(); i++) {
		    	Measurement m = measurements.get(i);
		    	newObservationCluster[i] = clusterLabel[m.virus];
		    	partition[clusterLabel[m.virus]].add(m);  // this is correct
		     }
			for(int i=0; i<numClusters; i++){
				needUpdateCluster[i] = 1;
			}
		}
		else{
	     for (int i=0; i < measurements.size(); i++) {
	    	 Measurement m = measurements.get(i);
	    	 newObservationCluster[i] = clusterLabel[m.virus];
	    	 partition[clusterLabel[m.virus]].add(m);  // this is correct
	    	 
	    	 if(newObservationCluster[i] != oldObservationCluster[i]){
	    		 needUpdateCluster[newObservationCluster[i]] = 1;
	    		 needUpdateCluster[oldObservationCluster[i]] = 1;
	    	 }
	     }
		}
		
		
		//for(int i=0; i < needUpdateCluster.length; i++){
		//	System.out.println("update cluster " +i + "="+ needUpdateCluster[i] );
		//}
	     //System.out.println("==================");
		
		//Method 1
		/*
		double logL=0;
		for(int i=0; i< measurements.size(); i++ ){
			if(oldContribution != null && needUpdateCluster[newObservationCluster[i]]==0){
				//System.out.println("don't need to update cluster");
				newContribution[i] = oldContribution[i];				
			}
			else{
				//recompute.
				newContribution[i] = computeContribution(measurements.get(i), precision, sd);
			}
			logL += newContribution[i];
			
		}
		*/
		
		
		//Method 2
		
		double logL=0;
   		for(int i=0; i < numClusters; i++ ){			
			if(partition[i].size() >0){
				if(needUpdateCluster[i] == 0 && oldClusterSum != null){
					//update the cluster i's contribution
					newClusterSum[i] += oldClusterSum[i];
				}
				else{
					newClusterSum[i] += computeLikelihoodBasedOnClusters(partition[i]);
				}
				logL += newClusterSum[i];
			}
			else{
				newClusterSum[i] = 0;
			}
    		
       	}//for	
		
		
		return(logL);
	}
        
        
    
    
    
    public double getClusterLogLikelihood(int[] clusterLabel, int numClusters) {
	
	//public double getClusterLogLikelihood(double[][] vLoc, double[] mu0_offset,
		//	int[] clusterLabel, int numClusters, int[] needUpdateCluster) {
		
		//double logL = computeLogLikelihood();
    	
		List<Measurement>[] partition = new List[clusterLabel.length]; // the size may be an overkill, but it's ok
		for(int i=0; i < clusterLabel.length; i++){
			partition[i] = new LinkedList();
		}
		//partition measurements into clusters using viruses
	//	int j=0;
	//	for(Measurement measurement: measurements){
		//	if(j < 5000){
			//	partition[0].add(measurement);
			//}
			//else{
			//	partition[1].add(measurement);
		//	}
		//	j++;
		//}

		
		
	     for (Measurement measurement : measurements) {
	    	 //System.out.println("partition #: " + clusterLabel[measurement.virus]);
	    	 partition[clusterLabel[measurement.virus]].add(measurement);  // this is correct
	     }
	     
		
		//System.exit(0);
		
		double logL=0;
		// TODO Auto-generated method stub
		for(int i=0; i < numClusters; i++ ){
		//for(int i=0; i < clusterLabel.length; i++ ){
			if(partition[i].size() >0){
				//System.out.println("partition size is = " + partition[i].size());
				logL += computeLikelihoodBasedOnClusters(partition[i]);
			}
			//if(needUpdateCluster[i] == 1){
				//update the cluster i's contribution
				//logClusterL[i] = 1;
			//}
			//calculate logClusterL[i]
			
			
			
			//logL += logClusterL[i]; 
		}
		
		
		
		return(logL);
	}
    
    

    
    

	private double computeLikelihoodBasedOnClusters(
			List<Measurement> linkedList) {
  	

        double precision = mdsPrecisionParameter.getParameterValue(0);
        double sd = 1.0 / Math.sqrt(precision);

        double logL = 0.0;
        int i = 0;

        for (Measurement measurement : linkedList) {
        		double curLogL=0;
            

                double expectation = calculateBaseline(measurement.virus, measurement.serum) - computeDistance(measurement.virus, measurement.serum);
                switch (measurement.type) {
                    case INTERVAL: {

                        double minTitre = measurement.log2Titre;
                        double maxTitre = measurement.log2Titre + intervalWidth;
                        curLogL = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
                    } break;
                    case POINT: {

                    	curLogL = computeMeasurementLikelihood(measurement.log2Titre, expectation, sd);
                    } break;
                    case THRESHOLD: {

                    	if(measurement.isLowerThreshold){
                    		curLogL = computeMeasurementThresholdLikelihood(measurement.log2Titre, expectation, sd);
                    	}
                    	else{
                    		curLogL = computeMeasurementUpperThresholdLikelihood(measurement.log2Titre, expectation, sd);                  		
                    	}
                    } break;
                    case MISSING:
                        break;
                }
            
            logL += curLogL;
            i++;
               // System.out.println("curLogL = " + curLogL);
        }

//        likelihoodKnown = true;

//        setLocationChangedFlags(false);
 //       setSerumEffectChangedFlags(false);
  //      setVirusEffectChangedFlags(false);
//System.out.println("logL = " + logL);
        return logL;
    }

	public int getNumObservations() {
		return measurements.size();
	}



	public double getNumSera() {
		return serumLocationsParameter.getParameterCount();
	}



	public MatrixParameter getSerumLocationsParameter() {
		return serumLocationsParameter;
	}




	public void printVirusLocations() {
		
		for(int i=0; i < virusLocationsParameter.getParameterCount(); i++){
		 Parameter v = virusLocationsParameter.getParameter(i);
		 System.out.println(v.getId() + " " + v.getParameterValue(0) + "," + v.getParameterValue(1));
		}
		
	}
	
	

}
