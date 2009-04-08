package dr.evomodel.coalescent;

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

public abstract class AbstractFixedGridLogger extends MCLogger{

	public static final String GRID_STOP_TIME = "gridStopTime";
	public static final String NUMBER_OF_INTERVALS = "numberOfIntervals";
	
	private double[] gridHeights;
	private int intervalNumber;
	
	public AbstractFixedGridLogger(LogFormatter formatter, int logEvery,
			double gridHeight, int intervalNumber) {
		super(formatter, logEvery, false);
		
		this.intervalNumber = intervalNumber;
		
		gridHeights = new double[intervalNumber + 1];
		for(int i = 1; i < gridHeights.length; i++){
			gridHeights[i] = gridHeights[i - 1] + gridHeight/(double)intervalNumber;
		}
	}
	
	public int getIntervalNumber(){
		return intervalNumber;
	}
	
	public double getGridStopTime(){
		return gridHeights[gridHeights.length - 1];
	}
	
	private String getProjectedPopSize(double[] bslPopSize,double[] groupHeights, 
			   double intervalStart,double intervalEnd){
		int startGroup = 0;
		int endGroup = 0;

		while(startGroup < groupHeights.length && intervalStart >= groupHeights[startGroup]){
			startGroup++;
		}
		while(endGroup < groupHeights.length && intervalEnd > groupHeights[endGroup]){
			endGroup++;
		}
			
		
		if(startGroup >= groupHeights.length){
			return "" + bslPopSize[bslPopSize.length - 1];
		}

		if(startGroup == endGroup){
			return "" + bslPopSize[startGroup];
		}

		double value = 0;

		if(endGroup == groupHeights.length){

			value += (groupHeights[startGroup] - intervalStart)*bslPopSize[startGroup];

			for(int i = startGroup + 1; i < groupHeights.length; i++){
				value += bslPopSize[i]*(groupHeights[i] - groupHeights[i-1]);
			}

			return "" + value/(double)(groupHeights[groupHeights.length - 1] - intervalStart);
		}

		value += bslPopSize[startGroup]*(groupHeights[startGroup] - intervalStart);
		value += bslPopSize[endGroup]*(intervalEnd - groupHeights[endGroup - 1]);

		for(int i = startGroup + 1; i < endGroup; i++){
			value += bslPopSize[i]*(groupHeights[i] - groupHeights[i-1]);
		}

		return "" + value/(double)(intervalEnd - intervalStart);

	}
	
	private void getStringGrid(String[] values){
		
		double[] bslPopSize = getPopSizes();
		double[] groupHeight = getCoalescentHeights();
				
		for(int i = 0; i < intervalNumber; i++){
			values[i + 1] = getProjectedPopSize(bslPopSize,groupHeight,gridHeights[i],gridHeights[i+1]);
		}
	}
	
	public void log(int state){
		if (logEvery <= 0 || ((state % logEvery) == 0)) {
			String[] values = new String[intervalNumber + 2];
		
			values[0] = Integer.toString(state);
			
			getStringGrid(values);
			
			values[values.length - 1] = "" + getAdditionalDensity();
						
			logValues(values);
		}
	}
	
	abstract double[] getPopSizes();
	abstract double[] getCoalescentHeights();
	abstract double   getAdditionalDensity();
}
