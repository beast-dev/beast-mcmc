package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import java.util.Comparator;

public class OrderDouble implements Comparator<OrderDouble>, Comparable<OrderDouble>{
	   private Integer index;
	   private double value;
	   private double value2;
	   
	   OrderDouble(){
	   }

	   OrderDouble(Integer i, double v){
	      index = i;
	      value = v;
	   }
	   
	   OrderDouble(Integer i, double v, double a){
		      index = i;
		      value = v;
		      value2 = a;
		   }

	   public Integer getIndex(){
	      return index;
	   }

	   public double getValue(){
	      return value;
	   }
	   public double getValue2(){
		   return value2;
	   }
	   // Overriding the compareTo method
	   public int compareTo(OrderDouble d){
	      return (this.index).compareTo(d.index);
	   }

	   // Overriding the compare method to sort the value
	   //Returns a negative integer, zero, or a positive integer as the first argument is less than, equal to, or greater than the second.
	   public int compare(OrderDouble d, OrderDouble d1){
		  if(d.value - d1.value > 0){
	      return 1;
		  }
		  else if(d.value - d1.value < 0){
			  return -1;
		  }
		  else{
			  return 0;
		  }
	   }
	}
