package dr.app.gui.chart;

/**
 * make x and y together as the key of map
 *
 * @author Walter Xie
 */
public class XY implements Comparable<XY> {
    final double x, y;

    XY(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public boolean equals(Object o){
        if(this==o)return true;
        if(o==null)return false;
        if(o instanceof XY){
            XY xy = (XY)o;
            return this.x==xy.x && this.y==xy.y;
        }
        return false;
    }

    public int hashCode(){
        return (int) (x * 31 + y);
    }

    //http://stackoverflow.com/questions/9307751/override-compareto-and-sort-using-two-strings
    @Override
    public int compareTo(final XY o) {
        int cmp = Double.compare(this.x, o.x);
        if (cmp == 0) cmp = Double.compare(this.y, o.y);
        return cmp;
    }

//    public String toString() {
//        return Double.toString(x) + "|" + Double.toString(y);
//    }
}
