package dr.inference.model;

/**
 * a class that stores the sum of the diagonal elements and all elements of a matrix
 */
class TreeVarianceSums {

    private double diagonalSum;
    private double totalSum;

    protected TreeVarianceSums(double diagonalSum, double totalSum) {

        this.diagonalSum = diagonalSum;
        this.totalSum = totalSum;
    }

    public double getDiagonalSum() {
        return diagonalSum;
    }

    public double getTotalSum() {
        return totalSum;
    }

    public void setDiagonalSum(double diagonalSum) {
        this.diagonalSum = diagonalSum;
    }

    public void setTotalSum(double totalSum) {
        this.totalSum = totalSum;
    }
}
