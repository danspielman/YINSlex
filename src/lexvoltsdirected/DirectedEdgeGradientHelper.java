package lexvoltsdirected;

import java.util.Random;

/**
 * Created by sushant on 1/31/15.
 */
public class DirectedEdgeGradientHelper {
    public static double tol = 1e-15; // Tolerance for comparisons

    int maxSize;

    // These four arrays store information about the end terminals for paths going through this vertex
    double[] srcV; // Voltages of src of such paths
    double[] srcD; // Distances to those terminals
    double[] destV; // Voltages of dests for such paths
    double[] destD; // Distances to those terminals
    int srcCount;
    int destCount;

    Random random;

    public DirectedEdgeGradientHelper(int maxSize, long seed) {
        this.maxSize = maxSize;
        srcV = new double[maxSize];
        srcD = new double[maxSize];
        destV = new double[maxSize];
        destD = new double[maxSize];
        srcCount = 0;
        destCount = 0;
        random = new Random(seed);
    }

    // Add a terminal on the left with voltage volt, at distance dis
    public void addSrc(double volt, double dis) {
        assert (srcCount < maxSize);
        srcV[srcCount] = volt;
        srcD[srcCount] = dis;
        srcCount++;
    }

    // Add a terminal on the right with voltage volt, at distance dis
    public void addDest(double volt, double dis) {
        assert (destCount < maxSize);
        destV[destCount] = volt;
        destD[destCount] = dis;
        destCount++;
    }

    // Resets the arrays
    // NOTE: MUST be called whenever we start working with a new edge
    private void reset() {
        srcCount = 0;
        destCount = 0;
    }

    private double srcGrad(int src, double edgeLength) {
        assert (src < srcCount);
        double grad = 0;
        for (int i = 0; i < destCount; i++) {
            double testGrad = (srcV[src] - destV[i]) / (srcD[src] + edgeLength + destD[i]);
            if (testGrad < 0)
                testGrad = 0;
            if (grad < testGrad)
                grad = testGrad;
        }
        return grad;
    }

    private double destGrad(int dest, double edgeLength) {
        assert (dest < destCount);
        double grad = 0;
        for (int i = 0; i < srcCount; i++) {
            double testGrad = (srcV[i] - destV[dest]) / (srcD[i] + edgeLength + destD[dest]);
            if (testGrad < 0)
                testGrad = 0;
            if (grad < testGrad)
                grad = testGrad;
        }
        return grad;
    }

    // Takes the length of the edge, and returns the gradient of the edge
    // Implements a randomized algorithm
    public double computeRandGrad(double length) {
        if ((srcCount == 0) || (destCount == 0)) {
            reset();
            return 0;
        }

        // Compute min voltage
        double minV = Double.POSITIVE_INFINITY;
        for (int i = 0; i < srcCount; i++) {
            if (minV > srcV[i])
                minV = srcV[i];
        }
        for (int i = 0; i < destCount; i++) {
            if (minV > destV[i])
                minV = destV[i];
        }

        // Subtract min voltage
        for (int i = 0; i < srcCount; i++)
            srcV[i] -= minV;
        for (int i = 0; i < destCount; i++)
            destV[i] -= minV;

        int randSrc = random.nextInt(srcCount);
        int randDest = random.nextInt(destCount);
        double srcRandGrad = srcGrad(randSrc, length);
        double destRandGrad = destGrad(randDest, length);
        double randGrad = srcRandGrad;
        if (randGrad < destRandGrad)
            randGrad = destRandGrad;
        double testGrad = (1 + tol) * randGrad;

        double vLow2 = Double.POSITIVE_INFINITY;
        for (int i = 0; i < destCount; i++) {
            double testVLow2 = destV[i] + testGrad * destD[i];
            if (vLow2 > testVLow2)
                vLow2 = testVLow2;
        }
        double vLow1 = vLow2 + testGrad * length;

        int newSrcCount = 0;
        for (int i = 0; i < srcCount; i++) {
            if (srcV[i] > vLow1 + testGrad * srcD[i] + tol) {
                srcV[newSrcCount] = srcV[i];
                srcD[newSrcCount] = srcD[i];
            }
            newSrcCount++;
        }
        if (newSrcCount == 0) {
            reset();
            return randGrad;
        }

        srcCount = newSrcCount;

        double vHigh1 = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < srcCount; i++) {
            double testVHigh1 = srcV[i] - testGrad * srcD[i];
            if (vHigh1 < testVHigh1)
                vHigh1 = testVHigh1;
        }
        double vHigh2 = vHigh1 - testGrad * length;

        int newDestCount = 0;
        for (int i = 0; i < destCount; i++) {
            if (destV[i] < vHigh2 - testGrad * destD[i] - tol) {
                destV[newDestCount] = destV[i];
                destD[newDestCount] = destD[i];
                newDestCount++;
            }
        }
        if (newDestCount == 0) {
            reset();
            return randGrad;
        }
        destCount = newDestCount;

        return computeRandGrad(length);
    }

}
