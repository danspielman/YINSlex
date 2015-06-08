package lexvolts;

import java.util.Random;

/**
 * Created by sushant on 1/18/15.
 */
public class EdgeGradientHelper {
    public static double tol = 1e-15; // Tolerance for comparisons

    int maxSize;

    double[] leftV;
    double[] leftD;
    double[] rightV;
    double[] rightD;
    int leftCount;
    int rightCount;

    Random random;

    public EdgeGradientHelper(int maxSize, long seed) {
        this.maxSize = maxSize;
        leftV = new double[maxSize];
        leftD = new double[maxSize];
        rightV = new double[maxSize];
        rightD = new double[maxSize];
        leftCount = 0;
        rightCount = 0;
        random = new Random(seed);
    }

    // Add a terminal on the left with voltage volt, at distance dis
    public void addLeftTerm(double volt, double dis) {
        assert (leftCount < maxSize);
        leftV[leftCount] = volt;
        leftD[leftCount] = dis;
        leftCount++;
    }

    // Add a terminal on the right with voltage volt, at distance dis
    public void addRightTerm(double volt, double dis) {
        assert (rightCount < maxSize);
        rightV[rightCount] = volt;
        rightD[rightCount] = dis;
        rightCount++;
    }

    // Resets the arrays
    // NOTE: MUST be called whenever we start working with a new edge
    private void reset() {
        leftCount = 0;
        rightCount = 0;
    }

    private double leftTermGrad(int term, double edgeLength) {
        assert (term < leftCount);
        double grad = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < rightCount; i++) {
            double testGrad = (rightV[i] - leftV[term]) / (rightD[i] + edgeLength + leftD[term]);
            if (testGrad < 0)
                testGrad = -testGrad;
            if (grad < testGrad)
                grad = testGrad;
        }
        return grad;
    }

//    private boolean robustGreaterEq(double d1, double d2){
//        double abs = Math.abs(d1) + Math.abs(d2);
////        double scale = Math.pow(2, StrictMath.getExponent(abs));
//        double scale = Double.longBitsToDouble(Double.doubleToLongBits(abs) & 0x7ff0000000000000L);
//        // This line just computes the exponent of the double corresponding to abs
////        System.out.println("scale : " + scale + " abs : " + abs);
//        if(d1 > (d2 + scale*tol))
//            return true;
//        else
//            return false;
//    }

    // Takes the length of the edge, and returns the gradient of the edge
    // Implements a randomized algorithm
    public double computeRandGrad(double length) {
        if ((leftCount == 0) || (rightCount == 0))
            throw new Error("Did you call computeGrad twice? You must re-populate the terminals after computing terminals once.");

        // Compute min voltage
        double minV = Double.POSITIVE_INFINITY;
        for (int i = 0; i < leftCount; i++) {
            if (minV > leftV[i])
                minV = leftV[i];
        }
        for (int i = 0; i < rightCount; i++) {
            if (minV > rightV[i])
                minV = rightV[i];
        }

        // Subtract min voltage
        for (int i = 0; i < leftCount; i++)
            leftV[i] -= minV;
        for (int i = 0; i < rightCount; i++)
            rightV[i] -= minV;

        int randLeftTerm = random.nextInt(leftCount);
        double leftRandGrad = leftTermGrad(randLeftTerm, length);
        double testGrad = (1 + tol) * leftRandGrad;

        boolean[] leftTermSteeperPath = new boolean[leftCount];
        for (int i = 0; i < leftCount; i++)
            leftTermSteeperPath[i] = false;
//

        double vLow2 = Double.POSITIVE_INFINITY;
        for (int i = 0; i < rightCount; i++) {
            double testVLow2 = rightV[i] + testGrad * rightD[i];
            if (vLow2 > testVLow2)
                vLow2 = testVLow2;
        }
        double vLow1 = vLow2 + testGrad * length;

        for (int i = 0; i < leftCount; i++) {
            if (leftV[i] > vLow1 + testGrad * leftD[i])
//            if(robustGreaterEq(leftV[i],vLow1 + testGrad * leftD[i]))
                leftTermSteeperPath[i] = true;
        }

        double vHigh2 = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < rightCount; i++) {
            double testVHigh2 = rightV[i] - testGrad * rightD[i];
            if (vHigh2 < testVHigh2)
                vHigh2 = testVHigh2;
        }
        double vHigh1 = vHigh2 - testGrad * length;

        for (int i = 0; i < leftCount; i++) {
            if (leftV[i] < vHigh1 - testGrad * leftD[i])
//            if(robustGreaterEq(vHigh1, leftV[i] + testGrad*leftD[i]))
                leftTermSteeperPath[i] = true;
        }

        // Test if each left terminal is STRICTLY ABOVE the gradient
        // if yes, put it at the top of the list
        // increment newLeftCount

        int newLeftCount = 0;
        for (int i = 0; i < leftCount; i++) {
            if (leftTermSteeperPath[i]) {
                leftV[newLeftCount] = leftV[i];
                leftD[newLeftCount] = leftD[i];
                newLeftCount++;
            }
        }
        if (newLeftCount == 0) {
            reset();
            return leftRandGrad;
        }

        leftCount = newLeftCount;

        // Swap left and right
        int maxSwap = leftCount;
        if (maxSwap < rightCount)
            maxSwap = rightCount;

        double swapV, swapD;
        for (int i = 0; i < maxSwap; i++) {
            swapD = leftD[i];
            leftD[i] = rightD[i];
            rightD[i] = swapD;

            swapV = leftV[i];
            leftV[i] = rightV[i];
            rightV[i] = swapV;
        }

        int swapCount = leftCount;
        leftCount = rightCount;
        rightCount = swapCount;

        return computeRandGrad(length);
    }

}
