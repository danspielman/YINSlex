package lexvoltsdirected;

import java.util.Random;

/**
 * Created by sushant on 1/31/15.
 */
public class RandLexVoltDirected {

    //    public static double eps = 1e-6; // 1+eps approximation is okay
    public static double testEps = 1e-10; // We test that the final voltage satisfies |voltage-testvoltage| <= testEps * scale, where
    // testvoltage is the voltage predicted by neighbors, and scale is the exponent of the magnitudes for these voltages
    public static double tol = 1e-12; // Relative tolerance for testing whether there is any edge with higher pressure
    //    public static double highTol = 1e-8; // Tolerance for some error checks
    private static boolean verbose = false;

    public int nv; // Number of total vertices
    public int ne; // Number of total edges

    public DirectedEdge[] edges;
    // This will contain the adj list of both incoming and outgoing edges
    // For every vertex, there is a bogus entry in edges in order to mark the end of the list
    public Vertex[] vertices;

    public int[] origNbrsStart;
    // NOTE: nbrsStart will have one additional entry marking the end of the array

    public int numTerminals;

    Random random;

    NodeHeap heap;

    //    private int numActEdges;
    int activeEdgesStart;

    DirectedEdgeGradientHelper directedEdgeGradientHelper;

    // ai contains the *dest* of each edge
    // aj contains the *src* of each edge
    // av contains the lengths of the edges
    public RandLexVoltDirected(int[] ai, int[] aj, double[] av, long randomSeed) {
        // Count number of vertices
        nv = -1;
        for (int i = 0; i < ai.length; i++) {
            if (nv < ai[i])
                nv = ai[i];
            if (nv < aj[i])
                nv = aj[i];
        }
        assert (nv != -1);
        nv = nv + 1;


        //initialize vertices array
        vertices = new Vertex[nv];
        for (int i = 0; i < nv; i++) {
            Vertex v = new Vertex();
            v.actDeg = 0;
            v.isTerminal = false;
            v.isOrigTerminal = false;
            vertices[i] = v;
        }

        // Count number of edges
        ne = 2 * aj.length;
        // We will store two copies of each directed edge, one in the src list, and one in the dest list

        // Initialize arrays
        origNbrsStart = new int[nv + 1];
        edges = new DirectedEdge[ne + nv + 2]; // The last entry will be the tail of the activeEdges list
        // There will be one dummy entry for each vertex marking the start of the adj list for that vertex
        // Note that for legacy reasons, we traverse the adj lists in the reverse order
        for (int i = 0; i < ne + nv + 2; i++) {
            edges[i] = new DirectedEdge(-1, -1, 0, true);
        }

        // Count degree for each vertex
        for (int i = 0; i < aj.length; i++) {
            vertices[ai[i]].actDeg++;
            vertices[aj[i]].actDeg++;
        }

        // Set the pointers for origNbrsStart
        origNbrsStart[0] = 0;
        for (int i = 0; i < nv; i++) {
            origNbrsStart[i + 1] = origNbrsStart[i] + vertices[i].actDeg + 1;
            assert (vertices[i].actDeg > 0);
        }
        assert (origNbrsStart[nv] == ne + nv);

        // Add the edges to the array
        for (int i = 0; i < nv; i++)
            vertices[i].actDeg = 0;

        for (int i = 0; i < aj.length; i++) {
            int dest = ai[i];
            int src = aj[i];
            double len = av[i];

            // Create the two edges
            DirectedEdge e1 = new DirectedEdge(src, dest, len, true); // in src's list
            int loc1 = (origNbrsStart[src] + 1) + vertices[src].actDeg;

            DirectedEdge e2 = new DirectedEdge(dest, src, len, false); // in dest's list
            int loc2 = (origNbrsStart[dest] + 1) + vertices[dest].actDeg;

            // back pointers
            e1.symEdge = loc2;
            e2.symEdge = loc1;

            // Set previous and next active edge
            e1.prevEdge = loc1 - 1;
            e1.nextEdge = -1;
            edges[loc1 - 1].nextEdge = loc1;

            e2.prevEdge = loc2 - 1;
            e2.nextEdge = -1;
            edges[loc2 - 1].nextEdge = loc2;

            // Add the two edges
            edges[loc1] = e1;
            edges[loc2] = e2;

            // Increment degrees
            vertices[src].actDeg++;
            vertices[dest].actDeg++;
        }

        for (int i = 0; i < nv; i++)
            assert (vertices[i].actDeg == origNbrsStart[i + 1] - origNbrsStart[i] - 1);

        numTerminals = 0;

        // Setting up active neighbors data structure
        for (int i = 0; i < nv; i++) {
            Vertex v = vertices[i];
            v.actNbrsEnd = origNbrsStart[i] + vertices[i].actDeg;
            v.deg = v.actDeg;
        }

        // New random number generator
        random = new Random(randomSeed);

        // Initialize data structures for Djikstra
        heap = new NodeHeap(nv);
//
//        // Initialize activeEdges
        edges[ne + nv] = new DirectedEdge(-1, -1, 0, true);
        activeEdgesStart = ne + nv;

        for (int i = 0; i < nv; i++) {
            int eptr = vertices[i].actNbrsEnd;
            vertices[i].actDeg = vertices[i].deg;
            for (int j = 0; j < vertices[i].actDeg; j++) {
                edges[eptr].nextActEdge = activeEdgesStart;
                edges[activeEdgesStart].prevActEdge = eptr;
                edges[eptr].prevActEdge = -1;
                activeEdgesStart = eptr;
                eptr = edges[eptr].prevEdge;
            }
        }


        directedEdgeGradientHelper = new DirectedEdgeGradientHelper(nv, randomSeed);
    }

    public void setTerminals(int[] where, double[] what) {
        for (int i = 0; i < where.length; i++) {
            setVoltage(where[i], what[i]);
            vertices[where[i]].isOrigTerminal = true;
        }
    }

    // NOTE: This removes ONLY active term-term edges
    // There could be non-active term-term edges
    public void setVoltage(int node, double voltage) {
        if (verbose)
            System.out.println("Setting vertex " + node + " with voltage " + voltage);

        Vertex v = vertices[node];
        if (v.isTerminal)
            throw new Error("attempt to set voltage of node already declared a terminal.");

        v.isTerminal = true;
        v.voltage = voltage;

        v.vhigh = voltage;
        v.vlow = voltage;
        v.lowTermDist = 0;
        v.highTermDist = 0;
        v.lowTermVolt = voltage;
        v.highTermVolt = voltage;


        // *******
        // IMP : This is the way to traverse edges at all pressure level
        int eptr = v.actNbrsEnd; // Start from the end of the list
        // NOTE: No -1 since we changed that convention

        int oldActDeg = v.actDeg;
        for (int i = 0; i < oldActDeg; i++) {
            DirectedEdge e1 = edges[eptr];
            Vertex nbr = vertices[e1.dest];
            assert (e1.src == node);
            if (nbr.isTerminal) {
                removeTermTermEdge(eptr);

                // Symmetric op
                removeTermTermEdge(e1.symEdge);
                assert (edges[e1.symEdge].dest == node);
                assert (e1.dest == edges[e1.symEdge].src);

                v.actDeg--;
                nbr.actDeg--;
                // WE REALLY HAVE TO MAKE SURE THE EDGE IS ACTIVE IN BOTH LISTS
            }

            eptr = e1.prevEdge;
        }

        numTerminals++;
    }

    // NOTE: You have to explicitly call removeTermTermEdge for the reverse edge
    private void removeTermTermEdge(int edge) {
        DirectedEdge e = edges[edge];
        Vertex node = vertices[e.src];
        int prev = e.prevEdge;
        int next = e.nextEdge;
        int prevAct = e.prevActEdge;
        int nextAct = e.nextActEdge;

        assert (prev != -1);
        edges[prev].nextEdge = next;
        if (next != -1) {
            edges[next].prevEdge = prev;
        } else {
            assert (edge == node.actNbrsEnd);
            node.actNbrsEnd = prev;
        }

        assert (nextAct != -1);
        edges[nextAct].prevActEdge = prevAct;
        if (prevAct != -1) {
            edges[prevAct].nextActEdge = nextAct;
        } else {
            assert (edge == activeEdgesStart);
            activeEdgesStart = nextAct;
        }

        assert (edge != ne + nv);

        // Decrement degree
        vertices[e.src].deg--;
    }

    // Currently this only returns the gradient of one tight path
    public void computeVoltages() {
        if (verbose)
            System.out.println("computeVoltages");

        activeEdgesStart = ne + nv;
        for (int i = 0; i < nv; i++) {
            int eptr = vertices[i].actNbrsEnd;
            vertices[i].actDeg = vertices[i].deg;
            for (int j = 0; j < vertices[i].actDeg; j++) {
                edges[eptr].nextActEdge = activeEdgesStart;
                edges[activeEdgesStart].prevActEdge = eptr;
                edges[eptr].prevActEdge = -1;
                activeEdgesStart = eptr;

                eptr = edges[eptr].prevEdge;
            }
        }
        setPathsAboveOrEqual2(ne + nv, -1);
        // System.out.println("Found a tight path with gradient : " + tightGrad);
        testFinalVoltages();
    }

    // activeEdgesEnd marks the end of the edges that are active at this pressure.
    // NOTE: The edge pointed to by activeEdgesEnd is NOT active!!!
    private double setPathsAboveOrEqual2(int activeEdgesEnd, double lowergrad) {
        if (verbose)
            System.out.println("setPathsAboveOrEqual2: lowergrad " + lowergrad + " starts");

        double returnGrad = -1;
        while (activeEdgesEnd != activeEdgesStart) {
            if (verbose)
                System.out.println("computeVoltages: lowergrad " + lowergrad + " while");

            zeroActiveDeg(activeEdgesEnd);
            countActiveDeg(activeEdgesEnd);
            int e1 = activeEdgesStart;
            DirectedEdge e1Edge = edges[e1];
            int boundaryEnd = identifyReachableEdges(e1);

            zeroActiveDeg(boundaryEnd);
            countActiveDeg(boundaryEnd);

            int e = sampleActiveEdge(boundaryEnd);
            DirectedEdge eEdge = edges[e];
            double grad = computeEdgeGradient(e);
//            if(!(grad > lowergrad))
//                throw new Error("This shouldn't happen");

            int newActiveEdgesEnd = newIdentifyEdgesStrictlyAbove(boundaryEnd, grad * (1 + tol));

            if (newActiveEdgesEnd == activeEdgesStart) {
                if (verbose)
                    System.out.println("Fixing edges with gradient " + grad);
                if (grad > tol)
                    fixVoltagesWithGradient(boundaryEnd, grad * (1 - tol));
                else {
                    fixVoltagesWithZeroGradient(boundaryEnd);
                }
                returnGrad = grad;
                // Note that inactive t-t edges may not have been removed
            } else {
                if (newActiveEdgesEnd == activeEdgesEnd)
                    System.out.println("Warning: the graph didn't shrink after sampling and picking grad = " + grad);
                returnGrad = setPathsAboveOrEqual2(newActiveEdgesEnd, grad);
            }

            removeAllTermTermEdges(boundaryEnd);

            if ((lowergrad != -1) && (boundaryEnd != activeEdgesStart)) {
                zeroActiveDeg(boundaryEnd);
                countActiveDeg(boundaryEnd);
                // Need to test edges after boundaryEnd and demote some if necessary
                newActiveEdgesEnd = newIdentifyEdgesStrictlyAbove(boundaryEnd, lowergrad * (1 + tol));
//                activeEdgesEnd = newActiveEdgesEnd;

                if (boundaryEnd != newActiveEdgesEnd) {
                    // We have found new inactive edges between boundaryEnd and newActiveEdgesEnd
                    // We need to DEMOTE these edges
                    if (boundaryEnd != activeEdgesEnd) {
                        if (activeEdgesStart == newActiveEdgesEnd) {
                            // No new active edges
                            edges[edges[activeEdgesEnd].prevActEdge].nextActEdge = newActiveEdgesEnd;
                            edges[newActiveEdgesEnd].prevActEdge = edges[activeEdgesEnd].prevActEdge;

                            edges[edges[boundaryEnd].prevActEdge].nextActEdge = activeEdgesEnd;
                            edges[activeEdgesEnd].prevActEdge = edges[boundaryEnd].prevActEdge;

                            activeEdgesStart = boundaryEnd;
                            edges[activeEdgesStart].prevActEdge = -1;
                        } else {
                            // No need to change start
                            edges[edges[activeEdgesEnd].prevActEdge].nextActEdge = newActiveEdgesEnd;
                            edges[edges[newActiveEdgesEnd].prevActEdge].nextActEdge = boundaryEnd;
                            edges[edges[boundaryEnd].prevActEdge].nextActEdge = activeEdgesEnd;

                            int swapedge = edges[activeEdgesEnd].prevActEdge;
                            edges[activeEdgesEnd].prevActEdge = edges[boundaryEnd].prevActEdge;
                            edges[boundaryEnd].prevActEdge = edges[newActiveEdgesEnd].prevActEdge;
                            edges[newActiveEdgesEnd].prevActEdge = swapedge;
                        }
                    }
                    activeEdgesEnd = newActiveEdgesEnd;
                }
            }
        }

        if (verbose)
            System.out.println("computeVoltages: lowergrad " + lowergrad + " ends");

        return returnGrad;
    }

    private void resetVHighVlow(int boundaryEnd) {
        int actPtr = edges[boundaryEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            Vertex v = vertices[e.src];
            if (!v.isTerminal) {
                v.vhigh = Double.NEGATIVE_INFINITY;
                v.vlow = Double.POSITIVE_INFINITY;
            }
            actPtr = edges[actPtr].prevActEdge;
        }
    }

    private void removeAllTermTermEdges(int activeEdgesEnd) {
        int actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            if (vertices[e.src].isTerminal && vertices[e.dest].isTerminal)
                removeTermTermEdge(actPtr);
            actPtr = e.prevActEdge;
        }
    }

    // activeEdgesEnd marks the end of the edges that are active at this pressure.
    // NOTE: The edge pointed to by activeEdgesEnd is NOT active!!!
    private double setPathsAboveOrEqual(int activeEdgesEnd, double lowergrad) {
        zeroActiveDeg(activeEdgesEnd);
        countActiveDeg(activeEdgesEnd);

        int e = sampleActiveEdge(activeEdgesEnd);
        double grad = computeEdgeGradient(e);

        resetVHighVlow(activeEdgesEnd);
        computeVHighVlow(activeEdgesEnd, grad);
        int newActiveEdgesEnd = newIdentifyEdgesStrictlyAbove(activeEdgesEnd, grad);

        double returnGrad;
        if (newActiveEdgesEnd == activeEdgesStart) {
            if (grad > 0)
                fixVoltagesWithGradient(activeEdgesEnd, grad * (1 - tol));
            else
                fixVoltagesWithZeroGradient(activeEdgesEnd);
            returnGrad = grad;
        } else {
            returnGrad = setPathsAboveOrEqual(newActiveEdgesEnd, grad);
            zeroActiveDeg(activeEdgesEnd);
            countActiveDeg(activeEdgesEnd);
        }

        return returnGrad;
    }

    private void fixVoltagesWithZeroGradient(int boundaryEnd) {

        // Since we must have gradient zero, we must assign voltages to non-terminals using paths from terminals

        // We first assign voltages using incoming paths to terminals
        // If there are multiple paths that lead from the non-terminals to terminals, we must pick the one with the lowest voltage
        heap.resetHeap();
        int actPtr = edges[boundaryEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            if (!e.outgoing) {
                int src = edges[actPtr].src;
                int dest = edges[actPtr].dest;
                if (vertices[src].isTerminal) {
                    assert (!vertices[dest].isTerminal);
                    heap.add(dest, vertices[src].voltage);
                }
            }
            actPtr = edges[actPtr].prevActEdge;
        }


        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();
            Vertex v = vertices[node];

            // Go over all active edges
            int eptr = v.actNbrsEnd;
            for (int i = 0; i < v.actDeg; i++) {
                assert (eptr != -1);
                DirectedEdge e = edges[eptr];
                assert (e.src == node);
                if (!e.outgoing) {
                    int nbr = e.dest;
                    if (!vertices[nbr].isTerminal && !heap.finishedThisTime(nbr)) {
                        heap.add(nbr, val);
                    }
                }
                eptr = e.prevEdge;
            }

            setVoltage(node, val);
        }

        // We now assign voltages using outgoing paths from terminals
        // If there are multiple outgoing paths from terminals to one vertex, we must pick the highest voltage
        heap.resetHeap();
        actPtr = edges[boundaryEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            if (e.outgoing) {
                int src = edges[actPtr].src;
                int dest = edges[actPtr].dest;
                if (vertices[src].isTerminal) {
                    assert (!vertices[dest].isTerminal);
                    heap.add(dest, -vertices[src].voltage);
                }
            }
            actPtr = edges[actPtr].prevActEdge;
        }


        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();
            Vertex v = vertices[node];

            // Go over all active edges
            int eptr = v.actNbrsEnd;
            for (int i = 0; i < v.actDeg; i++) {
                assert (eptr != -1);
                DirectedEdge e = edges[eptr];
                assert (e.src == node);
                if (e.outgoing) {
                    int nbr = e.dest;
//                    if(nbr==1073)
//                        System.out.println("Here");
                    if (!vertices[nbr].isTerminal && !heap.finishedThisTime(nbr))
                        heap.add(nbr, val);
                }
                eptr = e.prevEdge;
            }

            setVoltage(node, -val);
        }

    }

    private void zeroActiveDeg(int activeEdgesEnd) {
        int actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            vertices[e.src].actDeg = 0;
            actPtr = edges[actPtr].prevActEdge;
        }
    }

    private void countActiveDeg(int activeEdgesEnd) {
        int actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            vertices[e.src].actDeg++;
            actPtr = edges[actPtr].prevActEdge;
        }
    }

    // This method is identical to computeVHighVlow, except that it tracks lowTermVolt, lowTermDist, highTermVolt, highTermDist
    // Boundary specifies the edges we should initialize Dijkstra using
//    private void fixVoltagesWithGradient(LinkedList<Integer> boundary, double grad) {
    private void fixVoltagesWithGradient(int activeEdgesEnd, double grad) {
        assert (grad > 0);
        double vShift = getMinVoltage(activeEdgesEnd);

        // NOTE: Don't call with zero grad
        // Compute VLow

        // Use the boundary to Initialize the run for VLow
        heap.resetHeap();
        int actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            if (!e.outgoing) {
                int src = edges[actPtr].src;
                int dest = edges[actPtr].dest;
                if (vertices[src].isTerminal) {
                    assert (!vertices[dest].isTerminal);
                    boolean treeModified = heap.add(dest, (vertices[src].voltage - vShift) + grad * e.len);
                    if (treeModified) {
                        vertices[dest].lowTermDist = e.len;
                        vertices[dest].lowTermVolt = vertices[src].voltage - vShift;
                    }
                }
            }
            actPtr = edges[actPtr].prevActEdge;
        }

        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();
            Vertex v = vertices[node];

            v.vlow = val;
//            assert ((v.vlow - v.lowTermVolt - grad * v.lowTermDist) < tol);
//            assert ((v.vlow - v.lowTermVolt - grad * v.lowTermDist) > -tol);
            // Go over all active edges
            int eptr = v.actNbrsEnd;
            for (int i = 0; i < v.actDeg; i++) {
                assert (eptr != -1);
                DirectedEdge e = edges[eptr];
                assert (e.src == node);
                if (!e.outgoing) {
                    int nbr = e.dest;

                    if (!vertices[nbr].isTerminal && !heap.finishedThisTime(nbr)) {
                        boolean treeModified = heap.add(nbr, val + grad * e.len);
                        if (treeModified) {
                            vertices[nbr].lowTermVolt = vertices[node].lowTermVolt;
                            vertices[nbr].lowTermDist = vertices[node].lowTermDist + e.len;
                        }
                    }
                }
                eptr = e.prevEdge;
            }
        }

        // Compute VHigh

        // Use the boundary to Initialize the run for VHigh
        heap.resetHeap();
        actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            if (e.outgoing) {
                int src = edges[actPtr].src;
                int dest = edges[actPtr].dest;
                if (vertices[src].isTerminal) {
                    assert (!vertices[dest].isTerminal);
                    boolean treeModified = heap.add(dest, -(vertices[src].voltage - vShift) + grad * edges[actPtr].len);
                    if (treeModified) {
                        vertices[dest].highTermDist = e.len;
                        vertices[dest].highTermVolt = vertices[src].voltage - vShift;
                    }
                }
            }
            actPtr = edges[actPtr].prevActEdge;
        }

        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();
            Vertex v = vertices[node];

            v.vhigh = -val; // The value being computed is -vhigh

//            if ((v.vhigh - v.highTermVolt + grad * v.highTermDist) > tol) {
//                if (v.isTerminal)
//                    throw new Error("Terminal found!");
//                System.out.println("It seems there is some error in tree computation");
//            }
//            if ((v.vhigh - v.highTermVolt + grad * v.highTermDist) < -tol)
//                System.out.println("It seems there is some error in tree computation");

            // Go over all active edges
            int eptr = v.actNbrsEnd;
            for (int i = 0; i < v.actDeg; i++) {
                assert (eptr != -1);
                DirectedEdge e = edges[eptr];
                assert (e.src == node);
                if (e.outgoing) {
                    int nbr = e.dest;

                    if (!vertices[nbr].isTerminal && !heap.finishedThisTime(nbr)) {
                        boolean treeModified = heap.add(nbr, val + grad * e.len);
                        if (treeModified) {
                            vertices[nbr].highTermVolt = v.highTermVolt;
                            vertices[nbr].highTermDist = v.highTermDist + e.len;
                        }
                    }
                }
                eptr = e.prevEdge;
            }
        }

        // Assign vertices to be terminals if they satisfy requirements
        // TODO: The assignment to voltages happens earlier in the undirected code. Is there an issue?
        actPtr = edges[activeEdgesEnd].prevActEdge;
        boolean anySet = false;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            int node = e.src;
            Vertex v = vertices[node];

            if (!v.isTerminal && (v.vhigh >= v.vlow)) {
                setVoltage(node, vShift + (v.lowTermVolt * v.highTermDist + v.highTermVolt * v.lowTermDist) / (v.lowTermDist + v.highTermDist));
                anySet = true;
            }
            actPtr = edges[actPtr].prevActEdge;
        }
        if (!anySet)
            System.out.println("No terminals set!");

        // Check the assigned terminals
        actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            int node = e.src;
            Vertex v = vertices[node];

            if (v.isTerminal && !v.isOrigTerminal) {
                double inGrad = 0;
                double inVolt = 0;
                double indist = 0;
                double outGrad = 0;
                double outVolt = 0;
                double outdist = 0;

                for (int j = origNbrsStart[node] + 1; j < origNbrsStart[node + 1]; j++) {
                    DirectedEdge e1 = edges[j];
                    Vertex nbr = vertices[e1.dest];
                    if (!nbr.isTerminal)
                        continue;

                    if (e1.outgoing) {
                        double testgrad = (v.voltage - nbr.voltage) / e1.len;
                        if (testgrad > outGrad) {
                            outGrad = testgrad;
                            outVolt = nbr.voltage;
                            outdist = e1.len;
                        }
                    } else {
                        double testgrad = (nbr.voltage - v.voltage) / e1.len;
                        if (testgrad > inGrad) {
                            inGrad = testgrad;
                            inVolt = nbr.voltage;
                            indist = e1.len;
                        }
                    }
                }
//            System.out.println("Vertex " + i + " maxgrad: " + maxgrad + " mingrad: " + mingrad);
                // System.out.println((maxgrad+mingrad) + " " + tol);

//                if (Math.abs(inGrad) + Math.abs(outGrad) < highTol) {
//                    //we won't worry about very small gradients
//                    actPtr = edges[actPtr].prevActEdge;
//                    continue;
//                }
//
//                if ((inGrad - outGrad) > testEps * (inGrad + outGrad) || (inGrad - outGrad) < -testEps * (inGrad + outGrad))
//                    throw new Error("Error when setting gradients! inGrad = " + inGrad + " outGrad =" + outGrad);

/*                double testVolt = (indist * outVolt + outdist * inVolt) / (indist + outdist);
                double scale = Double.longBitsToDouble(Double.doubleToLongBits(Math.abs(testVolt) + Math.abs(v.voltage)) & 0x7ff0000000000000L);
                if ((testVolt - v.voltage) > testEps * scale || (testVolt - v.voltage) < -testEps * scale) {
                    System.out.println("Warning: Error when setting gradients! v.voltage: " + v.voltage + " testvoltage: " + testVolt);
                }*/
            }

            actPtr = edges[actPtr].prevActEdge;
        }
    }

    private int identifyReachableEdges(int edge) {
        int src = edges[edge].src;
        if (vertices[src].isTerminal) {
            src = edges[edge].dest;
        }
        if (vertices[src].isTerminal) {
            throw new Error("Both end points of e are terminals!");
        }

        heap.resetHeap(); // NOTE: We'll use the heap just to track the list of visited vertices
        heap.add(src, 0);
        int newActiveEdgesEnd = activeEdgesStart;
        boolean anyTermsFound = false;

        while (heap.hasMore()) {
            int node = heap.popMin();
            Vertex v = vertices[node];
            int eptr = v.actNbrsEnd;
            DirectedEdge e;

            for (int i = 0; i < v.actDeg; i++) {
                e = edges[eptr];
                int nextEptr = e.prevEdge;

                if (!heap.finishedThisTime(e.dest)) {
                    // move eptr to the head of active edges
                    if (eptr == newActiveEdgesEnd) {
                        newActiveEdgesEnd = edges[eptr].nextActEdge;
                    }
                    if (eptr != activeEdgesStart) {
                        int prevAct = edges[eptr].prevActEdge;
                        int nextAct = edges[eptr].nextActEdge;
                        assert (nextAct != -1);
                        assert (prevAct != -1);

                        edges[nextAct].prevActEdge = edges[eptr].prevActEdge;
                        edges[prevAct].nextActEdge = edges[eptr].nextActEdge;

                        edges[eptr].prevActEdge = -1;
                        edges[eptr].nextActEdge = activeEdgesStart;
                        edges[activeEdgesStart].prevActEdge = eptr;
                        activeEdgesStart = eptr;
                    }

                    // move symedge to the head of active edges
                    if (e.symEdge == newActiveEdgesEnd) {
                        newActiveEdgesEnd = edges[e.symEdge].nextActEdge;
                    }
                    if (e.symEdge != activeEdgesStart) {
                        int prevAct = edges[e.symEdge].prevActEdge;
                        int nextAct = edges[e.symEdge].nextActEdge;
                        assert (nextAct != -1);
                        assert (prevAct != -1);

                        edges[nextAct].prevActEdge = edges[e.symEdge].prevActEdge;
                        edges[prevAct].nextActEdge = edges[e.symEdge].nextActEdge;

                        edges[e.symEdge].prevActEdge = -1;
                        edges[e.symEdge].nextActEdge = activeEdgesStart;
                        edges[activeEdgesStart].prevActEdge = e.symEdge;
                        activeEdgesStart = e.symEdge;
                    }
                }

                assert (edges[eptr].src == node);
                int nbr = edges[eptr].dest;
                if (vertices[nbr].isTerminal) {
                    anyTermsFound = true;
                } else {
                    if (!heap.finishedThisTime(nbr))
                        heap.add(nbr, 0);
                }

                eptr = nextEptr;
            }
        }
        // At the end of this while loop, all the reachable edges have been added to the front of the activeEdges list, all before newActiveEdgesEnd

        if (!anyTermsFound)
            throw new Error("Ran into a component with no terminals!");

        return newActiveEdgesEnd;
    }


    private double getMinVoltage(int boundaryEnd) {
        double minV = Double.POSITIVE_INFINITY;

        int actPtr = edges[boundaryEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            Vertex v = vertices[e.src];
            if (v.isTerminal) {
                double volt = v.voltage;
                if (volt < minV)
                    minV = volt;
            }
            actPtr = edges[actPtr].prevActEdge;
        }

        return minV;
    }

    private void applyForwardTransformation(double shift, int boundaryEnd) {
        heap.resetHeap();
        // We are using the heap just to store the terminals in the current component

        int actPtr = edges[boundaryEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            Vertex v = vertices[e.src];
            if (v.isTerminal)
                heap.add(e.src, 0);
            actPtr = edges[actPtr].prevActEdge;
        }

        // Apply transformation to all inputs
        while (heap.hasMore()) {
            int node = heap.popMin();
            vertices[node].voltage += shift;
            vertices[node].vhigh += shift;
            vertices[node].vlow += shift;
        }
    }


    private void applyReverseTransformation(double shift, int boundaryEnd) {
        heap.resetHeap();
        // We are using the heap just to store the terminals in the current component

        int actPtr = edges[boundaryEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            Vertex v = vertices[e.src];
            if (v.isTerminal)
                heap.add(e.src, 0);
            actPtr = edges[actPtr].prevActEdge;
        }

        // Apply transformation to all inputs
        while (heap.hasMore()) {
            int node = heap.popMin();
            vertices[node].voltage -= shift;
            vertices[node].vhigh -= shift;
            vertices[node].vlow -= shift;
        }
    }

    private int newIdentifyEdgesStrictlyAbove(int activeEdgesEnd, double grad) {
        double vshift = getMinVoltage(activeEdgesEnd);
        applyReverseTransformation(vshift, activeEdgesEnd);
        computeVHighVlow(activeEdgesEnd, grad);

        int newActiveEdgesEnd = activeEdgesStart;
        int actPtr = activeEdgesStart;

        while (actPtr != activeEdgesEnd) {
            DirectedEdge e = edges[actPtr];
            int nextActPtr = edges[actPtr].nextActEdge;

            // Test if e is above grad
            boolean flag = false;
            if (e.outgoing) {
                if ((vertices[e.src].vhigh - vertices[e.dest].vlow) > e.len * grad)
                    flag = true;
            } else {
                if ((vertices[e.dest].vhigh - vertices[e.src].vlow) > e.len * grad)
                    flag = true;
            }
            if (flag) {
                // Move e to front of actEdges list
                assert (edges[actPtr].nextActEdge != -1);
                if (actPtr == newActiveEdgesEnd) {
                    newActiveEdgesEnd = nextActPtr;
                }
                if (actPtr != activeEdgesStart) {
                    int prevAct = edges[actPtr].prevActEdge;
                    int nextAct = edges[actPtr].nextActEdge;
                    assert (prevAct != -1);

                    edges[nextAct].prevActEdge = edges[actPtr].prevActEdge;
                    edges[prevAct].nextActEdge = edges[actPtr].nextActEdge;

                    edges[actPtr].prevActEdge = -1;
                    edges[actPtr].nextActEdge = activeEdgesStart;
                    edges[activeEdgesStart].prevActEdge = actPtr;
                    activeEdgesStart = actPtr;
                }

                // Move e to end of its source's adjacency list
                int prev = e.prevEdge;
                int next = e.nextEdge;
                Vertex v = vertices[e.src];
                if (actPtr != v.actNbrsEnd) {
                    // You have to move eptr only if it is not at actNbrsEnd
                    assert (next != -1);
                    edges[next].prevEdge = e.prevEdge;

                    assert (prev != -1);
                    edges[prev].nextEdge = e.nextEdge;

                    // Add e to the head of the list
                    e.prevEdge = v.actNbrsEnd;
                    e.nextEdge = -1;
                    edges[v.actNbrsEnd].nextEdge = actPtr;
                    v.actNbrsEnd = actPtr;
                }
            }

            actPtr = nextActPtr;
        }

        applyForwardTransformation(vshift, activeEdgesEnd);
        return newActiveEdgesEnd;
    }

    // Boundary specifies the edges we should initialize Dijkstra using
    private void computeVHighVlow(int activeEdgesEnd, double grad) {
        resetVHighVlow(activeEdgesEnd);

        // Compute VLow

        // Use the boundary to Initialize the run for VLow
        heap.resetHeap();
        int actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            if (!edges[actPtr].outgoing) {
                // For vlow, we only need incoming edges into terminals for initialization
                int src = edges[actPtr].src;
                int dest = edges[actPtr].dest;
                if (vertices[src].isTerminal) {
                    if (vertices[dest].isTerminal)
                        throw new Error("An active edge has both end points as terminals!");
                    heap.add(dest, vertices[src].voltage + grad * edges[actPtr].len);
                }
            }
            actPtr = edges[actPtr].prevActEdge;
        }

        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();
            Vertex v = vertices[node];

            v.vlow = val;
            // Go over all active edges
            int eptr = v.actNbrsEnd;
            for (int i = 0; i < v.actDeg; i++) {
                assert (eptr != -1);
                DirectedEdge e = edges[eptr];
                assert (e.src == node);
                int nbr = e.dest;

                if (!e.outgoing && !vertices[nbr].isTerminal && !heap.finishedThisTime(nbr)) {
                    heap.add(nbr, val + grad * e.len);
                }
                eptr = e.prevEdge;
            }
        }

        // Compute VHigh

        // Use the boundary to Initialize the run for VHigh
        heap.resetHeap();
        actPtr = edges[activeEdgesEnd].prevActEdge;
        while (actPtr != -1) {
            DirectedEdge e = edges[actPtr];
            if (e.outgoing) {
                // Only outgoing edges for vhigh
                int src = e.src;
                int dest = e.dest;
                if (vertices[src].isTerminal) {
                    assert (!vertices[dest].isTerminal);
                    heap.add(dest, -vertices[src].voltage + grad * edges[actPtr].len);
                }
            }
            actPtr = edges[actPtr].prevActEdge;
        }

        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();
            Vertex v = vertices[node];

            v.vhigh = -val; // The value being computed is -vhigh

            // Go over all active edges
            int eptr = v.actNbrsEnd;
            for (int i = 0; i < v.actDeg; i++) {
                assert (eptr != -1);
                DirectedEdge e = edges[eptr];
                assert (e.src == node);
                int nbr = e.dest;

                if (e.outgoing && !vertices[nbr].isTerminal && !heap.finishedThisTime(nbr)) {
                    heap.add(nbr, val + grad * e.len);
                }
                eptr = e.prevEdge;
            }
        }
    }

    private int sampleActiveEdge(int activeEdgesEnd) {
        int numActEdges = 0;
        int eptr = activeEdgesEnd;
        while (eptr != activeEdgesStart) {
            numActEdges++;
            eptr = edges[eptr].prevActEdge;
        }
        assert (numActEdges % 2 == 0);

        int randomIndex = random.nextInt(numActEdges);
        eptr = activeEdgesStart;
        for (int i = 0; i < randomIndex; i++) {
            eptr = edges[eptr].nextActEdge;
        }
        return eptr;
    }

    public double computeEdgeGradient(int e) {
        int src = edges[e].src;
        int dest = edges[e].dest;
        if (vertices[src].isTerminal && vertices[dest].isTerminal)
            throw new Error("Why are we computing the gradient of a term-term edge");
        //we should not compute the gradient of a term-term edge

        if (edges[e].outgoing) {
            // Perform Djikstra from the src of the edge to find src terms
            findSrcTermsFromVertex(src);

            // Perform Djikstra from the dest of the edge to find dest terms
            findDestTermsFromVertex(dest);
        } else {
            findSrcTermsFromVertex(dest);
            findDestTermsFromVertex(src);
        }

        double grad;
        grad = directedEdgeGradientHelper.computeRandGrad(edges[e].len);

        return grad;
    }

    // Performs Djikstra from the node v
    // Only taking incoming edges
    // Adds the terminals discovered to edgeGradientHelper as left terminals
    private void findSrcTermsFromVertex(int v) {
        if (vertices[v].isTerminal) {
            directedEdgeGradientHelper.addSrc(vertices[v].voltage, 0.0);
        } else {
            // Perform Djikstra only if v is not a terminal
            heap.resetHeap();
            heap.add(v, 0);

            while (heap.hasMore()) {
                double val = heap.minKey();
                int node = heap.popMin();

                if (vertices[node].isTerminal) {
                    directedEdgeGradientHelper.addSrc(vertices[node].voltage, val);
                } else {
                    int eptr = vertices[node].actNbrsEnd;
                    for (int i = 0; i < vertices[node].actDeg; i++) {
                        if (!edges[eptr].outgoing) {
                            // Take only incoming edges
                            int nbr = edges[eptr].dest;
                            if (!heap.finishedThisTime(nbr)) { // Only if the vertex isn't already finished in this iteration of Djikstra
                                heap.add(nbr, val + edges[eptr].len); // NOTE: Even if the node already exists, the heap handles it
                            }
                        }
                        eptr = edges[eptr].prevEdge;
                    }
                }
            }
        }
    }

    // Performs Djikstra from the node v
    // Only taking outgoing edges
    // Adds the terminals to directedEdgeGradientHelper as right terminals
    private void findDestTermsFromVertex(int v) {
        if (vertices[v].isTerminal) {
            directedEdgeGradientHelper.addDest(vertices[v].voltage, 0.0);
            return;
        }
        // Perform Djikstra only if v is not a terminal

        heap.resetHeap();
        heap.add(v, 0);

        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();

            if (vertices[node].isTerminal) {
                directedEdgeGradientHelper.addDest(vertices[node].voltage, val);
            } else {
                int eptr = vertices[node].actNbrsEnd;
                for (int i = 0; i < vertices[node].actDeg; i++) {
                    if (edges[eptr].outgoing) {
                        // Only if it's an outgoing edge
                        int nbr = edges[eptr].dest;
                        if (!heap.finishedThisTime(nbr)) { // Only if the vertex isn't already finished in this iteration of Djikstra
                            heap.add(nbr, val + edges[eptr].len); // NOTE: Even if the node already exists, the heap handles it
                        }
                    }
                    eptr = edges[eptr].prevEdge;
                }
            }
        }
    }

    // Code to test final voltages
    public void testFinalVoltages() {
        if (numTerminals < nv)
            throw new Error("Test only after all voltages have been computed");
        for (int i = 0; i < nv; i++) {
            Vertex v = vertices[i];
            if (v.isOrigTerminal)
                continue;
            double ingrad = 0;
            double outgrad = 0;
            double inVolt = v.voltage;
            double inDist = 0;
            double outVolt = v.voltage;
            double outDist = 0;

            for (int j = origNbrsStart[i] + 1; j < origNbrsStart[i + 1]; j++) {
                DirectedEdge e = edges[j];
                Vertex nbr = vertices[e.dest];
                if (e.outgoing) {
                    // Outgoing edge
                    double grad = (v.voltage - nbr.voltage) / e.len;
                    if (grad > outgrad) {
                        outgrad = grad;
                        outVolt = nbr.voltage;
                        outDist = e.len;
                    }
                } else {
                    // incoming edge
                    double grad = (nbr.voltage - v.voltage) / e.len;
                    if (grad > ingrad) {
                        ingrad = grad;
                        inVolt = nbr.voltage;
                        inDist = e.len;
                    }
                }
            }

            double testVolt = (outDist * inVolt + inDist * outVolt) / (outDist + inDist);
            double scale = Double.longBitsToDouble(Double.doubleToLongBits(Math.abs(testVolt) + Math.abs(v.voltage)) & 0x7ff0000000000000L);
//            double scale = Math.pow(2, StrictMath.getExponent(Math.abs(testVolt) + Math.abs(v.voltage)));
/*            if ((testVolt - v.voltage) > testEps * scale || (testVolt - v.voltage) < -testEps * scale) {
                System.out.println("Warning: There seems to be some error in the voltages v.voltage: " + v.voltage + " testvoltage: " + testVolt);
            }*/
//            if (Math.abs(ingrad) + Math.abs(outgrad) < highTol) {
//                //we won't worry about very small gradients
//                continue;
//            }
//
//            if ((outgrad - ingrad) > testEps * (ingrad + outgrad) || (outgrad - ingrad) < -testEps * (ingrad + outgrad)) {
//                throw new Error("There seems to be some error in final gradients! inGrad = " + ingrad + " outgrad =" + outgrad);
//            }
        }
        /*System.out.println("Tested final voltages");*/
    }
}
