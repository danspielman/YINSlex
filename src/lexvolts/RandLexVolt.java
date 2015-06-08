package lexvolts;

/**
 * Created by spielman on 9/1/14.
 *
 * Uses EdgeListGraph, and Sushant's randomized ideas to make a fast algorithm
 */

import java.util.*;

public class RandLexVolt {

    public static double testEps = 1e-10; // We test that the maxgrad + mingrad at each node is at most testEps*(maxgrad - mingrad)

    // will modify this
    public EdgeListGraph g;

    public int[] terminalList;
    public boolean[] isTerminal;
    public boolean[] isOrigTerminal;
    public double[] voltages;
    public int numTerminals;
    public int numOrigTerminals;
    double[] vlow;
    double[] vhigh;
    Random random;

    private ReusableNodeHeap heap;

    private EdgeGradientHelper edgeGradientHelper;

    protected static int threadThreshold = 10;

//    public static Logger logger = new Logger();
    //  logger.start("rlv.log");

    public RandLexVolt(RandLexVolt orig) {
        this.g = new EdgeListGraph(orig.g);
        this.terminalList = orig.terminalList.clone();
        this.isTerminal = orig.isTerminal.clone();
        this.voltages = orig.voltages.clone();
        this.numTerminals = orig.numTerminals;
        this.numOrigTerminals = orig.numOrigTerminals;
        this.isOrigTerminal = orig.isOrigTerminal.clone();
//        this.vlow = orig.vlow.clone();
//        this.vhigh = orig.vhigh.clone();
        this.random = new Random(orig.random.nextLong());

        this.heap = new ReusableNodeHeap(this.g.nv);
        this.edgeGradientHelper = new EdgeGradientHelper(orig.g.nv, orig.random.nextLong());
    }

    // makes a copy of the input graph, because it modifies it
    public RandLexVolt(EdgeListGraph inputg) {
        long randomSeed = 0;
        g = new EdgeListGraph(inputg);

        terminalList = new int[g.nv];
        isTerminal = new boolean[g.nv];
        isOrigTerminal = new boolean[g.nv];
        for (int i = 0; i < g.nv; i++) {
            isTerminal[i] = false;
            isOrigTerminal[i] = false;
        }

        voltages = new double[g.nv];
        numTerminals = 0;

        random = new Random(randomSeed);
        // Specify the seed here

        this.heap = new ReusableNodeHeap(this.g.nv);
        this.edgeGradientHelper = new EdgeGradientHelper(this.g.nv, this.random.nextLong());
    }


    // makes a copy of the input graph, because it modifies it
    public RandLexVolt(EdgeListGraph inputg, long randomSeed) {
        g = new EdgeListGraph(inputg);

        terminalList = new int[g.nv];
        isTerminal = new boolean[g.nv];
        isOrigTerminal = new boolean[g.nv];
        for (int i = 0; i < g.nv; i++) {
            isTerminal[i] = false;
            isOrigTerminal[i] = false;
        }

        voltages = new double[g.nv];
        numTerminals = 0;

        random = new Random(randomSeed);
        // Specify the seed here

        this.heap = new ReusableNodeHeap(this.g.nv);
        this.edgeGradientHelper = new EdgeGradientHelper(this.g.nv, this.random.nextLong());
    }

    public void setTerminals(int[] where, double[] what) {
        for (int i = 0; i < where.length; i++) {
            setVoltage(where[i], what[i]);
            isOrigTerminal[where[i]] = true;
        }
        numOrigTerminals = numTerminals;
    }

    // make node a terminal, set its voltage, and remove its terminal nbrs from the graph
    // The terminal neighbors are moved to the end of the adjacency list for that vertex, and the degree of the vertex is reduced.
    // this is done inefficiently right now.  It could be sped up with back pointers for edges.
    public void setVoltage(int node, double voltage) {
        if (isTerminal[node])
            throw new Error("attempt to set voltage of node already declared a terminal.");

        isTerminal[node] = true;
        terminalList[numTerminals++] = node;
        voltages[node] = voltage;

        for (int i = g.nbrsStart[node]; i < g.nbrsStart[node] + g.deg[node]; i++) {
            int nbr = g.ai[i];
            if (isTerminal[nbr]) {
                int j = g.nbrsStart[node] + g.deg[node] - 1;
                int swap = g.ai[j];
                g.ai[j] = nbr;
                g.ai[i] = swap;

                double tmp = g.av[j];
                g.av[j] = g.av[i];
                g.av[i] = tmp;

                g.deg[node]--;
                i = i - 1;

                // unfortunately, we now have to do the symmetric op, which is slow if the vertex
                // has high degree in this graph data structure

                for (int k = g.nbrsStart[nbr]; k < g.nbrsStart[nbr] + g.deg[nbr]; k++) {
                    if (g.ai[k] == node) {
                        int k2 = g.nbrsStart[nbr] + g.deg[nbr] - 1;
                        g.ai[k] = g.ai[k2];
                        g.ai[k2] = node;

                        tmp = g.av[k];
                        g.av[k] = g.av[k2];
                        g.av[k2] = tmp;

                        g.deg[nbr]--;

                    }
                }
            }
        }
    }


    // Returns the value of the steepest gradient in the instance
    public double steepestGradient() {
        if (numTerminals == 0)
            throw new Error("Computing gradient in a zero terminal instance!");

        if (numTerminals == 1)
            return 0.0;

        int term = randActiveTerminal();
        double grad = gradientFromTerminal(term);
        System.out.println("Steepest gradient from one terminal : " + grad);

        boolean[] flag = edgesStrictlyAbovePressure(grad);

        boolean anyAbove = false;
        for (int i = 0; i < flag.length; i++) {
            if (flag[i]) {
                anyAbove = true;
                break;
            }
        }

        if (!anyAbove) {
            System.out.println("Nothing else left! Double checking");
            double grad2 = gradientFromTerminal(term);
            if (grad != grad2) {
                throw new Error("Gradients no longer match!");
            } else {
                System.out.println("Gradient checked! val : " + grad2);
            }
            return grad;
        }
        else {
            Vector<Subgraph> subs = edgeInducedComponentSubgraphs(flag);
            System.out.println("Calling subgraphs");
            double maxgrad = grad;
            for (Subgraph sub : subs) {
                double candidateGrad = sub.rlv.steepestGradient();
                System.out.println("Subgraph returned gradient " + candidateGrad);
                if (candidateGrad > maxgrad)
                    maxgrad = candidateGrad;
            }
            if (maxgrad == grad) {
                throw new Error("Subgraphs returned lower gradient");
            }
            return maxgrad;
        }

    }

    // Sushant: This method uses the new way to assign voltages to vertices on tight paths
    public void computeVoltages5() {
        boolean connected = isGraphConnected();

        boolean[] flag = new boolean[g.ne]; // flag false for term-term edges
        for (int i = 0; i < g.ne; i++) {
            if (!isTerminal[g.ai[i]] || !isTerminal[g.aj[i]])
                flag[i] = true;
        }

        Vector<Subgraph> subs = edgeInducedComponentSubgraphs(flag);
        for (Subgraph sub : subs) {
            if (sub.rlv.numTerminals == 0)
                throw new Error("Graph has a component with 0 terminals!");
            // System.out.println("Calling subgraph with " + sub.rlv.numTerminals + " terminals");
            sub.rlv.workOnSubgraphComponentsWithNewFixMethod(0.0);
            if (sub.rlv.numTerminals < sub.rlv.g.nv) {
                throw new Error("recursion failed. Lowergradient :" + "0.0");
            }
            setTerminalsFromSubgraph(sub);
        }
     /*   testFinalVoltages();
        testFinalVoltagesScaled();*/

//        int term = randActiveTerminal();
//        System.out.println("Has active terminals initially");
        if (!connected)
            System.out.println("Original graph disconnected");

        if (numTerminals < g.nv) {
            throw new Error("Unfinished work!");
//            workOnSubgraphComponents(0.0);
        }
    }


    // Works on this graph, and assigns path until the gradient falls below lowergrad
    // Assumes that there actually exists a path steeper than lowergrad
    /* a more agressive version of workOnSubgraph
     and keeps working as long as it has edges above the lowergrad value supplied */
    // This procedure must be called without terminal-terminal edges
    public void workOnSubgraphComponentsWithNewFixMethod(double lowergrad) {

        if (numTerminals == 0)
            throw new Error("strange! Number of terminals " + numTerminals);

        if (numTerminals == 1) {
            setNodesFromTerminalWithZeroGradient(terminalList[0]);
            return;
        }

        int edge = sampleActiveEdge();
        double grad = computeEdgeGradient(edge);

        if (grad == 0) {
            int term = randActiveTerminal();
            setNodesFromTerminalWithZeroGradient(term);
//            System.out.println("Fixed voltages for " + (numTerminals - numOldTerms) + " terminals");
        } else {
            // Non-zero gradient

            boolean[] flag = edgesStrictlyAbovePressure(grad);
            //      System.out.println("Called edgesAbovePressure");

            // If none strictly above pressure, then set with this.
            boolean anyAbove = false;
            boolean allAbove = true;
            for (int i = 0; i < flag.length; i++) {
                if (flag[i]) anyAbove = true;
                allAbove = flag[i] && allAbove;
            }

            if (allAbove) {
                throw new Error("recursion failed: did not reduce graph size");
            }

            if (!anyAbove) {
                //System.out.println("Setting with gradient " + grad);
                    setVoltagesWithGrad(grad);
            } else {

                Vector<Subgraph> subs = edgeInducedComponentSubgraphs(flag);
                for (Subgraph sub : subs) {
                    // System.out.println("Calling subgraph with " + sub.rlv.numTerminals + " terminals");
                    sub.rlv.workOnSubgraphComponentsWithNewFixMethod(grad);

                    if (sub.rlv.numTerminals <= sub.rlv.numOrigTerminals) {
                        throw new Error("recursion failed. Lowergradient :" + lowergrad);
                    }

                    setTerminalsFromSubgraph(sub);
                }

            }
        }

        boolean[] flag2 = edgesAboveOrEqualPressure(lowergrad);
        boolean anyAbove = false;
        for (int i = 0; i < flag2.length; i++)
            if (flag2[i]) anyAbove = true;

        if (anyAbove) {
//            System.out.println("Recursing on graph again with " + numTerminals + " terminals");
            Subgraph sub = edgeInducedSubgraph(flag2);
            sub.rlv.workOnSubgraphComponentsWithNewFixMethod(lowergrad);
            setTerminalsFromSubgraph(sub);
        }

    }


    public WorkOnSubgraphComponentsThread makework(double lowergrad, Subgraph sub) {
        return new WorkOnSubgraphComponentsThread(lowergrad, sub);
    }

    // a copy of workOnSubgraphComponentsWithNewFixMethod,
    // designed for threading

    public class WorkOnSubgraphComponentsThread extends java.util.concurrent.RecursiveTask<Void> {


        double lowergrad;
        Subgraph subgraph;

        public WorkOnSubgraphComponentsThread(double lowergrad, Subgraph subgraph) {
            this.lowergrad = lowergrad;
            this.subgraph = subgraph;
        }

        protected Void compute() {

            RandLexVolt rl = this.subgraph.rlv;

            if (rl.g.ne < threadThreshold) {
                rl.workOnSubgraphComponentsWithNewFixMethod(lowergrad);

            } else {

                if (rl.numTerminals == 0)
                    throw new Error("strange! Number of terminals " + rl.numTerminals);

                if (rl.numTerminals == 1) {
                    rl.setNodesFromTerminalWithZeroGradient(rl.terminalList[0]);
                    return null;
                }

                int edge = rl.sampleActiveEdge();
                double grad = rl.computeEdgeGradient(edge);

                if (grad == 0) {
                    int term = rl.randActiveTerminal();
                    rl.setNodesFromTerminalWithZeroGradient(term);
                    //            System.out.println("Fixed voltages for " + (numTerminals - numOldTerms) + " terminals");
                } else {
                    // Non-zero gradient

                    boolean[] flag = rl.edgesStrictlyAbovePressure(grad);
                    //      System.out.println("Called edgesAbovePressure");

                    // If none strictly above pressure, then set with this.
                    boolean anyAbove = false;
                    boolean allAbove = true;
                    for (int i = 0; i < flag.length; i++) {
                        if (flag[i]) anyAbove = true;
                        allAbove = flag[i] && allAbove;
                    }

                    if (allAbove) {
                        throw new Error("recursion failed: did not reduce graph size");
                    }

                    if (!anyAbove) {
                        //System.out.println("Setting with gradient " + grad);
                        rl.setVoltagesWithGrad(grad);
                    } else {

                        Vector<Subgraph> subs = rl.edgeInducedComponentSubgraphs(flag);

 /*                       if ((subs.size() > 1) && (rl.g.ne > 100)) {
                            System.out.print("split : ");
                            for (Subgraph sub : subs)
                                System.out.print(sub.rlv.g.ne + " ");
                            System.out.print("\n");
                        }*/


                        List<WorkOnSubgraphComponentsThread> subtasks =
                                new ArrayList<WorkOnSubgraphComponentsThread>();

                        for (Subgraph sub : subs) {
//                            if ((sub.rlv.numTerminals > 100)) {
//                                System.out.println("Calling subgraph with " + sub.rlv.numTerminals + " terminals");
//                            }

                            if (sub.rlv.g.ne < threadThreshold) {
                                sub.rlv.workOnSubgraphComponentsWithNewFixMethod(grad);
                            } else {
                                WorkOnSubgraphComponentsThread task = new WorkOnSubgraphComponentsThread(grad, sub);
                                subtasks.add(task);
                                task.fork();
                            }


                        }

                        for (WorkOnSubgraphComponentsThread task : subtasks) {
                            task.join();


                        }

                        for (Subgraph sub : subs) {
                            if (sub.rlv.numTerminals <= sub.rlv.numOrigTerminals) {
                                throw new Error("recursion failed. Lowergradient :" + lowergrad);
                            }

                            rl.setTerminalsFromSubgraph(sub);
                        }


                    }
                }


                boolean[] flag2 = rl.edgesAboveOrEqualPressure(lowergrad);
                boolean anyAbove = false;
                for (int i = 0; i < flag2.length; i++)
                    if (flag2[i]) anyAbove = true;

                if (anyAbove) {
//            System.out.println("Recursing on graph again with " + numTerminals + " terminals");
                    Subgraph sub = rl.edgeInducedSubgraph(flag2);


                    if (sub.rlv.g.ne < threadThreshold) {
                        sub.rlv.workOnSubgraphComponentsWithNewFixMethod(lowergrad);
                    } else {
                        WorkOnSubgraphComponentsThread task = new WorkOnSubgraphComponentsThread(lowergrad, sub);
                        task.fork();
                        task.join();
                    }
                    rl.setTerminalsFromSubgraph(sub);
                }

            }
            return null;
        }

    }



    Vector<Subgraph> edgeInducedComponentSubgraphs(boolean[] flag) {
        int[] componentIndex = identifyComponents(flag);

        int largestComponentIndex = -1;
        for (int i = 0; i < g.ne; i++)
            if (largestComponentIndex < componentIndex[i])
                largestComponentIndex = componentIndex[i];
        int numberOfComponents = largestComponentIndex + 1;

        Vector<Subgraph> subgraphs = new Vector<Subgraph>();

        // STARTING CODE TO IDENTIFY SUBGRAPHS IN LINEAR TIME

        // First create a vector of lists of the edges for each of the components
        Vector<LinkedList<Integer>> edgeLists = new Vector<LinkedList<Integer>>();
        for (int componentNumber = 0; componentNumber < numberOfComponents; componentNumber++) {
            edgeLists.add(new LinkedList<Integer>());
        }

        for (int i = 0; i < g.ne; i++) {
            if (flag[i]) {
                // Add edge i to list of edges for componentIndex[i]
                if (isTerminal[g.ai[i]] && isTerminal[g.aj[i]])
                    throw new Error("Getting terminal-terminal edges");
                edgeLists.elementAt(componentIndex[i]).add(i);
            }
        }


        for (int componentNumber = 0; componentNumber < numberOfComponents; componentNumber++) {


            /// BEGIN PASTED CODE

            // this maps nodes in the orig to nodes in the subgraph
            /* int[] nodeProj = new int[g.nv];
            for (int i = 0; i < g.nv; i++)
                nodeProj[i] = -1; */
            HashMap<Integer, Integer> nodeProj = new HashMap<Integer, Integer>();
            LinkedList<Integer> vertexList = new LinkedList<Integer>();

            int nvSub = 0;
            int neSub = 0;
            int prevNode = -1;
            int numSubTerminals = 0;

            // find the nodes that appear, and name them
            Iterator<Integer> edgeIter = edgeLists.elementAt(componentNumber).iterator();
            while (edgeIter.hasNext()) {
                int edge = edgeIter.next();
                neSub++;
                int node = g.aj[edge];
                if (node != prevNode) {
                    assert (node > prevNode);
                    prevNode = node;
                    vertexList.add(node);
                    //nodeProj[node] = nvSub;
                    nodeProj.put(node, nvSub);
                    nvSub++;
                    if (isTerminal[node])
                        numSubTerminals++;
                }
            }


            int[] mapBack = new int[nvSub];

            int[] subTerminals = new int[numSubTerminals];
            double[] subVoltages = new double[numSubTerminals];

            int subTerminalInd = 0;

            Iterator<Integer> vertexIter = vertexList.iterator();
            int nodeIndex = -1;
            while (vertexIter.hasNext()) {
                int node = vertexIter.next();
                nodeIndex++;
                mapBack[nodeIndex] = node;
                if (isTerminal[node]) {
                    subVoltages[subTerminalInd] = voltages[node];
                    subTerminals[subTerminalInd++] = nodeIndex;
                }
            }

            // now shift voltages to improve accuracy
            double shift = getMinVoltage(subVoltages);
            for(int i = 0; i<subVoltages.length; i++) {
                subVoltages[i] = subVoltages[i]-shift;
            }

            int[] subi = new int[neSub];
            int[] subj = new int[neSub];
            double[] subv = new double[neSub];

            int ePtr = 0;

            edgeIter = edgeLists.elementAt(componentNumber).iterator();

            while (edgeIter.hasNext()) {
                int edge = edgeIter.next();
                subi[ePtr] = nodeProj.get(g.ai[edge]);
                subj[ePtr] = nodeProj.get(g.aj[edge]);
                subv[ePtr] = g.av[edge];
                ePtr++;
            }

            EdgeListGraph subg = new EdgeListGraph(subi, subj, subv);
            RandLexVolt rlv = new RandLexVolt(subg, random.nextLong());
            rlv.setTerminals(subTerminals, subVoltages);

            Subgraph sg = new Subgraph();
            sg.mapBack = mapBack;
            sg.rlv = rlv;
            sg.shift = shift;


            subgraphs.add(sg);
            // END OF PASTED CODE


        }

        /// End of the code modification

        /*

        for(int componentNumber=0;componentNumber < numberOfComponents; componentNumber++){
            boolean[] componentFlag = new boolean[g.ne];
            boolean nonEmptySubgraph = false;

            for(int i=0;i<g.ne;i++) {
                componentFlag[i] = false;
            }
            Iterator<Integer> edgeIter = edgeLists.elementAt(componentNumber).iterator();
            while(edgeIter.hasNext()) {
                int edge = edgeIter.next();
                componentFlag[edge] = true;
                nonEmptySubgraph = true;
            }
            if(nonEmptySubgraph)
                subgraphs.add(edgeInducedSubgraph(componentFlag));
        } */

//        logger.write("Number of components : " + subgraphs.size());

        return subgraphs;
    }

    // Returns a vector of integers, one for each edge, identifying the component of the edge
    int[] identifyComponents(boolean[] flag) {
        int[] vertexComponentIndex = new int[g.nv];
        boolean[] completed = new boolean[g.nv];
        // completed indicates if you've already processed the vertex

        for (int i = 0; i < g.nv; i++) {
            vertexComponentIndex[i] = -1;
            completed[i] = false;
        }

        // Push all non-terminals connected to "active" edges (with flag[e] = true) on to a stack
        int componentNumber = -1;
        Stack<Integer> nodeStack = new Stack<Integer>();
        for (int i = 0; i < g.ne; i++) {
            if (flag[i] && !isTerminal[g.aj[i]]) {
                nodeStack.push(g.aj[i]);
            }
            if (flag[i] && !isTerminal[g.ai[i]]) {
                nodeStack.push(g.ai[i]);
            }
        }

        while (!nodeStack.isEmpty()) {
            int node = nodeStack.pop();
            if (!completed[node]) {
                if (vertexComponentIndex[node] == -1) {
                    componentNumber++;
                    vertexComponentIndex[node] = componentNumber;
                }
                for (int i = g.nbrsStart[node]; i < g.nbrsStart[node] + g.deg[node]; i++) {
                    int nbr = g.ai[i];
                    if (flag[i] && (!isTerminal[nbr])) {
                        if (vertexComponentIndex[nbr] != -1 && vertexComponentIndex[nbr] != vertexComponentIndex[node])
                            throw new Error("Sounds like an error");
                        if ((vertexComponentIndex[nbr] == -1)) {
//                          if (flag[i] && (vertexComponentIndex[nbr] == -1)) {
                            // Edge is in the graph, and the neighbor hasn't been seen so far
                            vertexComponentIndex[nbr] = vertexComponentIndex[node];
                            nodeStack.push(nbr);
                        }
                    }
                }

                completed[node] = true;
            }
        }

        // Now vertexComponentIndex contains the index number of the component for each of the vertices
        // However, we do not assign component numbers to terminals. Moreover, components are not allowed to cross terminals

        int[] edgeComponentIndex = new int[g.ne];
        for (int i = 0; i < g.ne; i++) {
            if (flag[i]) {
                if (!(((vertexComponentIndex[g.ai[i]] == -1) && isTerminal[g.ai[i]]) || ((vertexComponentIndex[g.aj[i]] == -1) && isTerminal[g.aj[i]]) || (vertexComponentIndex[g.ai[i]] == vertexComponentIndex[g.aj[i]])))
                    throw new Error("Incorrectly identified components!");
                edgeComponentIndex[i] = vertexComponentIndex[g.ai[i]];
                if (edgeComponentIndex[i] == -1)
                    edgeComponentIndex[i] = vertexComponentIndex[g.aj[i]];
                if (edgeComponentIndex[i] == -1)
                    throw new Error("Active edge between two terminals??");
            }
        }
/*
        for(int i=0;i<g.ne;i++){
            if(flag[i]){
                componentIndex[g.ai[i]]=0;
                componentIndex[g.aj[i]]=0;
            }
        } */
        return edgeComponentIndex;
    }


    void setNodesFromTerminalWithZeroGradient(int term) {
        Stack<Integer> nodeStack = new Stack<Integer>();
        nodeStack.push(term);
        double volt = voltages[term];

        while (!nodeStack.isEmpty()) {
            int node = nodeStack.pop();
            if (isTerminal[node]) {
                if (voltages[node] != volt) {
                    System.out.println(" Setting voltage : " + volt + " terminal with voltage: " + voltages[node]);
                    throw new Error("Why did we get a gradient zero?");
                }
            }
            if (!isTerminal[node]) {
                setVoltage(node, volt);
            }
            for (int i = g.nbrsStart[node]; i < g.nbrsStart[node] + g.deg[node]; i++) {
                int nbr = g.ai[i];
                if (!isTerminal[nbr])
                    nodeStack.push(nbr);
            }
        }
    }

    // find random terminal that has non-inf edges
    // only useful if can avoid a max in voltages,
    // and if the graph is connected
    // not used right now
    public int randActiveTerminal() {

        int numActive = 0;

        boolean[] isActive = new boolean[numTerminals];
        for (int i = 0; i < numTerminals; i++)
            isActive[i] = false;

        for (int i = 0; i < numTerminals; i++) {
            int node = terminalList[i];

//            if (g.deg[node] > 0) {
//                numActive++;
//                isActive[i] = true;
//
//            }
            for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                if (!isTerminal[g.ai[j]]) {
                    isActive[i] = true;
                }
            }

        }

        for (int i = 0; i < numTerminals; i++) {
            if (isActive[i])
                numActive++;
        }

        int node = -1;

        if (numActive == 0)
            throw new Error("No active terminals!");

        int ind = random.nextInt(numActive);
        int count = 0;

        for (int i = 0; i < numTerminals; i++) {
            if (isActive[i]) {
                if (count == ind) {
                    node = terminalList[i];
                    break;
                }
                count++;
            }
        }

        return node;

    }

    // compute the max gradient up or down from node, which should be a terminal
    public double gradientFromTerminal(int startNode) {

        // This procedure computes the gradient from a terminal.
        // It DOES take into account terminals that are directly connected to this terminal.
        // However, it does not allow the paths to run through terminals.
        // It runs Djikstra from node without passing through terminals

        double[] dist = new double[g.nv];
        boolean[] finished = new boolean[g.nv];

        for (int i = 0; i < g.nv; i++) {
            dist[i] = Double.POSITIVE_INFINITY;
            finished[i] = false;
        }

        dist[startNode] = 0;

        NodeHeap pq = new NodeHeap(g.nv);
        for (int i = g.nbrsStart[startNode]; i < g.nbrsStart[startNode] + g.deg[startNode]; i++) {
            pq.add(g.ai[i], g.av[i]);
        }
        finished[startNode] = true;

        while (pq.hasMore()) {
            int node = pq.popMin();
            finished[node] = true;
            double val = pq.keys[node];
            dist[node] = val;
            //logger.write(node + " : " + val);

            if (!isTerminal[node]) {
                for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                    int nbr = g.ai[j];
                    if (!finished[nbr]) {
                        double newVal = val + g.av[j];
                        pq.add(nbr, newVal);
                    }
                }
            }
        }


        double grad = 0;
        double nodeVolt = voltages[startNode];
        for (int i = 0; i < numTerminals; i++) {
            if (dist[terminalList[i]] > 0 && dist[terminalList[i]] != Double.POSITIVE_INFINITY) {
                double gtest = (voltages[terminalList[i]] - nodeVolt) / dist[terminalList[i]];
                if (gtest < 0) gtest = -gtest;
                if (grad < gtest) grad = gtest;
            }
        }

        return grad;
    }



    // return a boolean flag for whether or not edge exceeds pressure pr
    // it would be better to do this asym, using a graph for vlow and graph-transpose for vhigh
    // somehow, this managed to get some edges with a component that has just one terminal!?
    public boolean[] edgesStrictlyAbovePressure(double grad) {
        boolean[] flag = new boolean[g.ne];

        compLowAndHigh(grad);

        // we need a symmetric graph!
        // make sure not to grab terminal-terminal edges!!!!!

        for (int i = 0; i < g.ne; i++) {
            flag[i] = false;
            if (!(isTerminal[g.ai[i]] && isTerminal[g.aj[i]])) {
                flag[i] = ((vhigh[g.aj[i]] - vlow[g.ai[i]]) > grad * g.av[i] + 1e-10) ||
                        ((vhigh[g.ai[i]] - vlow[g.aj[i]]) > grad * g.av[i] + 1e-10);
            }
        }

        return flag;
    }

    public boolean[] edgesAboveOrEqualPressure(double grad) {
        boolean[] flag = new boolean[g.ne];

        compLowAndHigh(grad);

        // we need a symmetric graph!
        // make sure not to grab terminal-terminal edges!!!!!

        for (int i = 0; i < g.ne; i++) {
            flag[i] = false;
            if (!(isTerminal[g.ai[i]] && isTerminal[g.aj[i]])) {
                flag[i] = ((vhigh[g.aj[i]] - vlow[g.ai[i]]) >= (1 - 1e-8) * grad * g.av[i]) ||
                        ((vhigh[g.ai[i]] - vlow[g.aj[i]]) >= (1 - 1e-8) * grad * g.av[i]);
            }
        }

        return flag;
    }

    public void compLowAndHigh(double grad) {
//        double mult = ((double) 1) / pr;
        int[] where = new int[numTerminals];
        double[] what = new double[numTerminals];
        for (int i = 0; i < numTerminals; i++) {
            where[i] = terminalList[i];
            what[i] = voltages[terminalList[i]];
        }

        // SUSHANT: Changed the following from Djikstra, so that paths are not allowed to go through terminals
        vlow = noTerminalDijkstra(where, what, grad);

        for (int i = 0; i < numTerminals; i++) {
            what[i] = -what[i];
        }

        // SUSHANT: Changed the following from Djikstra, so that paths are not allowed to go through terminals
        vhigh = noTerminalDijkstra(where, what, grad);
        for (int i = 0; i < g.nv; i++) {
            vhigh[i] = -vhigh[i];
        }
    }

    // SUSHANT: Wrote this new version of Djikstra that doesn't allow paths to go through terminal vertices
    // Used to compute vhigh and vlow
    // The values returned at terminals is exactly the original voltage
    // Takes a POSITIVE parameter grad, which is the multiplicative factor for the distances
    public double[] noTerminalDijkstra(int[] where, double[] what, double grad) {
        if (grad < 0)
            throw new Error("Negative gradient!");

        double[] minVal = new double[g.nv];
        boolean[] finished = new boolean[g.nv];

        for (int i = 0; i < g.nv; i++) {
            minVal[i] = Double.POSITIVE_INFINITY;
            finished[i] = false;
        }

        //Logger logger = new Logger();
        //logger.start("dijkstra.log");

        NodeHeap pq = new NodeHeap(g.nv);
        for (int i = 0; i < where.length; i++) {
            if (!isTerminal[where[i]])
                throw new Error("Starting Djikstra from non-terminal!");
            pq.add(where[i], what[i]);
            //logger.write("(" + where[i] + ", " + what[i] + "), ");
        }

        while (pq.hasMore()) {
            int node = pq.popMin();
            finished[node] = true;
            double val = pq.keys[node];
            minVal[node] = val;
            //logger.write(node + " : " + val);

            for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                int nbr = g.ai[j];
                if (!finished[nbr] && !isTerminal[nbr]) {
                    double newVal = val + grad * g.av[j];
                    pq.add(nbr, newVal);
                }
            }
        }

        return minVal;
    }


    public void setVoltagesWithGrad(double grad) {
        double tol = 1e-10;
        grad = (1 - tol) * grad;

        boolean somethingSet = false;

        // First compute VHigh and Vlow while storing the voltages and the distances to the terminals assigning vhigh and vlow
        int[] where = new int[numTerminals];
        double[] what = new double[numTerminals];
        for (int i = 0; i < numTerminals; i++) {
            where[i] = terminalList[i];
            what[i] = voltages[terminalList[i]];
        }

        double[] lowTermVolt = new double[g.nv];
        double[] lowTermDist = new double[g.nv];
        boolean[] finished = new boolean[g.nv];

        for (int i = 0; i < g.nv; i++) {
            finished[i] = false;
        }

        NodeHeap pq = new NodeHeap(g.nv);
        for (int i = 0; i < where.length; i++) {
            if (!isTerminal[where[i]])
                throw new Error("Starting Djikstra from non-terminal!");
            int node = where[i];
            for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                int nbr = g.ai[j];
                boolean treeModified = pq.add(nbr, what[i] + grad * g.av[j]);
                if (treeModified) {
                    lowTermDist[nbr] = g.av[j];
                    lowTermVolt[nbr] = what[i];
                }
            }
//            pq.add(where[i], what[i]);
//            lowTermVolt[where[i]] = what[i];
//            lowTermDist[where[i]] = 0;
            //logger.write("(" + where[i] + ", " + what[i] + "), ");
        }

        while (pq.hasMore()) {
            int node = pq.popMin();
            finished[node] = true;
            double val = pq.keys[node];
            vlow[node] = val;
            assert ((val - lowTermVolt[node] - grad * lowTermDist[node]) < tol);
            assert ((val - lowTermVolt[node] - grad * lowTermDist[node]) > -tol);

            for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                int nbr = g.ai[j];
                if (!finished[nbr] && !isTerminal[nbr]) {
                    double newVal = val + grad * g.av[j];
                    boolean treeModified = pq.add(nbr, newVal);
                    if (treeModified) {
                        lowTermDist[nbr] = lowTermDist[node] + g.av[j];
                        lowTermVolt[nbr] = lowTermVolt[node];
                    }
                }
            }
        }

        double[] highTermVolt = new double[g.nv];
        double[] highTermDist = new double[g.nv];

        for (int i = 0; i < g.nv; i++) {
            finished[i] = false;
        }

        pq = new NodeHeap(g.nv);
        for (int i = 0; i < where.length; i++) {
            if (!isTerminal[where[i]])
                throw new Error("Starting Djikstra from non-terminal!");
            int node = where[i];
            for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                int nbr = g.ai[j];
                boolean treeModified = pq.add(nbr, -what[i] + grad * g.av[j]);
                if (treeModified) {
                    highTermDist[nbr] = g.av[j];
                    highTermVolt[nbr] = what[i];
                }
            }
//            pq.add(where[i], -what[i]);
//            highTermVolt[where[i]] = what[i];
//            highTermDist[where[i]] = 0;
            //logger.write("(" + where[i] + ", " + what[i] + "), ");
        }

        while (pq.hasMore()) {
            int node = pq.popMin();
            finished[node] = true;
            double val = pq.keys[node];
            vhigh[node] = -val;

            assert ((vhigh[node] - highTermVolt[node] + grad * highTermDist[node]) < tol);
            assert ((vhigh[node] - highTermVolt[node] + grad * highTermDist[node]) > -tol);

            if (vhigh[node] >= vlow[node]) {
                somethingSet = true;
                if (!isTerminal[node])
                    setVoltage(node, (lowTermVolt[node] * highTermDist[node] + highTermVolt[node] * lowTermDist[node]) / (lowTermDist[node] + highTermDist[node]));
            }

            for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                int nbr = g.ai[j];
                if (!finished[nbr] && !isTerminal[nbr]) {
                    double newVal = val + grad * g.av[j];
                    boolean treeModified = pq.add(nbr, newVal);
                    if (treeModified) {
                        highTermDist[nbr] = highTermDist[node] + g.av[j];
                        highTermVolt[nbr] = highTermVolt[node];
                    }
                }
            }
        }


        ////

//        // find every edge that is tight: if low end is not a terminal, make it one
//        for (int i = 0; i < g.ne; i++) {
//            if ((vhigh[g.aj[i]] - vlow[g.ai[i]]) > (1 - 1e-8) * grad * g.av[i]) {
//                if (!isTerminal[g.ai[i]]) {
//                    if ((vhigh[g.aj[i]] - vlow[g.ai[i]]) > (1 + 1e-8) * grad * g.av[i]) {
//                        System.out.println("Voltages: " + vhigh[g.aj[i]] + " " + vlow[g.ai[i]]);
//                        throw new Error("Fixing edge with the wrong gradient! ratio = " + ((vhigh[g.aj[i]] - vlow[g.ai[i]]) / (grad * g.av[i])));
//                    }
//                    int node = g.ai[i];
//                    if ((vhigh[node] - vlow[node] > 1e-8) || (vlow[node] - vhigh[node] > 1e-8))
//                        throw new Error("Too much diff betwen vhigh and vlow");
//
//                    setVoltage(node, (vlow[node] + vhigh[node]) / 2);
//
//                    // Testing if all neighbors have lower gradient
//                    for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
//                        int nbr = g.ai[j];
//                        if (isTerminal[nbr]) {
//                            if (voltages[nbr] - voltages[node] > (1 + 1e-8) * grad * g.av[j])
//                                throw new Error("Too much gradient!");
//                            if (voltages[node] - voltages[nbr] > (1 + 1e-8) * grad * g.av[j])
//                                throw new Error("Too much gradient!");
//                        }
//                    }
//
//                    somethingSet = true;
//                }
//            }
//        }

        if (!somethingSet) {
            throw new Error("Nothing Set!");
            // System.out.println("Setting all nodes in component");
            // setNodesFromTerminalWithZeroGradient(terminalList[0]);
        }
    }


    // use highv and lowv to set terminals, as appropriate
    // needs the gradient used to compute lowv and highv
    // Warning: Needs grad to be the gradient of the tightest path
    public void setFromLowHigh(double grad) {
        boolean somethingSet = false;

        // find every edge that is tight: if low end is not a terminal, make it one
        for (int i = 0; i < g.ne; i++) {
            if ((vhigh[g.aj[i]] - vlow[g.ai[i]]) > (1 - 1e-8) * grad * g.av[i]) {
                if (!isTerminal[g.ai[i]]) {
                    if ((vhigh[g.aj[i]] - vlow[g.ai[i]]) > grad * g.av[i] + 1e-10) {
                        System.out.println("Voltages: " + vhigh[g.aj[i]] + " " + vlow[g.ai[i]]);
                        throw new Error("Fixing edge with the wrong gradient! ratio = " + ((vhigh[g.aj[i]] - vlow[g.ai[i]]) / (grad * g.av[i])));
                    }
                    int node = g.ai[i];
                    if ((vhigh[node] - vlow[node] > 1e-8) || (vlow[node] - vhigh[node] > 1e-8))
                        throw new Error("Too much diff betwen vhigh and vlow");

                    setVoltage(node, (vlow[node] + vhigh[node]) / 2);

                    // Testing if all neighbors have lower gradient
                    for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                        int nbr = g.ai[j];
                        if (isTerminal[nbr]) {
                            if (voltages[nbr] - voltages[node] > grad * g.av[j] + 1e-10)
                                throw new Error("Too much gradient!");
                            if (voltages[node] - voltages[nbr] > grad * g.av[j] + 1e-10)
                                throw new Error("Too much gradient!");
                        }
                    }

                    somethingSet = true;
                }
            }
        }

        if (!somethingSet) {
            throw new Error("Nothing Set!");
            // System.out.println("Setting all nodes in component");
            // setNodesFromTerminalWithZeroGradient(terminalList[0]);
        }
    }


    /* a subgraph isn't really just a subgraph:
      it also needs a mapping to vertices in the original.
     */
    public static class Subgraph {
        public RandLexVolt rlv;
        public int[] mapBack;
        public double shift;

        public Subgraph() {
        }
    }

    // return the edge-induced subgraph with edges indicated by flag
    // assumes that the graph indicated by flag is symmetric, otherwise could fail
    public Subgraph edgeInducedSubgraph(boolean[] flag) {

        // this maps nodes in the orig to nodes in the subgraph
        int[] nodeProj = new int[g.nv];
        for (int i = 0; i < g.nv; i++)
            nodeProj[i] = -1;

        int nvSub = 0;
        int neSub = 0;

        int prevNode = -1;

        int numSubTerminals = 0;

        // find the nodes that appear, and name them
        for (int i = 0; i < g.ne; i++) {
            if (flag[i]) {
                neSub++;
                int node = g.aj[i];
                if (node != prevNode) {
                    prevNode = node;
                    nodeProj[node] = nvSub++;
                    if (isTerminal[node])
                        numSubTerminals++;
                }
            }
        }

        int[] mapBack = new int[nvSub];

        int[] subTerminals = new int[numSubTerminals];
        double[] subVoltages = new double[numSubTerminals];

        int subTerminalInd = 0;

        for (int i = 0; i < g.nv; i++) {
            if (nodeProj[i] >= 0) {
                mapBack[nodeProj[i]] = i;
                if (isTerminal[i]) {
                    subVoltages[subTerminalInd] = voltages[i];
                    subTerminals[subTerminalInd++] = nodeProj[i];
                }
            }
        }

        // now shift voltages to improve accuracy
        double shift = getMinVoltage(subVoltages);
        for(int i = 0; i<subVoltages.length; i++) {
            subVoltages[i] = subVoltages[i]-shift;
        }

        int[] subi = new int[neSub];
        int[] subj = new int[neSub];
        double[] subv = new double[neSub];

        int ePtr = 0;

        for (int i = 0; i < g.ne; i++) {
            if (flag[i]) {
                subi[ePtr] = nodeProj[g.ai[i]];
                subj[ePtr] = nodeProj[g.aj[i]];
                subv[ePtr] = g.av[i];
                ePtr++;
            }
        }


        EdgeListGraph subg = new EdgeListGraph(subi, subj, subv);
        RandLexVolt rlv = new RandLexVolt(subg,random.nextLong());
        rlv.setTerminals(subTerminals, subVoltages);

        Subgraph sg = new Subgraph();
        sg.mapBack = mapBack;
        sg.rlv = rlv;
        sg.shift = shift;

        return sg;
    }

    public void setTerminalsFromSubgraph(Subgraph sub) {
        for (int i = sub.rlv.numOrigTerminals; i < sub.rlv.numTerminals; i++) {
            int node = sub.mapBack[sub.rlv.terminalList[i]];
            setVoltage(node, sub.rlv.voltages[sub.rlv.terminalList[i]]+sub.shift);

        }
    }

    public boolean isGraphConnected() {
        Stack<Integer> s = new Stack<Integer>();

        boolean[] seen = new boolean[g.nv];
        boolean[] completed = new boolean[g.nv];
        for (int i = 0; i < g.nv; i++) {
            seen[i] = false;
            completed[i] = false;
        }

        s.push(0);
        seen[0] = true;

        while (!s.isEmpty()) {
            int node = s.pop();
            if (!completed[node]) {
                for (int i = g.nbrsStart[node]; i < g.nbrsStart[node] + g.deg[node]; i++) {
                    int nbr = g.ai[i];
                    if (!seen[nbr]) {
                        s.push(nbr);
                        seen[nbr] = true;
                    }
                }
                completed[node] = true;
            }
        }

        for (int i = 0; i < g.nv; i++)
            if (!completed[i])
                return false;

        return true;
    }

    public void testFinalVoltages() {
        if (numTerminals < g.nv)
            throw new Error("Test only after all voltages have been computed");
        double maxError = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < g.nv; i++) {
            if (isOrigTerminal[i])
                continue;
            double maxgrad = 0;
            double mingrad = 0;
            for (int j = g.nbrsStart[i]; j < g.nbrsStart[i + 1]; j++) {
                int nbr = g.ai[j];
                double grad = (voltages[nbr] - voltages[i]) / g.av[j];
                if (grad > maxgrad)
                    maxgrad = grad;
                if (grad < mingrad)
                    mingrad = grad;
            }
//            System.out.println("Vertex " + i + " maxgrad: " + maxgrad + " mingrad: " + mingrad);
            if ((maxgrad + mingrad) > maxError || (maxgrad + mingrad) < -maxError)
                maxError = Math.abs(maxgrad + mingrad);
        }

 /*       if (maxError > testEps)
            System.out.println("Error in gradients! Max + min gradient = " + maxError );
        else
        System.out.println("Tested final voltages");*/

    }

    /* ** NEW SECTION : SAMPLING EDGES AND COMPUTING THEIR GRADIENTS */

    private int sampleActiveEdge() {
        //we don't think the subgraph should contain any term-term edges
        int randEdgeInd = random.nextInt(g.ne);
        return randEdgeInd;
    }

    public double computeEdgeGradient(int e) {
        int src = g.aj[e];
        int dest = g.ai[e];
        if (isTerminal[src] && isTerminal[dest])
            throw new Error("Why are we computing the gradient of a term-term edge");
        //we should not compute the gradient of a term-term edge

        // Perform Djikstra from the src of the edge
        findTermsFromLeftVertex(src);

        // Perform Djikstra from the dest of the edge
        findTermsFromRightVertex(dest);

//        double grad = edgeGradientHelper(edges[e].len, leftV, leftD, rightV, rightD);
        double grad;
//        try {
        grad = edgeGradientHelper.computeRandGrad(g.av[e]);
//        } catch(Error e){
//            System.out.println("Heap Usage count : " + heap.getUsageCount());
//            throw e;
//        }

        return grad;
    }

    // Performs Djikstra from the node v
    // Adds the terminals discovered to edgeGradientHelper as left terminals
    private void findTermsFromLeftVertex(int v) {
        if (isTerminal[v]) {
            edgeGradientHelper.addLeftTerm(voltages[v], 0.0);
        } else {
            // Perform Djikstra only if v is not a terminal
            assert(heap != null);
            heap.resetHeap();
            heap.add(v, 0);

            while (heap.hasMore()) {
                double val = heap.minKey();
                int node = heap.popMin();

                if (isTerminal[node]) {
                    edgeGradientHelper.addLeftTerm(voltages[node], val);
                } else {

                    // ****** IMP ******
                    // THIS SEEMS TO BE THE SAFEST WAY TO TRAVERSE THE ADJ LIST
                    // This traverses only active edges
                    // ****** IMP ******
//                    int eptr = vertices[node].actNbrsEnd;
//                    for (int i = 0; i < vertices[node].actDeg; i++) {
//                        assert (edges[eptr].src == node);
//                        int nbr = edges[eptr].dest;
//
//                        if (!heap.finishedThisTime(nbr)) { // Only if the vertex isn't already finished in this iteration of Djikstra
//                            heap.add(nbr, val + edges[eptr].len); // NOTE: Even if the node already exists, the heap handles it
//                        }
//
//                        eptr = edges[eptr].prevEdge;
//                    }

                    //TODO: I may well have refactored this wrong!
                    for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                        int nbr = g.ai[j];
                        if (!heap.finishedThisTime(nbr)) {
                            heap.add(nbr, val + g.av[j]);
                        }
                    }
                }
            }
        }
    }

    // Performs Djikstra from the node v
    // Adds the terminals to edgeGradientHelper as right terminals
    private void findTermsFromRightVertex(int v) {
        if (isTerminal[v]) {
            edgeGradientHelper.addRightTerm(voltages[v], 0.0);
            return;
        }
        // Perform Djikstra only if v is not a terminal

        heap.resetHeap();
        heap.add(v, 0);

        while (heap.hasMore()) {
            double val = heap.minKey();
            int node = heap.popMin();

            if (isTerminal[node]) {
                edgeGradientHelper.addRightTerm(voltages[node], val);
            } else {

//                // ****** IMP ******
//                // THIS SEEMS TO BE THE SAFEST WAY TO TRAVERSE THE ADJ LIST
//                // This traverses only active edges
//                // ****** IMP ******
//                int eptr = vertices[node].actNbrsEnd;
//                for (int i = 0; i < vertices[node].actDeg; i++) {
//                    assert (edges[eptr].src == node);
//                    int nbr = edges[eptr].dest;
//
//                    if (!heap.finishedThisTime(nbr)) { // Only if the vertex isn't already finished in this iteration of Djikstra
//                        heap.add(nbr, val + edges[eptr].len); // NOTE: Even if the node already exists, the heap handles it
//                    }
//
//                    eptr = edges[eptr].prevEdge;
//                }

                //TODO: I may well have refactored this wrong!
                for (int j = g.nbrsStart[node]; j < g.nbrsStart[node] + g.deg[node]; j++) {
                    int nbr = g.ai[j];
                    if (!heap.finishedThisTime(nbr)) {
                        heap.add(nbr, val + g.av[j]);
                    }
                }
            }
        }
    }

    // Code to test final voltages
    public void testFinalVoltagesScaled() {
            System.out.println("Testing final voltages with scaling");
        if (numTerminals < g.nv)
            throw new Error("Test only after all voltages have been computed");
        for (int i = 0; i < g.nv; i++) {
            if (isOrigTerminal[i])
                continue;
            double maxgrad = Double.NEGATIVE_INFINITY;
            double maxvolt = voltages[i];
            double maxdist = 0;
            double mingrad = Double.POSITIVE_INFINITY;
            double minvolt = voltages[i];
            double mindist = 0;


            for (int j = g.nbrsStart[i]; j < g.nbrsStart[i + 1]; j++) {
                int nbr = g.ai[j];
                double grad = (voltages[nbr] - voltages[i]) / g.av[j];
                if (grad > maxgrad) {
                    maxgrad = grad;
                    maxvolt = voltages[nbr];
                    maxdist = g.av[j];
                }
                if (grad < mingrad) {
                    mingrad = grad;
                    minvolt = voltages[nbr];
                    mindist = g.av[j];
                }
            }

    /*        double testVolt = (mindist * maxvolt + maxdist * minvolt) / (mindist + maxdist);
            double scale = Double.longBitsToDouble(Double.doubleToLongBits(Math.abs(testVolt) + Math.abs(voltages[i])) & 0x7ff0000000000000L);
            if ((testVolt - voltages[i]) > testEps * scale || (testVolt - voltages[i]) < -testEps * scale) {
                System.out.println("Warning: There seems to be some error in the final voltages v.voltage: " + voltages[i] + " testvoltage: " + testVolt);
            }*/

        }
        /*System.out.println("Tested final voltages w scaling");*/
    }

    // -- voltage shift for increased accuracy
    private static double getMinVoltage(double[] subVoltages) {
        double minV = Double.POSITIVE_INFINITY;
        for(int i = 0; i<subVoltages.length; i++) {
            if (minV > subVoltages[i]) minV = subVoltages[i];
        }
        return minV;
    }

}
