import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Stack;

/**
 * 
 * @author Lexi Reicks and Jacob Cram
 *
 */
public class WGraph {
	/**
	 * Number of vertices in the graph
	 */
	public int V;
	/**
	 * Number of edges in the graph
	 */
	public int E;
	/**
	 * Height of the graph
	 */
	public int H;
	/**
	 * Width of the graph
	 */
	public int W;
	/**
	 * Adjacency list holding the edges
	 */
	public LinkedList<Edge> adj[];
	/**
	 * List of vertices
	 */
	public Vertex[] vertices;
	
	class Vertex {
		private final int index;
		private final int row;
		private final int col;
		
		public Vertex(int index, int row, int col) {
			this.index = index;
			this.row = row;
			this.col = col;
		}
		
		public int index() {
			return index;
		}
		
		public int row() {
			return row;
		}
		
		public int col() {
			return col;
		}
		
		public int[] coords() {
			return new int[] { row, col };
		}
	}
	
	public class Edge { 
	    private final int src;
	    private final int dest;
	    private int weight;

	    public Edge(int src, int dest, int weight) {
	        this.src = src;
	        this.dest = dest;
	        this.weight = weight;
	    }

	    public int weight() {
	        return weight;
	    }

	    public int src() {
	        return this.src;
	    }
	    
	    public int dest() {
	    	return this.dest;
	    }
	}
	
	public class IndexedMinPQ {
		private int[] heap;
		private Integer[] keys;
		private int[] index;
		private int size;
		private int max;
		
		public IndexedMinPQ(int max) {
			this.max = max;
			this.heap = new int[max + 1];
			this.keys = new Integer[max + 1];
			this.index = new int[max + 1];
			this.size = 0;
			
			for (int i = 0; i < max; i++) {
				this.index[i] = -1;
			}
		}

		/**
		 * Adds a string s with priority p to the priority queue
		 * @param s
		 * @param p
		 */
		public void add(int i, int p) {
	        if (i < 0 || i >= max) throw new IllegalArgumentException();
	        if (contains(i)) throw new IllegalArgumentException("index is already in the priority queue");
			size++;
			index[i] = size;
			// add a new priority node to the end
			heap[size] = i;
			keys[i] = p;
			// heapify from the bottom up
			heapifyUp(size);
		}

		/**
		 * returns a string whose priority is maximum
		 * @return
		 */
		public int returnMin() {
			return this.heap[1];
		}

		/**
		 * returns a string whose priority is minimum and removes if form the queue
		 * @return
		 */
		public int extractMin() {
	        if (size == 0) throw new NoSuchElementException("Priority queue underflow");
	        int min = heap[1];
	        swap(1, size--);
	        heapifyDown(1);
	        index[min] = -1;
	        keys[min] = null;
	        heap[size + 1] = -1;
	        return min;
		}

		/**
		 * Removes the element from the priority queue whose array index is i
		 * @param i
		 * @return
		 */
		public int remove(int i) {
	        int index = this.index[i];
	        swap(index, size--);
	        heapifyUp(index);
	        heapifyDown(index);
	        keys[i] = null;
	        this.index[i] = -1;
			
			return index;
		}

		/**
		 * Decrements the priority of the ith element by k
		 * @param i
		 * @param k
		 */
		public void decrementPriority(int i, int k) {
	        if (i < 0 || i >= max) throw new IllegalArgumentException();
	        if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
	        if (keys[i].compareTo(k) <= 0)
	            throw new IllegalArgumentException("Calling decreaseKey() with given argument would not strictly decrease the key");
	        keys[i] = k;
	        heapifyUp(index[i]);
		}
		
		public boolean contains(int i) {
	        if (i < 0 || i >= max) throw new IllegalArgumentException();
	        return index[i] != -1;
		}
		
		/**
		 * Returns the key(A[i])
		 * where A is the array used to represent the priority queue
		 * @param i - index of array
		 * @return key
		 */
		public int getKey(int i) {
			return this.keys[i];
		}
		
		/***
		 * Returns true if and only if the queue is empty
		 * @return
		 */
		public boolean isEmpty() {
			return size == 0;
		}
		
		/**
		 * Returns the size of the heap
		 * @return size of the array
		 */
		public int size() {
			return size;
		}
		
		private boolean greater(int i, int j) {
			return keys[heap[i]].compareTo(keys[heap[j]]) > 0;
		}
		
	    private void heapifyUp(int k) {
	        while (k > 1 && greater(k / 2, k)) {
	            swap(k, k / 2);
	            k = k / 2;
	        }
	    }

	    private void heapifyDown(int k) {
	        while (2 * k <= size) {
	            int j = 2 * k;
	            if (j < size && greater(j, j + 1)) j++;
	            if (!greater(k, j)) break;
	            swap(k, j);
	            k = j;
	        }
	    }
		
		private void swap(int x, int y) {
	        int swap = heap[x];
	        heap[x] = heap[y];
	        heap[y] = swap;
	        index[heap[x]] = x;
	        index[heap[y]] = y;
		}
	}
	
	class Dijkstra {
		private int[] distTo;
	    private Edge[] edgeTo;
	    private IndexedMinPQ pq;
	    Dijkstra(int ux, int uy) {
	    	int src = findVertex(ux, uy);
			
			if (src < 0) {
				throw new IllegalArgumentException("vertex (" + ux + ", " + uy + ") does not exist in the graph");
			}
			
			distTo = new int[V];
			edgeTo = new Edge[V];
			
			for (int i = 0; i < V(); i++) {
				distTo[i] = Integer.MAX_VALUE;
			}
			
			distTo[src] = 0;
			
			pq = new IndexedMinPQ(V);
			pq.add(src, distTo[src]);
			
			while (!pq.isEmpty()) {
				int v = pq.extractMin();
				
				for (Edge e: adj(v)) {
					updatePQ(e);
				}
			}
		}
	    
	    private void updatePQ(Edge e) {
	        int v = e.src(), w = e.dest();
	        if (distTo[w] > distTo[v] + e.weight()) {
	            distTo[w] = distTo[v] + e.weight();
	            edgeTo[w] = e;
	            if (pq.contains(w)) {
	            	pq.decrementPriority(w, distTo[w]);
	            } else {
	            	pq.add(w, distTo[w]);
	            }
	        }
	    }
		
		public Stack<Edge> pathTo(int v) {
			if (!(distTo[v] < Integer.MAX_VALUE)) return null;
			Stack<Edge> path = new Stack<Edge>();
			for (Edge e = edgeTo[v]; e != null; e = edgeTo[e.src()]) {
				path.add(e);
			}
			return path;
		}
	}
	
	/**
	 * File:
	 * 
	 * First line - contains a number indicating the number of vertices in the graph
	 * Second line - contains a number indicating the number of edges in the graph
	 * 
	 * All other lines- 5 numbers 
	 * source vertex coordinates (first two numbers),
	 * destination vertex coordinates (third and fourth numbers)
	 * and weight of the edge connecting the source vertex to the destination (assume direction of edge from source to destination)
	 */
	public WGraph(String FName) throws NumberFormatException, IOException {
		File f = new File(FName);
		
		this.parseFile(f);
	}
	
	public void parseFile(File f) throws FileNotFoundException {
		Scanner br = new Scanner(f);
		
		this.V = br.nextInt();
		this.E = br.nextInt();
		
		this.adj = new LinkedList[V];
		this.vertices = new Vertex[this.V];
		
		for (int i = 0; i < this.V; i++) {
			this.adj[i] = new LinkedList<Edge>();
			this.vertices[i] = null;
		}
		
		while (br.hasNextLine()) {
			int srcX = br.nextInt();
			int srcY = br.nextInt();
			
			int srcIndex = findVertex(srcX, srcY);
			
			int destX = br.nextInt();
			int destY = br.nextInt();
			
			int destIndex = findVertex(destX, destY);
		
			int edgeWeight = br.nextInt();
			
			addEdge(srcIndex, destIndex, edgeWeight);
		}

		br.close();
	}
	
	/**
	 * Given vertices u and v, find the shortest path from u to v. We will refer to 
	 * such path as V2V(u, v)
	 * 
	 * 	V2V(u, v) = argmin cost(π | src(π) = u and dest(π) = v)
	 * 
	 * @param ux - vertex u - x coord
	 * @param uy - vertex u - y coord
	 * @param vx - vertex v - x coord
	 * @param vy - vertex v - y coord
	 * @return arraylist contains even number of integers,
	 * 		   for any even i,
	 * 		   i-th and i+1-th integers in the array represent
	 * 		   the x-coordinate and y-coordinate of the i/2-th vertex
	 * 		   in the returned path (path is an ordered sequence of vertices)
	 */
	public ArrayList<Integer> V2V(int ux, int uy, int vx, int vy) {
		ArrayList<Integer> output = new ArrayList<Integer>();
		
		if (!validateVertex(ux, uy) || !validateVertex(vx, vy)) {
			return output;
		}
		
		output.add(ux);
		output.add(uy);
		Dijkstra d = new Dijkstra(ux, uy);
		
		int v = findVertex(vx, vy);
		Stack<Edge> edges = d.pathTo(v);
		while (!edges.isEmpty()) {
			Edge e = edges.pop();
			int[] dest = findCoords(e.dest());
			output.add(dest[0]);
			output.add(dest[1]);
		}
		return output;
	}
	
	/**
	 * Given a vertex u and a set of vertices S, find the shortest path from u to 
	 * some vertex in S. We will refer to such a path as V2S(u, S)
	 * 
	 * 	V2S(u, S) = argmin cost(π | src(π) = u and dest(π) ∈ S
	 * 
	 * @param ux - x coord of vertex u
	 * @param uy - y coord of vertex y
	 * @param 	S - set a vertices
	 * 	      	The arraylist S contains an even number of integers
	 * 		  	for any even i,
	 * 			i-th and i+1-th integers in the array represent
	 * 			the x-coord and y-coord of i/2-th vertex
	 * 			in the set S
	 * @return arraylist contains even number of integers,
	 * 		   for any even i,
	 * 		   i-th and i+1-th integers in the array represent
	 * 		   the x-coordinate and y-coordinate of the i/2-th vertex
	 * 		   in the returned path (path is an ordered sequence of vertices)
	 */
	public ArrayList<Integer> V2S(int ux, int uy, ArrayList<Integer> S) {
		ArrayList<Integer> output = new ArrayList<Integer>();
		
		if (!validateVertex(ux, uy)) {
			return output;
		}
		
		output.add(ux);
		output.add(uy);
		int[] dist = new int[S.size() / 2];
		int min = Integer.MAX_VALUE;
		int minIndex = -1;
		
		for (int i = 0; i < S.size() / 2; i++) {
			 dist[i] = Integer.MAX_VALUE;
		}
		
		Dijkstra d = new Dijkstra(ux, uy);
		
		for (int i = 0; i < S.size() / 2; i++) {
			int v = findVertex(S.get(i * 2), S.get(i * 2 + 1));
			dist[i] = d.distTo[v];
			
			if (dist[i] < min) {
				min = dist[i];
				minIndex = v;
			}
		}
		
		if (minIndex == -1) {
			return output;
		}
		
		Stack<Edge> path = d.pathTo(minIndex);
		
		while (!path.isEmpty()) {
			Edge e = path.pop();
			int[] dest = findCoords(e.dest());
			output.add(dest[0]);
			output.add(dest[1]);
		}
		
		return output;
	}
	
	/**
	 * 
	 * Given a set of vertices S1 and a set of vertices S2, find the shortest path
	 * from some vertex in S1 to some vertex in S2. We will refer to such a path as S2S(S1, S2)
	 * 
	 * 	S2S(S1, S2) = argmin cost(π | src(π) ∈ S1 and dest(π) ∈ S2
	 * 
	 * @param S1 - represents set of vertices
	 * @param S2 - represents set of vertices
	 * @return
	 */
	public ArrayList<Integer> S2S(ArrayList<Integer> S1, ArrayList<Integer> S2) {
		ArrayList<Integer> output = new ArrayList<Integer>();
		Dijkstra d = null;
		int min = Integer.MAX_VALUE;
		
		Dijkstra minD = null;
		int S2Vertex = -1;
		int S1Vertex = -1;
		
		for (int i = 0; i < S1.size() / 2; i++) {
			int src = findVertex(S1.get(i * 2), S1.get(i * 2 + 1));
			d = new Dijkstra(row(src), col(src));
			
			for (int j = 0; j < S2.size() / 2; j++) {
				int dest = findVertex(S2.get(j * 2), S2.get(j * 2 + 1));
				int dist = d.distTo[dest];
				
				if (dist < min) {
					min = dist;
					minD = d;
					S1Vertex = src;
					S2Vertex = dest;
				}
			}
		}
		
		if (minD == null) {
			return output;
		}
		
		Stack<Edge> path = minD.pathTo(S2Vertex);
		
		output.add(this.vertices[S1Vertex].row);
		output.add(this.vertices[S1Vertex].col);
		
		while (!path.isEmpty()) {
			Edge e = path.pop();
			Vertex dest = this.vertices[e.dest()];
			output.add(dest.row);
			output.add(dest.col);
		}
		
		return output;
	}
	
	/**
	 * Getter methods
	 */
	
	public LinkedList<Edge> adj(int v) {
		return adj[v];
	}
	
	public int E() {
		return E;
	}
	
	public int V() {
		return V;
	}
	
	/**
	 * Helper Methods
	 */
	public void addEdge(int src, int dest, int weight) {
		Edge toAdd = new Edge(src, dest, weight);
		this.adj[src].add(toAdd);
	}
	
	public boolean validateVertex(int ux, int uy) {
		for (Vertex v: this.vertices) {
			if (v.row() == ux && v.col() == uy) {
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Gets the coordinates of a vertex
	 * @param index
	 * @return
	 */
	public int[] findCoords(int index) {
		Vertex v = this.vertices[index];
		return v.coords();
	}
	
	/**
	 * Gets the row of a vertex
	 * @param index
	 * @return
	 */
	private int row(int index) {
		return this.vertices[index].row;
	}
	
	/**
	 * Gets the column of a vertex
	 * @param index
	 * @return
	 */
	private int col(int index) {
		return this.vertices[index].col;
	}
	
	/**
	 * Finds a vertex's index by row/col
	 * @param row
	 * @param col
	 * @return
	 */
	public int findVertex(int row, int col) {
		int vertex = -1;

		for (int i = 0; i < this.V; i++) {
			if (this.vertices[i] == null) {
				this.vertices[i] = new Vertex(i, row, col);
				vertex = i;
				return vertex;
			}
			if (this.vertices[i].row == row && this.vertices[i].col == col) {
				vertex = i;
				return vertex;
			}
		}
		
		return vertex;
	}
}
