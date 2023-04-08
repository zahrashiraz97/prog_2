
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <string>


using namespace std;

long num_comparisons = 0;

void increment_comparison()
{
	num_comparisons++;
}

// kruskal_v data structure

// path compression + rank by union
class kruskal_v {
	int* parent;
	int* rank;

public:
	kruskal_v(int n)
	{
		parent = new int[n];
		rank = new int[n];

		for (int i = 0; i < n; i++) {
			parent[i] = -1;
			rank[i] = 1;
		}
	}

	
	int find(int i)    // Find function
	{
        // increment_comparison();
		if (parent[i] == -1)
			return i;

		return parent[i] = find(parent[i]);
	}

	
	void unite(int u, int v)   // Union function
	{
		int s = find(u);
		int t = find(v);

        increment_comparison();

		if (s != t) {
            increment_comparison();
			if (rank[s] < rank[t]) {
				parent[s] = t;
			}
			else if (rank[s] > rank[t]) {
                increment_comparison() ;
				parent[t] = s;
			}
			else {
				parent[t] = s;
				rank[s] += 1;
			}
		}
	}
};

class Graph {
	vector<vector<int>> edgelist;
	int V;

public:
	Graph() {}
    Graph(int V) { this->V = V; }

	// Function to add edge in a graph
	vector<vector<int>> addEdge(int x, int y, int w)
	{
		edgelist.push_back({ w, x, y });
        return edgelist;
	}

	void secondkruskals_mst()
	{
		// Sort all edges
		sort(edgelist.begin(), edgelist.end());

		// Initialize the kruskal_v
		kruskal_v s(V);
		int ans = 0;
		cerr << "Following are the edges in the "
				"constructed MST"
			<< endl;
		for (auto edge : edgelist) {
			int w = edge[0];
			int x = edge[1];
			int y = edge[2];

			// Take this edge in MST if it does
			// not forms a cycle
			if (s.find(x) != s.find(y)) {
				s.unite(x, y);
				ans += w;
				cout << x << " -- " << y << " == " << w
					<< endl;
			}
		}
		cout << "Minimum Cost Spanning Tree: " << ans<<endl;
	}
};

// // Driver code
// int main()
// {
// 	Graph g(5);
// 	g.addEdge(1, 2, 1);
// 	g.addEdge(2, 3, 2);
// 	g.addEdge(3, 4, 1);
// 	g.addEdge(4, 5, 1);
// 	g.addEdge(1, 5, 1);

// 	// Function call
// 	g.kruskals_mst();

// 	return 0;
// }

int main() {

    char c;
	long n = 0;
	long m = 0;
	int flag = 0;
	int edge_counter = 0;
     Graph g;
     int i = 0;
	while((c = cin.peek())) {
		if (c == 'c')
		{
	    	cin.ignore(numeric_limits<streamsize>::max(), '\n');
		}
		else if (c == 'n')
		{
	    	cin.ignore(numeric_limits<streamsize>::max(), '\n');
		}
	    else if ((c == 'g') && (flag != 1))
	    {
	    	cin>>c;
	    	cin>>n;
	    	cin>>m;
	    	flag = 1;
	    	edge_counter = m;
            g = Graph(n);
	    }
	    else if ((c == 'e') && (flag == 1))
	    {
	    	if (edge_counter > 0)
	    	{
	    		long source;
	    		long target;
	    		long weight;
	    		cin>>c;
	    		cin>>source;
	    		cin>>target;
	    		cin>>weight;
                g.addEdge(source, target, weight);
                i++;
	    		edge_counter--;
	    		if (edge_counter == 0)
	    		{
	    			break;
	    		}
	    	}
	    	else
	    	{
	    		cin.ignore(numeric_limits<streamsize>::max(), '\n');

	    	}

	    }
	    else
	    {
	    	cin.ignore(numeric_limits<streamsize>::max(), '\n');
	    }

	}


    auto start_time = chrono::high_resolution_clock::now();
    	g.secondkruskals_mst();

    auto end_time = chrono::high_resolution_clock::now();

    cout<<endl <<"number of comparisons: "<<num_comparisons<<endl;


    return 0;

}



