// Kruskal's algorithm in C++ with sorted edges

#include <iostream>
#include <vector>
using namespace std;

#define edge pair<int, int>

long num_comparisons = 0;
long num_edge = 0;


void increment_comparison()
{
	num_comparisons++;
	
}

class Graph {
   private:
  vector<pair<int, edge> > G;  // graph
  vector<pair<int, edge> > T;  // mst
  int *parent;
  int V;  // number of vertices in graph
  public:
  Graph(int V) {
  parent = new int[V];

  for (int i = 0; i < V; i++)
    parent[i] = i;

}

// function to add edges
  void AddWeightedEdge(int u, int v, int w){
  G.push_back(make_pair(w, edge(u, v)));
}

// find set function to remove cycles from mst
  int find_set(int i){
  
  increment_comparison();

  // If i is the parent of itself
  if (i == parent[i])
    return i;
  else
    // Else if i is not the parent of itself
    // Then then some other vertex is parent of this set,
    // so recursively calling find_set on its parent
    return find_set(parent[i]);
}
// union_set for uniting disjoint sets to find mst
  void union_set(int u, int v){
  parent[u] = parent[v];
}

// kruskal algorithm function
  void kruskal(){
  int i, source, target;
  sort(G.begin(), G.end());  // sorting in ascending order

  // traversing through the graph
  for (i = 0; i < G.size(); i++) {
    // G[i] has two values as it is a pair.
    //Here we are trying to find source and target of the current edge for mst
    source = find_set(G[i].second.first);
    target = find_set(G[i].second.second);
    if (source != target) {
      T.push_back(G[i]);  // add to tree
      num_edge ++ ;
      union_set(source, target);
    }
  }
}
  void print(){
  long total_weight = 0 ;



  for (int i = 0; i < T.size(); i++) {
    // cout << "e "<<T[i].second.first<<" "<<T[i].second.second << " "<< T[i].first<< endl;

    total_weight += T[i].first;
    
  }

  cout<<"total weight\t"<<total_weight<<endl;
}
};


int main() {
    char c;
	long n = 0; // number of vertices
	long m = 0; // number of edges
	int flag = 0;
	int edge_counter = 0; // for counting number of edges
    
    Graph g(n);

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
                g.AddWeightedEdge(source, target, weight);
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
    	g.kruskal();

    auto end_time = chrono::high_resolution_clock::now();

    cerr<< "g "<<n<<" "<<num_edge<<endl;

    cerr << "runtime\t" << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()<< endl;
    cerr << "number of comparisons\t" << num_comparisons<<endl ;
    // cerr<< "number of edges "<<num_edge<<endl;
  g.print();

  if(num_edge!= n-1) {
		 cerr<<"Error! Disconnected graph provided"<<endl;
	}
  return 0;
}