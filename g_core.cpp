#include <algorithm>

#include "g.h"

using namespace std;

g::g( int vsize ) :
  n( vsize )
{
  set_up();
} 

// copy constructor
g::g( const g &otherG ){
  n = otherG.order();
  set_up();
  for( int i = 0; i < n; i++ ){
    for( int j = 0; j < n; j++ ){
      if( otherG.is_edge(i,j) ){
	add_edge(i,j);
      }
    }
  }
}

// Destructor
g::~g(){
  delete[] vdegree;
}
 
void g::set_up(){
  oldN = n;
  arraySize = n / intSize;
  if( ( n % intSize ) != 0 ) arraySize++;
  gA.resize( n );
  gB.resize( n );
  /*  for( int i = 0; i < n; i++ ){
    gA[i].resize( arraySize, 0 );
    gB[i].resize( arraySize, 0 );
    }*/
  numEdges = 0;
  vdegree = new int[n];

  for( int i = 0; i < n; i++ ){
    vdegree[i] = 0;

    for( int b = i + 1; b < n; b++ ){
      set_insert( b, gB[i] );
    }
  }

  trisCalced = false;
}


int g::order() const{
  return n;
}


int g::num_edges(){
  return numEdges;
}

vset g::get_adj( int v ){
  return gA[v];
}


void g::add_edge( int u, int v ){
  if( u < n && u >= 0 && v < n && v >= 0 ){
    if( !is_edge(u,v) ){
      numEdges++;
      set_insert( u, gA[v] );
      set_insert( v, gA[u] );
      vdegree[u]++;
      vdegree[v]++;
    }
  }
  else{
    cerr << "Error: Invalid vertices " << u << " and " << v << endl;
    cerr << "Size of graph = " << n << endl;
  }
}

void g::remove_edge( int u, int v ){
  if( u < n && u >= 0 && v < n && v >= 0
      && in_set(u, gA[v]) && in_set(v, gA[u]) ){
    if( is_edge(u,v) ){
      numEdges--;
      set_delete( u, gA[v] );
      set_delete( v, gA[u] );
      vdegree[u]--;
      vdegree[v]--;
    }
  }
}


bool g::is_edge( int u, int v ) const{
  return in_set( u, gA[v] );
}

int g::degree( int v ){
  return vdegree[v];
}

int g::min_degree(){
  int min = n+1;
  for( int i = 0; i < n; i++ ){
    if( vdegree[i] < min ){
      min = vdegree[i];
    }
  }
  return min;
}

vector<int> g::neighbors(int v){
  vector<int> neighborhood;
  for( int i = 0; i < n; i++ ){
    if( in_set(i,gA[v] )){
      neighborhood.push_back(i);
    }
  }
  return neighborhood;
}


void g::make_complement(){
  for( int i = 0; i < n; i++ ){
    gA[i]=set_complement( gA[i], n );
    set_delete( i, gA[i] );

  }
  numEdges = n*(n-1)/2 - numEdges;
}


void g::max_clique_backtrack( int l, int k ){
  // get new optimum value
  if( l > optSize ){
    optSize = l;
    optClique = X;
  }

  if( optSize == k ){
    return;
  }
  
  // calculates new tree level
  if( l == 0 ){
    gC[l] = gV;
  }else{
    // gets only the vertices we want
    gC[l] = gC[l-1] & gA[X[l-1]] & gB[X[l-1]];
  }
  vset cl = gC[l];

  // Bounding/Tree Pruning
  int m;
  m = l + set_order( gC[l] );
  //  cout << "gC[l] " << gC[l] << " " << set_order( gC[l] ) << endl;
  
  for( int x = 0; x < n; x++ ){
    //    cout << "x = " << x << endl;
    // bound
    if( m <= optSize ) return;
    if( in_set( x, cl )){
      X[l] = x;
      max_clique_backtrack( l + 1, k );
    }
      //    cl >>= 1;
  }
}


vector<int> g::max_clique( bool print, int k ){
  // set all global variables needed
  // (only use them for max-clique, so no use putting them in set_up)
  gC.resize( n + 1 );
  /* gV.resize( arraySize );
  for( int i = 0; i < arraySize; i++ ){
    gV[i] = empty;
  }*/
  gV = empty;
  X.resize( n );
  
  if( k < 1 ){
    k = n;
  }
 
  for( int i = 0; i < n+1; i++ ){
    /*  gC[i].resize( arraySize );
    for( int j = 0; j < arraySize; j++ ){
      gC[i][j] = empty;
      }*/
    gC[i] = empty;
  }
  for( int i = 0; i < n; i++ ){
    X[i] = 0;
    set_insert( i, gV );
  }
  optSize = 0;
  
  max_clique_backtrack( 0, k );
  if(print){
    cout << "Size: " << optSize << endl;
  }
  vector<int> realClique;
  for( int i = 0; i < optSize; i++ ){
    realClique.push_back( optClique[i] );
  }

  return realClique;
}


vector<int> g::max_independent_set( bool print, int k ){
  make_complement();
  vector<int> realClique = max_clique( print, k );
  make_complement();
  return realClique;
}



bool g::has_clique( int k, bool is ){
  bool hasClique = false;
  if( is ){
    make_complement();
  }
  max_clique( false, k );
  if( optSize >= k ){
    hasClique = true;
  }
  if( is ){
    make_complement();
  }
  return hasClique;
}


bool g::has_c( int c ){
  switch(c){
  case 4:{
    for( int i = 0; i < n; i++ ){
      for( int j = i+1; j < n; j++ ){
	if( is_edge( i, j ) ){
	  for( int k = j+1; k < n; k++ ){
	    if( is_edge( i, k ) ){
	      for( int l = 0; l<n; l++ ){
		if( l!=i &&is_edge( j,l ) && is_edge( k,l) ){
		  return true;
		}
	      }
	    }
	  }
	}
      }
    }
    return false;
  }
  default:
    cerr << "Error: Counting/Removing C" << c << " is not supported" << endl;
    return false;
  }
}


void g::make_complete(){
  for( int i = 0; i < n-1; i++ ){
    for( int j = i+1; j < n; j++ ){
      add_edge(i,j);
    }
  }
}


// Makes the graph a random graph from Erdos-Renyi algorithm
// @return the number of edges added
int g::make_rand_er( float sigma ){
  // set seed for random variable
  srand((unsigned)time(0));
  float r;

  for( int x = 0; x < n - 1; x++ ){
    for( int y = x + 1; y < n; y++ ){
     
      // Uses random value to determine if edge or not edge
      r = rand() / (float)RAND_MAX;
      if( r >= ( 1 - sigma )){
	add_edge( x, y );
      }
    }
  }
  return numEdges;
}


void g::read_g6( string g6 ){
  int numEntries = n*(n-1)/2;
  int c = 0;
  int c_count = 1;
  int i = 0;
  int j = 0;
  int b = 0;
  for( int i = 0; i < n; i++ ){
    for( int j = 0; j < i; j++ ){
      if( b == 0 ){
	c = g6[c_count];
	c_count++;
	c = c - 63;
      }
      if( c & ( 1 << ( 5 - b ) ) ){
	add_edge( i, j );
      }
      b = ( b + 1 ) % 6;
    }
  }
}


string g::to_g6(){
  const int powers[] = {32,16,8,4,2,1};
  string g6_string;
  char start = 0;

  if( n <= 62 ){
    start = n+63;
  }

  g6_string.push_back( start );

  int bitSize = n*(n-1)/2;
  int count = 0;
  int cur = 63;
  int pow = 0;

  for( int i = 0; i < n; i++ ){
    for( int j = 0; j < i; j++ ){
      pow = count % 6;
      if( is_edge( j, i ) ){
	cur = cur + powers[pow];
      }
      if( pow == 5 || j == n - 2 ){
	g6_string.push_back( cur );
	cur = 63; 
      }
      count++;
    }
  }
  
  return g6_string;
}

void g::print_g6( ostream *o ){
  *o << to_g6() << endl;
}


void g::get_tri_stats(){
  int edges_in_tri = 0;
  tri_stats = new int[4];
  for( int i = 0; i < 4; i++ ){
    tri_stats[i] = 0;
  }
  for( int i = 0; i < n - 2; i++ ){
    for (int j = i+1; j < n-1; j++ ){
      for( int k = j+1; k < n; k++ ){
	edges_in_tri = 0;
	if( is_edge( i, j ) ){
	  edges_in_tri++;
	}
	if( is_edge( i, k ) ){
	  edges_in_tri++;
	}
	if( is_edge( j, k ) ){
	  edges_in_tri++;
	}
	tri_stats[ edges_in_tri ]++;
      }
    }
  }
  trisCalced = true;
}

int g::get_num_tri_edges( int e ){
  if( !trisCalced ){
    get_tri_stats();
  }
  if( e >= 0 && e <= 3 ){
    return tri_stats[ e ];
  }
}


void g::remove_vs( vector<int> cuts, int k ){
  int x;
  int removed = 0;
  
  // ges the array of cuts and sorts them
  // int elements = sizeof(cuts) / sizeof(cuts[0]);
  sort(cuts.begin(),cuts.end());
  for( int i = 0; i < k; i++ ){
    cuts[i] = cuts[i] - i;
  }
  
  // remove each vertex from the graph
  for( int j = 0; j < k; j++ ){
    x = cuts[j];
    //cout << "Removing " << x << "..." << endl;
    for( int i = 0; i < n; i++ ){
      if( i != x ){
	set_cut( x, gA[i] );
	if( i > x ){
	  gA[i-1] = gA[i];
	  vdegree[i-1] = vdegree[i];
	}
      }
    }
    n--;
    gA[n] = empty;
  }
}


int g::get_p3s( int * tab, int p ){
  int count = 0;
  vector<int*> paths;
  
  for( int i = 0; i < p; i++ ){
    tab[i] = 1;
  }

  // matrix for endpoints of p3
  int endpoints[n][n];
  for( int i = 0; i < n - 1; i++ ){
    for( int j = 0; j < n; j++ ){
      endpoints[i][j] = 0;
      endpoints[j][i] = 0;
    }
  }

  for( int i = 0; i < n - 2; i++ ){
    for( int j = i+1; j < n - 1; j++ ){
      for( int k = j+1; k < n; k++ ){
	if(is_edge( i, j ) && is_edge( j, k )){
	  endpoints[i][k] = 1;
	  endpoints[k][i] = 1;
	}
	if(is_edge( j, i ) && is_edge( i, k )){
	  endpoints[j][k] = 1;
	  endpoints[k][j] = 1;
	}
	if(is_edge( i, k ) && is_edge( k, j )){
	  endpoints[i][j] = 1;
	  endpoints[j][i] = 1;
	}
      }
    }
  }

  // find all legit cones with no p3 endpoints
  for( uint64_t i = 0; i < p; i++ ){
    for( int x = 0; x < n-1 && tab[i]; x++ ){
      for( int y = 0; y < n && tab[i]; y++ ){
	if( endpoints[x][y] ){
	  if( in_set( x, i ) && in_set( y, i ) ){
	    tab[i] = 0;
	    count++;
	  }
	}
      }
    }
  }
  return p - count;
}


int g::get_p3s2( bool * tab, int p ){

  // vector<int*> paths;
  for( int i = 0; i < p; i++ ){
    tab[i] = 1;
  }

  bcur_tab = tab;

  // matrix for endpoints of p3
  int endpoints[n][n];
  for( int i = 0; i < n - 1; i++ ){
    for( int j = 0; j < n; j++ ){
      endpoints[i][j] = 0;
      endpoints[j][i] = 0;
    }
  }

  int s;

  for( int i = 0; i < n - 2; i++ ){
    for( int j = i+1; j < n - 1; j++ ){
      for( int k = j+1; k < n; k++ ){
	if(is_edge( i, j ) && is_edge( j, k )){
	  endpoints[i][k] = 1;
	  endpoints[k][i] = 1;
	  s = (( 0 | 1 << i ) | 1 << k);
	  bcur_tab[s] = 0;
	  recursive_p3( s, -1 );
	}
	if(is_edge( j, i ) && is_edge( i, k )){
	  endpoints[j][k] = 1;
	  endpoints[k][j] = 1;
	  s = (( 0 | 1 << k ) | 1 << j);
	  bcur_tab[s] = 0;
	  recursive_p3( s, -1 );
	}
	if(is_edge( i, k ) && is_edge( k, j )){
	  endpoints[i][j] = 1;
	  endpoints[j][i] = 1;
	  s= (( 0 | 1 << i ) | 1 << j);
	  bcur_tab[s] = 0;
	  recursive_p3( s, -1 );
	}
      }
    }
  }
  int count = 0;

  for( int i = 0; i < p; i++ ){
    if( !bcur_tab[i] ){
      count++;
    }
  }
  return p - count;
}


void g::recursive_p3( uint64_t cur_cl, int cur_v ){
  int x;
  bool good = false;
  for( int i = cur_v + 1; i < n; i++ ){
    if( !in_set( i, cur_cl ) ){
      x = i;
      good = true;
      break;
    }
  }
  if( good ){
    recursive_p3( cur_cl, x );
    set_insert( x, cur_cl );
  
    if( bcur_tab[ cur_cl ] == 1 ){
      bcur_tab[ cur_cl ] = 0;
      recursive_p3( cur_cl, x );
    }
    else
      return;
  }
  else{
    return;
  }
}


// get all vertices that are adjacent to v \in S
void g::get_closures( uint64_t * tab, int p ){
  for( uint64_t s = 0; s < p; s++ ){
    tab[s] = 0;
    for( uint64_t v = 0; v < n; v++ ){
      if( s & gA[v] ){
	tab[s] = tab[s] | (1 << v );
      }
    }
  }
}

void g::get_closures2( uint32_t * tab, int p ){
  cur_utab = tab;
  for( uint32_t s = 0; s < p; s++ ){
    cur_utab[s] = 0;
  }
  for( int i = 0; i < n; i++ ){
    cur_utab[ 0 | 1 << i ] = gA[i];
  }
  
  for( uint32_t s = 3; s < p; s++ ){
    recursive_clos( s );
  }

  cur_utab[ p - 1 ] = p-1;
}

int g::recursive_clos( uint32_t s ){
  if( cur_utab[s] == 0 && s != 0 && s != (( 0 | 1 << n ) - 1)){
    int c = 0;
    uint32_t sc = s;
    uint32_t ss = s;
    while( !(sc & 1) ){
      c++;
      sc = sc >> 1;
    }
    set_delete( c, ss );
    cur_utab[s] = (cur_utab[ 0 | 1 << c ] | recursive_clos( ss ));// & ( (( 1 << n ) - 1) & ~s);
  }
  
  return cur_utab[s];

  //  return cur_tab[s];
}




bool g::glue_graphs( int d, int num_e, g * y, vector<int>cones ){
  // add all edges of original neighborhood
  for( int i = 1; i <= d; i++ ){
    add_edge( 0, i );
  }

  // add all edges of graph X
  int cur = 1;
  for( int i = 0; i < num_e; i++ ){
    add_edge( cur, cur + 1 );
    cur += 2;
  }

  // add edges in between cones
  for( int i = 0; i < d; i++ ){
    int c = cones[i];
    for( int j = 0; j < y->order(); j++ ){
      if( in_set( j, c ) ){
	add_edge( 1 + i, 1 + d + j );
      }
    }
  }

  // add edges of Y
  int diff = d + 1;
  for( int i = 0; i < y->order(); i++ ){
    gA[ i + diff ] |= (y->get_adj(i) << diff);
    vdegree[i + diff] += y->degree(i);
  }

  return true;
}


bool g::v_extend( int d, g * y, vset cone ){
  // add all edges of original neighborhood
  

  // add edges in between cones
  for( int j = 0; j < y->order(); j++ ){
    if( in_set( j, cone ) ){
      add_edge( 0, 1 + j );
    }
  }

  // add edges of Y
  int diff = 1;
  for( int i = 0; i < y->order(); i++ ){
    gA[ i + diff ] |= (y->get_adj(i) << diff);
    vdegree[i + diff] += y->degree(i);
  }

  return true;
}


bool g::glue_graphs( g * x, g * y, vector<int> cones, int p ){

  for( int i = 1; i < x->order() + 1; i++ ){
    add_edge( 0, i );
  }
  
  int diff = 1;

  for( int i = 0; i < x->order() - 1; i++ ){
    for( int j = i + 1; j < x->order(); j++ ){
      if( x->is_edge(i, j) ){
	add_edge( diff + i, diff + j );
      }
    }
  }
  
  diff = diff + x->order();

  for( int i = 0; i < y->order() - 1; i++ ){
    for( int j = i + 1; j < y->order(); j++ ){
      if( y->is_edge(i, j) ){
	add_edge( diff + i, diff + j );
      }
    }
  }  

  // add edges in between cones
  for( int i = 0; i < x->order(); i++ ){
    int c = cones[i];
    for( int j = 0; j < y->order(); j++ ){
      if( in_set( j, c ) ){
	add_edge( 1 + i, 1 + x->order() + j );
      }
    }
  }

  return true;
}

