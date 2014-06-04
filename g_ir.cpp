#include "g.h"

vset g::private_neighbors( int v, vset S ){
  vset pn = gA[v] | (shifter << v);
  int count = 0;
  for( int i = 0; i < n; i++ ){
    if( in_set(i,S) && i != v ){
      count++;
      pn = pn & ~(gA[i] | (shifter << i));
    }
  }
  return pn;
}

int * ir_num;
int max_possible;
vector<vset> * irs;
int max_ir;

bool g::not_isolate( int v, vset S ){
  return gA[v] & S;
}

/* algorithm for enumerating all irredundant sets
*/
void g::all_ir_bt( vset cur_ir, int l, int v ){
  if( l > 0 ){
    ir_num[l-1]++;
    irs[l-1].push_back(cur_ir);
    if( l > max_ir )
      max_ir = l;
  }

  if( l == max_possible) return;

  // new tree level
  vset candidates = 0;
  if( l == 0 )
    candidates = gV;
  else{
    for( int i = v + 1; i < n; i++){
      if( in_set(i,gV) ){
	set_insert(i,cur_ir);
	bool good = true;
	for( int j = 0; j < n; j++ ){
	  if( in_set(j,cur_ir)){
	    if( private_neighbors(j,cur_ir) == 0){
	      good = false;
	      break;
	    }
	  }
	}
	set_delete(i,cur_ir);
	if( good )
	  set_insert(i,candidates);
      }
    }
  }

  for( int x = 0; x < n; x++){
    if( in_set(x,candidates))
      all_ir_bt( ( cur_ir | 1 << x ), l + 1, x );
  }
}

void g::all_irs(int max_ir_possible){
  max_possible = max_ir_possible;
  

  int k = max_possible;
  irs = new vector<vset>[k];
  ir_num = new int[k];
  for( int i = 0; i < k; i++ )
    ir_num[i] = 0;

  
  all_ir_bt(0,0,0);
}

int g::max_irs(int max_ir_possible ){
  gV = (1 << n) - 1;
  max_ir = 0;
  all_irs(max_ir_possible);

  delete [] irs;
  delete [] ir_num;

  return max_ir;
}

int g::max_irs(int max_ir_possible, vset of_these){
  gV = of_these;
  max_ir = 0;
  all_irs(max_ir_possible);
  return max_ir;
}

int g::max_irs(int max_ir_possible, vector<vset>& max_irs){
  max_ir = 0;
  all_irs(max_ir_possible);
  for( vector<vset>::iterator it = irs[max_ir-1].begin(); 
    it != irs[max_ir-1].end(); it++){
    max_irs.push_back(*it);
  }
  
  delete [] irs;
  delete [] ir_num;

  return max_ir;
}


int * ir_tab;
int g::get_irs(int * tab, vector<vset>& max_irs, int max_ir_possible){
  max_ir = 0;
  all_irs(max_ir_possible);

  int p = 1 << n;
  ir_tab = new int[p];



  int ks = max_ir_possible;
  int ir_n;
  vset ir_set;

  for( int i = 0; i < p; i++ ){
    tab[i] = 0;
    ir_tab[i] = 0;
  }

  for( vector<vset>::iterator it = irs[max_ir-1].begin(); 
    it != irs[max_ir-1].end(); it++){
    max_irs.push_back(*it);
  }

  for( int k = ks; k > 0; k-- ){
    ir_n = ir_num[k-1];
    //    cout << cl_num << endl;
    for( int i = 0; i < ir_n; i++ ){
      ir_set = irs[k-1][i];
      ir_tab[ ir_set ] = k;
      recursive_ir( ir_set, k, -1 );
      //     cout << endl;
    }
  }


  for( int i = 0; i < p; i++ ){
    int s = (( 1 << n ) - 1) & (~i);
    tab[s] = ir_tab[i];
  }

  delete [] irs;
  delete [] ir_num;
  delete [] ir_tab;

  return max_ir;
}

void g::recursive_ir( uint64_t cur_cl, int k, int cur_v ){
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
    recursive_ir( cur_cl, k, x );
    set_insert( x, cur_cl );
    if( ir_tab[ cur_cl ] == 0 ){
      ir_tab[ cur_cl ] = k;
      recursive_ir( cur_cl, k, x );
    }
    else
      return;
  }
  else{
    return;
  }
}

int g::get_independent_sets( int * tab, int max_is ){
  //is_sets;
  int p = 1 << n;
  all_independent_sets( max_is, p );
  int cl_num;
  vset clique;

  for( int k = 0; k < max_is; k++ ){
    cl_num = clique_num[k];
    for( int cl = 0; cl < cl_num; cl++ ){
      clique = cliques[k][cl];
      vec_is_sets.push_back(clique);
    }
  }

  int k = vec_is_sets.size();
  /*tab_is = new int[k];
  for( int i = 0; i < k; i++ ){
    tab_is[i] = is_sets[i];
  }*/

  delete [] cliques;
  delete [] clique_num;

  return k;

}


vset cur_p5;
int g::get_p5s( bool * tab, int p ){

  // vector<int*> paths;
  for( int i = 0; i < p; i++ ){
    tab[i] = 1;
  }

  bcur_tab = tab;

  // matrix for endpoints of p5
  int endpoints[n][n];
  for( int i = 0; i < n - 1; i++ ){
    for( int j = 0; j < n; j++ ){
      endpoints[i][j] = 0;
      endpoints[j][i] = 0;
    }
  }

  vset s;

  int enda;
  int endb;
  int mid;

  for( int i = 0; i < n - 2; i++ ){
    for( int j = i+1; j < n - 1; j++ ){
      for( int k = j+1; k < n; k++ ){
	bool is_p3i = true;
	if(is_edge( i, j ) && is_edge( j, k ) && !is_edge(i,k)){
	  enda = i; mid = j; endb = k;
	}
	else if(is_edge( i, j ) && is_edge( k, i ) && !is_edge(j,k)){
	  enda = j; mid = i; endb = k;
	}
	else if(is_edge( i, k ) && is_edge( k, j ) && !is_edge(i,j)){
	  enda = i; mid = k; endb = j;
	}
	else is_p3i = false;
	if( is_p3i ){
	  vset p5_a_ends = gA[enda] & ~(gA[mid] | 1 << mid) & ~(gA[endb] | 1 << endb);
	  vset p5_b_ends = gA[endb] & ~(gA[mid] | 1 << mid) & ~(gA[enda] | 1 << enda);
	  for( int a = 0; a < n; a++ ){
	    if( in_set(a,p5_a_ends) ){
	      for( int b = 0; b < n; b++ ){
		if( in_set(b,p5_b_ends) ){
		  if( !is_edge(a,b) && a!=b ){
		    endpoints[a][b] = 1;
		    endpoints[b][a] = 1;
		    s = (( 0 | 1 << a ) | 1 << b);
		    bcur_tab[s] = 0;
		    cur_p5 = ((( s | 1 << i ) | 1 << j ) | 1 << k);
		    recursive_p5(s,-1);
		  }
		}
	      }
	    }
	  }
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


void g::recursive_p5( vset cur_cl, int cur_v ){
  int x;
  bool good = false;
  for( int i = cur_v + 1; i < n; i++ ){
    if( !in_set( i, cur_cl ) && !in_set( i, cur_p5) ){
      x = i;
      good = true;
      break;
    }
  }
  if( good ){
    recursive_p5( cur_cl, x );
    set_insert( x, cur_cl );
    bcur_tab[ cur_cl ] = 0;
    recursive_p5( cur_cl, x );
  }
  else{
    return;
  }
}
