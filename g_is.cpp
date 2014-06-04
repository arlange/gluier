#include "g.h"



struct g::set_list{
  int v;
  set_list * next;
  set_list ( int x, set_list* n ){
    v = x;
    next = n;
  }
};

void g::all_cliques_backtrack( uint64_t cur_cl, int l, int k, int v ){

  if( l > 0 ){
    clique_num[l-1]++;
    cliques[l-1].push_back( cur_cl );
  }

  if( l == k ){
    return;
  }
  
  // calculates new tree level
  if( l == 0 ){
    gC[l] = gV;
  }else{
    // gets only the vertices we want
    gC[l] = gC[l-1] & gA[v] & gB[v];
  }
  vset cl = gC[l];
  for( int x = 0; x < n; x++ ){
    if( cl & 1 ){
      all_cliques_backtrack( ( cur_cl | 1 << x ), l + 1, k, x );
    }
    cl >>= 1;
  }
}



void g::all_cliques( int k, int p ){
  // set all global variables needed
  // (only use them for max-clique, so no use putting them in set_up)
  gC.resize( n + 1 );
  /*gV.resize( arraySize );
  for( int i = 0; i < arraySize; i++ ){
    gV[i] = empty;
    }*/
  gV = empty;
  
  if( k < 1 ){
    k = n;
  }
 
  for( int i = 0; i < n+1; i++ ){
    /*gC[i].resize( arraySize );
    for( int j = 0; j < arraySize; j++ ){
      gC[i][j] = empty;
      }*/
    gC[i] = empty;
  }
  
  gV = p - 1;
  cliques = new vector<uint64_t>[ k ];
  clique_num = new int[k];
  for( int i = 0; i < k; i++ ){
    clique_num[i] = 0;
  }
  all_cliques_backtrack( 0, 0, k, 0 );
}

void g::all_independent_sets( int k, int p ){
  make_complement();
  all_cliques( k, p );
  make_complement();
}

void g::get_independences( int * tab, int p, int max_is){

  all_independent_sets( max_is, p );
  /*  for( int i = 0; i < max_is; i++ ){
    cout << cliques[i].size() << endl;
  }*/

  for( uint64_t i = 0; i < p; i++ ){
    int k = max_is;
    uint64_t s = (~i) & ( ( 1 << n ) - 1 );
    tab[i] = k;
    bool no_is = true;
    while( no_is && k > 0 ){
      //for( vector< uint64_t >::iterator it = cliques[k-1].begin();
      //   it != cliques[k-1].end() && no_is; it++ ){
	//	cerr << s << " " << *it << endl;
      if( set_order(s) >= k ){
      for( int cl = 0; cl < clique_num[k-1]; cl++ ){
	if( ( s | cliques[k-1][cl] ) == s ){
	  no_is = false;
	}
      }
      }
      if( no_is ){
	tab[i]--;
	k--;
      }
    }
  }

  delete [] cliques;
  delete [] clique_num;
}




void g::get_independences2( int * tab, int p, int max_is){
  all_independent_sets( max_is, p );
  
  bool no_is = true;
  int ks = max_is;

  int o;
  tab[ p - 1 ] = 0;
  set_list * slist = new set_list( -1, NULL );
  for( int i = 0; i < p-1; i++ ){
    tab[i] = 1;
    slist = new set_list( i, slist );
  }

  //  set_list * slist;
  set_list * head = slist;
  set_list * prev = head;
  uint64_t clique;

  int i, s, cl_num;

  for( int k = ks; k > 1; k-- ){
    cl_num = clique_num[k-1];
    // cerr << " k = " << k << ", num = " << cl_num << endl;

    for( int cl = 0; cl < cl_num; cl++ ){
      clique = cliques[k-1][cl];

      // start at beginning
      slist = head;

      // traverse the whole list
      while( slist != NULL ){
	i = slist->v;
	s = (~i) & ( ( 1 << n ) - 1);
	if( ( s | clique ) == s ){
	  if( slist == head ){
	    head = slist->next;
	    delete slist;
	    slist = head;
	  }
	  else{
	    prev->next = slist->next;
	    delete slist;
	    slist = prev->next;
	  }
	  tab[i] = k;
	}
	else{
	  prev = slist;
	  slist = slist->next;
        }
      }
      if( head == NULL )
	break;
    }
    if( head == NULL )
      break;
    //  }
  }

  slist = head;
  while( slist != NULL ){
    prev = slist->next;
    delete slist;
    slist = prev;
  }

  delete [] cliques;
  delete [] clique_num;
}


void g::get_independences4( uint8_t * tab, int p, int max_is){
  int ks = max_is;
  cur_tab = new uint8_t[p];
  all_independent_sets( max_is, p );

  for( int i = 0; i < p; i++ ){
    tab[i] = 0;
    cur_tab[i] = 0;
  }

  int cl_num;
  uint64_t clique;

  for( int k = ks; k > 0; k-- ){
    //    cout << k << " ";
    cl_num = clique_num[k-1];
    //    cout << cl_num << endl;
    for( int cl = 0; cl < cl_num; cl++ ){
      clique = cliques[k-1][cl];
      
      cur_tab[ clique ] = k;
      recursive_is( clique, k, -1 );
      //     cout << endl;
    }
  }

  for( int i = 0; i < p; i++ ){
    int s = (( 1 << n ) - 1) & (~i);
    tab[s] = cur_tab[i];
  }

  delete [] cliques;
  delete [] clique_num;
  delete [] cur_tab;

}


void g::recursive_is( uint64_t cur_cl, int k, int cur_v ){
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
    recursive_is( cur_cl, k, x );
    set_insert( x, cur_cl );
    if( cur_tab[ cur_cl ] == 0 ){
      cur_tab[ cur_cl ] = k;
      recursive_is( cur_cl, k, x );
    }
    else
      return;
  }
  else{
    return;
  }
}



void g::get_independences3( int * tab, int p, int max_is){
  all_independent_sets( max_is, p );
  
  bool no_is = true;
  int ks = max_is;

  // make lists of each order
  set_list * slists[ks];
  for( int i = 0; i < ks; i++ ){
    slists[i] = new set_list( -1, NULL );
  }

  int o;
  tab[ p - 1 ] = 0;
  // set_list * slist = new set_list( -1, NULL );
  for( int i = 0; i < p-1; i++ ){
    tab[i] = 1;
    //  slist = new set_list( i, slist );
    int s = (~i) & ( ( 1 << n ) - 1);
    o = set_order(s);
    if( o > ks ){
      o = ks;
    }
    slists[o-1] = new set_list( i, slists[o-1] );
  }

  set_list * slist;
  set_list * head;// = slist;
  set_list * prev;// = head;

  int i, s, cl_num;

  for( int k_1 = max_is; k_1 >= 2; k_1-- ){
    head = slists[k_1-1];
    prev = head;

    for( int k = k_1; k > 1; k-- ){
    bool no_is = true;
    cl_num = clique_num[k-1];
    for( int cl = 0; cl < cl_num; cl++ ){
      uint64_t clique = cliques[k-1][cl];

      // start at beginning
      slist = head;

      // traverse the whole list
      while( slist != NULL ){
	i = slist->v;
	s = (~i) & ( ( 1 << n ) - 1);
	if( ( s | clique ) == s ){
	  if( slist == head ){
	    head = slist->next;
	    delete slist;
	    slist = head;
	  }
	  else{
	    prev->next = slist->next;
	    delete slist;
	    slist = prev->next;
	  }
	  tab[i] = k;
	}
	else{
	  prev = slist;
	  slist = slist->next;
        }
      }
      if( head == NULL ){
	break;
      }
    }
    if( head == NULL ){
      break;
    }
      }
  }

  delete [] cliques;
  delete [] clique_num;
}
