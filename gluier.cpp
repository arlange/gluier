#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <list>
#include <map>
#include <limits>
#include <string>

#include "g.h"
#include "/home/alex/software/nauty24r2/gtools.h"

using namespace std;

bool * p5_tab;
uint8_t * tab1;
uint32_t * tab2;
int num_cones;

int ** valid_cones;
int * tnum_cones;

int * cone_size_counts;
int * max_ir_counts;

int ir, d, xe, e_min, e_max, y_p;
g * y_g;
bool check_all;
int y_order;

list<string> * g_list;
map<string,int> * g_map;

int g_count;

double avg_time, min_time, max_time, cur_time;
double total_time;
double avg_numg, min_numg, max_numg, cur_numg;
double avg_conesize, min_conesize, max_conesize, cur_conesize;

char opt;

///// OPTIONS /////

bool mixed;  // m
bool fix_degree;  
bool just_count;  // u
bool std_not_file; // s
bool no_log; // l
bool display_updates; //p

bool func_v_extend; // v
bool func_drop_v; // d
bool func_filter; // f
bool func_glue; // g

bool use_map; // a

bool store_and_sort;

//////////////////

int filter_num_ir;
int filter_num_tri;
int filter_num_c6;


int g_out;

const double DBL_MAX = numeric_limits<double>::max( );


/* converts a string of a graph to its g6 canonical labeling*/
string convert_gstring( string g_string ){
  int m = 1;
  int len = g_string.length();
  char * cstr = new char[len];
  for( int i = 0; i < len; i++ ){
    cstr[i] = g_string[i];
  }
  int n = graphsize( cstr );
  graph inG[m*n];
  graph outG[m*n];
  stringtograph( cstr, inG, m );
  fcanonise(inG, m, n, outG, NULL, 0 );
  char * out_cstr = ntog6( outG, m, n );
  string can_g( out_cstr );
  delete [] cstr;
  return can_g.substr(0,can_g.size()-1);
}

void add_to_map( string g_string ){
  g_map->insert(make_pair(convert_gstring(g_string),1));
}


void print_from_map( ostream * out ){
  g_out = 0;
  for( map<string,int>::iterator it = g_map->begin();
       it != g_map->end(); it++ ){
    *out << it->first << endl;
    g_out++;
  }
}

void add_to_list( string g_string ){
  g_list->push_back( convert_gstring( g_string ) );
}


/* adds converted string to the list */
void add_graph( string g_string ){
  if( use_map )
    add_to_map( g_string );
  else
    add_to_list( g_string );
}


/* sorts the list and counts the unique strings */
void count_graphs( ){
  g_out = 0;
  if( !g_list->empty() ){
    g_list->sort();
    string prev = g_list->front();
    g_out++;

    list<string>::iterator it = g_list->begin();
    it++;

    for( it; it != g_list->end(); it++ ){
      if( it->compare( prev ) != 0  ){
	g_out++;
      }
      prev = *it;
    }
  }
}

/* sorts the list and prints the unique strings */
void print_graphs( ostream * out ){
  if( use_map ){
    print_from_map( out );
  }
  else{
    g_out = 0;
    if( !g_list->empty() ){
      g_list->sort();
      string prev = g_list->front();
      *out << prev << endl;
      g_out++;
      
      list<string>::iterator it = g_list->begin();
      it++;
      
      for( it; it != g_list->end(); it++ ){
	if( it->compare( prev ) != 0  ){
	  g_out++;
	  *out << *it << endl;
	}
	prev = *it;
      }
    }
  }
}

/* prints the log, based on global data previously collected*/
void print_log( ostream * out ){
  *out << "Total in = " << g_count << endl;
  *out << "Total out = " << g_out << endl;
  *out << endl;
  *out << "Time" << endl;
  *out << "Total = " << total_time << endl;
  *out << "Avg = " << avg_time << endl;
  *out << "Min = " << min_time << endl;
  *out << "Max = " << max_time << endl;
  *out << endl;

  if( func_glue ){
    *out << "Number of Graphs Glued" << endl;
    *out << "Avg = " << avg_numg << endl;
    *out << "Min = " << min_numg << endl;
    *out << "Max = " << max_numg << endl;
    *out << endl;
    *out << "Number of Cones per Glue" << endl;
    *out << "Avg = " << avg_conesize << endl;
    *out << "Min = " << min_conesize << endl;
    *out << "Max = " << max_conesize << endl;
  }

  if( func_filter ){
    *out << endl;
    *out << "*** Filter Stats ***" << endl;
    *out << "Number with bad IR = " <<  filter_num_ir << endl;
    *out << "Number with triangle = " <<  filter_num_tri << endl;
    if( !mixed ) *out << "Number with induced C6 = " <<  filter_num_c6 << endl;
  }

  if( opt == 'v' ){
    *out << endl;
    *out << "Cone counts:" << endl;
    for( int i = 0; i <= ir; i++ ){
      *out << i << ": " << cone_size_counts[i] << endl;
    }

    *out << endl;
    *out << "Max IR counts:" << endl;
    for( int i = 0; i <= ir; i++ ){
      if( max_ir_counts[i] > 0 )
	*out << i << ": " << max_ir_counts[i] << endl;
    }
  }

}





/* extends a graph by one vertex, making sure 
   that the neighborhood (cone) of the new vertex:
     1. does not cause a triangle
     2. does not create an irredundant set of order IR
*/
void v_extend(){
  int * cur_ir_vs;
  vset cur_ir_set;
  int maxir;
  int p = 1 << y_order;

  // to avoid triangles, neighborhood of v has to be an independent set
  int * is_sets;
  int k = y_g->get_independent_sets( is_sets, y_order );

  // get all max irredundant sets
  int * ir_tab = new int[p];
  vector<vset> max_irs;
  maxir = y_g->get_irs( ir_tab, max_irs, ir + 1 );
  max_ir_counts[maxir]++;

  vset s;

  // if the max irredundant set is less than ir-1, then adding one vertex does not 
  // yield an irredundant set of order ir
  bool can_cause_ir = (maxir >= ir-1);
  if( !can_cause_ir ){
    g good(y_order + 1);
    good.v_extend( y_order, y_g, empty );
    add_graph(good.to_g6());
  }

  // loop through all independent sets -- each one is a possible cone
  for( int i = 0; i < k; i++ ){
    bool good_cone = true;
    s = y_g->vec_is_sets[i];

    if( p5_tab[s] || mixed ){
 
      if( can_cause_ir ){

	// check all irredundant sets
	for( vector<vset>::iterator it = max_irs.begin();
	     it != max_irs.end(); it++){
	  bool killed_pn = false;
	  
	  cur_ir_set = *it;
 
	  // diff is the difference between s (the IS/cone) and things not in the IRS
	  vset diff = s & ~(cur_ir_set);
	  // similar definition as diff
	  vset same = s & (cur_ir_set);
	  
	  // check if the cone "hits" the IR
	  bool hits_ir = false;
	  if( same != empty ){
	    hits_ir = true;
	  }
	  //bool hits_ir = (s & (*it)); //(does same thing)
	  
	  // if it does hit, we are done
	  if( hits_ir ){
	    killed_pn = true;
	  }
	  else{
	    // get vertices of current irredundant set
	    cur_ir_vs = new int[maxir];
	    vset cur_ir_pns[maxir];
	    
	    int cur_i = 0;
	    // get actual vertices of the irredundant set,
	    // and each one of their private neighborhoods
	    for( int j = 0; j < y_order; j++){
	      if( in_set(j, cur_ir_set )){
		cur_ir_vs[cur_i] = j;
		cur_ir_pns[cur_i] = y_g->private_neighbors(j,cur_ir_set);
		cur_i ++;
	      }
	    }
	    
	    // check to make sure v won't have any private neighbors
	    bool v_not_allowed = false;
	    for( int j = 0; j < maxir; j++ ){
	      if( (cur_ir_pns[j] | s) == s ){
		v_not_allowed = true;
		break;
	      }
	    }
	    
	    /* if v can't be in the new irredundant set, then some other
	       vertex of the PN that is not in the IR set might be able to use 
	       v has a PN */
	    bool s_okay = true;
	    if( v_not_allowed ){
	      for( int t = 0; t < y_order; t++ ){
		if( in_set(t,s) ){
		  vset s_copy = cur_ir_set;
		  set_insert(t,s_copy);
		  bool t_killed_one = false;
		  for( int k = 0; k < maxir; k++ ){
		    if( y_g->private_neighbors( cur_ir_vs[k], s_copy) == empty ){
		      t_killed_one = true;
		      break;
		    }
		  }
		  
		  if( !t_killed_one ){
		    s_okay = false;
		    break;
		  }
		}
	      }
	    }
	    else s_okay = false;
	    
	    if( s_okay ) killed_pn = true;
	    
	    delete [] cur_ir_vs;
	  }
	  if( !killed_pn ){
	    good_cone = false;
	    break;
	  }
	  else{
	    
	  }
	}
      }
      if( good_cone ){
	g good(y_order + 1);
	good.v_extend( y_order, y_g, s );
	cone_size_counts[set_order(s)]++;
	add_graph(good.to_g6());
      }
    }
  }
  
  delete [] ir_tab;
}



void glue(int d){


}


/* 
   Remove each vertex to create a set of (3,ir)-graphs that
   have one less vertex
 */
void drop_v(){
  int num_vs = y_g->order();
  cerr << "Cone size: " << endl;
  for( int i = 0; i < num_vs; i++ ){
    cerr << y_g->degree(i) << endl;
    vector< int > cut;
    cut.push_back(i);
    g y_copy = *y_g;
    y_copy.remove_vs( cut, 1 );
    add_graph( y_copy.to_g6() );
  }
}


/*
  Check graph to make sure it has no IR of order ir,
  as well as no IR set of order 3 in the complement
*/
void filter_g(){
  bool good = true;
  bool go_to_s = !mixed;

  int n = y_g->order();

  // check if ir < 
  int max_ir = y_g->max_irs(ir + 1);
 
  if( max_ir < ir ){
    // check for any triangles
    for( int a = 0; a < y_g->order(); a++ ){
      for( int b = a+1; b < y_g->order(); b++ ){
	if( y_g->is_edge(a,b) ){
	  for( int c = b+1; c < y_g->order(); c++ ){
	    if( y_g->is_edge(a,c) && y_g->is_edge(b,c) ){
	      good = false;
	      filter_num_tri++;
	      break;
	    }
	  }
	}
      }
    }
    // if irredundant (go_to_s), check for induced C6
    if( good && go_to_s ){
      y_p = 1 << y_order;
      p5_tab = new bool[y_p];
      y_g->get_p5s(p5_tab, y_p);
      for( int i = 0; i < n; i++){
	if(!p5_tab[y_g->get_adj(i)]){
	  good = false;
	  filter_num_c6++;
	  break;
	}
      }

      delete [] p5_tab;
    }
    
    if( good ){
      add_graph(y_g->to_g6());
    }
  }
  else{
   filter_num_ir++;
  }
}



int main( int argc, char *argv[] ){
  clock_t start, stop;
  clock_t total_start, total_stop;
  bool all = false;

  total_start = clock();

  // possible different options
  mixed = false;           // m
  fix_degree = false;      
  just_count = false;      // u
  std_not_file = false;    // s
  no_log = true;           // l
  func_v_extend = false;   // v
  func_drop_v = false;     // d
  func_filter = false;     // f
  func_glue = false;       // g
  display_updates = false; // p
  use_map = false;         // a
  store_and_sort = true;

  bool bad_args = false;

  string in_file, out_file, log_file;
  char * opts;
  char cur_opt;
  int degree;

  // get all the options from the command line and
  // set the flags accordingly
  bool get_more_opts = true;
  int cur_opt_i = 0;
  while( get_more_opts ){
    cur_opt_i++;
    opts = argv[cur_opt_i];
    if( opts[0] == '-' ){
      opts++;
      while( *opts != '\0' && get_more_opts ){
	switch(*opts){
	case 'm': mixed = true; break;
	case 'u': just_count = true; break;
	  //case 's': std_not_file = true; break;
	case 'l': no_log = false; break;
	case 'v': func_v_extend = true; break;
	case 'd': func_drop_v = true; break;
	case 'f': func_filter = true; break;
	case 'g': func_glue = true; break;
	case 'p': display_updates = true; break;
	case 'o': store_and_sort = false; break;
	case 'a': use_map = true; break;
	default: 
	  cerr << "Error: -" << *opts << " is not an option." << endl;
	  bad_args = true;
	  get_more_opts = false;
	}
	opts++;
      }
    }
    else{
      get_more_opts = false;
    }
    
    if( get_more_opts && ((cur_opt_i + 1 == argc) 
			  || (cur_opt_i + 2 == argc)) ){ 
      // didn't put anything except opts
      cerr << "Error: Need arguments for IR and in_file" << endl;
      get_more_opts = false;
      bad_args = true;
    }
  }

  g_list = new list<string>();
  g_map = new map<string,int>();

  int num_funcs = func_v_extend + func_drop_v + func_filter + func_glue;
  if( num_funcs > 1 && !bad_args ){
    cerr << "Error: -vdf are incompatible." << endl;
    bad_args = true;
  }
  else if( num_funcs < 1 ){
    func_v_extend = true;  // default to vertex extend
  }

  if( !bad_args ){
    ir = atoi( argv[ cur_opt_i ] );
    in_file = argv[ cur_opt_i + 1];
    cur_opt_i += 2;
    if( !just_count ){
      if( cur_opt_i != argc ){
	out_file = argv[cur_opt_i];
	cur_opt_i++;
      }else{
	bad_args = true;
	cerr << "Error: need output filename" << endl;
      }
    }
    if( !no_log && !bad_args){
      if( cur_opt_i != argc )
	log_file = argv[cur_opt_i];
      else{
	bad_args = true;
	cerr << "Error: need log filename" << endl;
      }
    }
  }

  if( bad_args ){
    cerr << "BAD ARGS" << endl;
    return 0;
  }

  if( mixed ){
    cerr << "Mixed Ramsey Numbers" << endl;
  }else{
    cerr << "Irredundant Ramsey Numbers" << endl;
  }

  /////// done setting up options and needed files

  // set up log variables
  avg_time = 0; min_time = DBL_MAX; max_time = 0;
  avg_numg = 0; min_numg = DBL_MAX; max_numg = 0;
  avg_conesize = 0; min_conesize = DBL_MAX; max_conesize = 0;

  cone_size_counts = new int[ir+1];
  max_ir_counts = new int[ir+1];
  for( int i = 0; i <= ir; i++ ) cone_size_counts[i] = 0;
  for( int i = 0; i <= ir; i++ ) max_ir_counts[i] = 0;

  check_all = true;
    
  vector<string> y_graphs;
  string g_string;


  // get input/output files
  ifstream ifs ( in_file.c_str() );
  if( !ifs.is_open() ){
    cerr << "Error opening file " + in_file << endl;
    return 0;
  }

  ofstream ofs, log;

  if( !just_count ){
    ofs.open( out_file.c_str() );
    if( !ofs.is_open() ){
      cerr << "Error opening " << out_file << endl;
      return 0;
    }
  }

  if( !no_log ){
    log.open( log_file.c_str() );
    if( !log.is_open() ){
      cerr << "Error opening " << out_file << endl;
      return 0;
    } //  */
  }
  // done with files, let's start!!!


  cerr << "****************" << endl;
  cerr << "ir = " << ir << endl;

  filter_num_ir = 0; filter_num_tri = 0; filter_num_c6 = 0;
  
  g_count = 0;

  // read through each g6 string of input file
  while( getline( ifs, g_string ) ){ 
    g_count++;

    start = clock();
    
    g y( g_string[0] - 63 );
    y.read_g6( g_string );
    y_order = y.order();

    // reset currents
    cur_time = 0; cur_numg = 0; cur_conesize = 0;

    if( func_filter ){
      y_g = &y;
      filter_g();
    }

    else if( func_v_extend ){
      y_p = 1 << y_order;
      p5_tab = new bool[y_p];
      y_g = &y;
      y_g->get_p5s(p5_tab, y_p);
      v_extend();
      delete [] p5_tab;
    }
    else if( func_drop_v ){
      y_g = &y;
      drop_v();
    }

    stop = clock();

    cur_time = ((float)(stop - start))/((float)CLOCKS_PER_SEC);
    avg_time += cur_time;
    avg_numg += cur_numg;
    avg_conesize += cur_conesize;
    if( cur_time < min_time )
      min_time = cur_time;
    if( cur_time > max_time )
      max_time = cur_time;
    if( cur_numg < min_numg )
      min_numg = cur_numg;
    if( cur_numg > max_numg )
      max_numg = cur_numg;
   
    if( display_updates ){
      if( ( g_count - 1 ) % 1000 == 0 )
	cerr << g_count << ": " << cur_time << endl;
    }
    //cerr << ((float)(stop - start))/((float)CLOCKS_PER_SEC) << endl;
 
  }

  avg_time = avg_time / g_count;
  avg_numg = avg_numg / g_count;
  avg_conesize = avg_conesize / g_count;
    
  if( just_count )
    count_graphs();
  else
    print_graphs( &ofs );

  total_stop = clock();
  total_time = ((float)(total_stop - total_start))/((float)CLOCKS_PER_SEC);

  if( !no_log )
    print_log( &log );

  ifs.close();
  if( !just_count ) ofs.close();
  if( !no_log ) log.close();

  delete [] cone_size_counts;
  delete [] max_ir_counts;

  delete g_list;
  delete g_map;

  cerr << endl;
  cerr << "Graph count = " << g_out << endl;
  cerr << endl;
  return 0;
}
