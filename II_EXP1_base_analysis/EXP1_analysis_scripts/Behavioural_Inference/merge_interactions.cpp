//to compile with option -std=c++11 


using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <vector>

struct collision_list {
  NumericVector ant1 ;
  NumericVector ant2 ;
  StringVector type ;
  NumericVector frame  ;
  NumericVector time ;
  NumericVector space    ;
  NumericVector x1     ;
  NumericVector x2      ;
  NumericVector y1       ;
  NumericVector y2        ;
  NumericVector angle1     ;
  NumericVector angle2;
  NumericVector pair;
  StringVector zone1   ;
  StringVector zone2    ;
};  

struct interaction_developed {
  int ant1;
  int ant2;
  double start;
  double end;
  int start_frame;
  int end_frame;
  double space;
  vector <string> type;
  vector <double> x1;
  vector <double> y1;
  vector <double> angle1;
  vector <string> zone1;
  vector <double> x2;
  vector <double> y2;
  vector <double> angle2;
  vector <string> zone2;
  bool interaction_finished;
};

struct summarised_interaction {
  int ant1;
  int ant2;
  double start;
  double end;
  double start_frame;
  double end_frame;
  double space;
  string type;
  double mean_x1;
  double mean_y1;
  double mean_angle1;
  double mean_x2;
  double mean_y2;
  double mean_angle2;
  string zone1;
  string zone2;
  int detections;
};

typedef vector <interaction_developed>  developed_interactions; 
typedef vector <summarised_interaction> summarised_interactions; 


double get_Average ( vector <double> v ){
  return(std::accumulate(v.begin(), v.end(), 0.0) / v.size());
};

double get_Circular_Mean ( vector <double> v ){
  double sum_sinus (0.0);double sum_cosinus (0.0);
  for (int i(0); i < v.size();i++){ //read angles
    sum_sinus   = sum_sinus   + sin(v[i]);
    sum_cosinus = sum_cosinus + cos(v[i]);
  }
  return(atan2 (sum_sinus,sum_cosinus));  
};

string simplify_type_list (vector <string> v){
  //initialize a new vector that will contain unique elements from StringVector v
  // where unique elements are separated by a comma within elements of v
  vector <string> unique_types;
  
  //also define a vector of string for each element of v
  vector <string> unique_types_per_element;
  
  // initalise variable deciding whether or not to add element
  bool add_it;
  // now loop over v
  for (int idx(0); idx < v.size(); idx++){
    //first clear unique_types_per_element
    unique_types_per_element.clear();
    //initialize first empty slot
    string current_string("");
    
    //then use ss to read the idxth element of v and separate it by commas
    std::stringstream ss(v[idx]);
    
    while (getline(ss, current_string, ',')) {
      //cout <<"current_string = " << current_string << endl;
      unique_types_per_element.push_back(current_string);
    }
    
    
   // for (string i; ss >> i;) {
     // current_string = current_string+i;    
    //  if (ss.peek() == ',')
    //    ss.ignore();
    //  cout << "current string = " << current_string << endl;
    //  unique_types_per_element.push_back(current_string);
    //  current_string = "";  
    //}
    
    //now loop over elements of unique_types_per_element, check if they are already listed in unique_types; and if not, add it
    
    for (int element(0);element<unique_types_per_element.size();element++){
      add_it = 1;
      if (unique_types.size()>0){
        for ( int elt(0); elt < unique_types.size(); elt++ ){
          if (  unique_types_per_element[element] == unique_types[elt]    ){
            add_it = 0;
          }
        }
      }
      if (add_it){
        unique_types.push_back(unique_types_per_element[element]);
      }
    }
  }
  //then let's sort the elements of unique types
  sort(unique_types.begin(), unique_types.end());
  
  //finally, define an output strung that contains all the types of unique_types separated by a comma
  string output (unique_types[0]);
  if (unique_types.size()>1){
    for (int elt(1); elt < unique_types.size(); elt++){
      output = output + "," + unique_types[elt];
    }
  }
  return(output);
};

interaction_developed define_new_interaction(int coll, collision_list all_collisions){
  interaction_developed new_interaction;
  new_interaction.ant1  =  all_collisions.ant1[coll];
  new_interaction.ant2  =   all_collisions.ant2[coll];
  new_interaction.start = all_collisions.time[coll];
  new_interaction.end   = all_collisions.time[coll];
  new_interaction.start_frame = all_collisions.frame[coll];
  new_interaction.end_frame   = all_collisions.frame[coll];
  new_interaction.space = all_collisions.space[coll];
  new_interaction.type.push_back(string(all_collisions.type[coll]));
  new_interaction.x1.push_back(all_collisions.x1[coll]);    
  new_interaction.y1.push_back(all_collisions.y1[coll]);    
  new_interaction.angle1.push_back(all_collisions.angle1[coll]);    
  new_interaction.zone1.push_back(string(all_collisions.zone1[coll]));    
  new_interaction.x2.push_back(all_collisions.x2[coll]);    
  new_interaction.y2.push_back(all_collisions.y2[coll]);    
  new_interaction.angle2.push_back(all_collisions.angle2[coll]);    
  new_interaction.zone2.push_back(string(all_collisions.zone2[coll]));
  new_interaction.interaction_finished = 0;
  return new_interaction;
};

void update_interaction(interaction_developed& current_interaction, int coll, collision_list all_collisions){
  current_interaction.end     =                all_collisions.time  [coll];
  current_interaction.end_frame     =                all_collisions.frame  [coll];
  current_interaction.type   .push_back(string(all_collisions.type  [coll]));
  current_interaction.x1     .push_back(       all_collisions.x1    [coll]);    
  current_interaction.y1     .push_back(       all_collisions.y1    [coll]);    
  current_interaction.angle1 .push_back(       all_collisions.angle1[coll]);    
  current_interaction.zone1  .push_back(string(all_collisions.zone1 [coll]));    
  current_interaction.x2     .push_back(       all_collisions.x2    [coll]);    
  current_interaction.y2     .push_back(       all_collisions.y2    [coll]);    
  current_interaction.angle2 .push_back(       all_collisions.angle2[coll]);    
  current_interaction.zone2  .push_back(string(all_collisions.zone2 [coll])); 
};

summarised_interaction close_interaction(interaction_developed& current_interaction){
  //declare the current interaction as closed
  current_interaction.interaction_finished = 1;
  
  //create an interaction_summary
  summarised_interaction interaction_summary;
  
  //fill in elements of interaction_summary
  interaction_summary.ant1         =                     current_interaction.ant1;
  interaction_summary.ant2         =                     current_interaction.ant2;
  interaction_summary.start        =                     current_interaction.start;
  interaction_summary.end          =                     current_interaction.end;
  interaction_summary.start_frame  =                     current_interaction.start_frame;
  interaction_summary.end_frame    =                     current_interaction.end_frame;
  
  interaction_summary.space        =                     current_interaction.space;
  interaction_summary.type         = simplify_type_list (current_interaction.type);
  interaction_summary.mean_x1      = get_Average        (current_interaction.x1);
  interaction_summary.mean_y1      = get_Average        (current_interaction.y1);
  interaction_summary.mean_angle1  = get_Circular_Mean  (current_interaction.angle1);
  interaction_summary.mean_x2      = get_Average        (current_interaction.x2);
  interaction_summary.mean_y2      = get_Average        (current_interaction.y2);
  interaction_summary.mean_angle2  = get_Circular_Mean  (current_interaction.angle2);
  interaction_summary.zone1        = simplify_type_list (current_interaction.zone1) ;
  interaction_summary.zone2        = simplify_type_list (current_interaction.zone2) ;
  interaction_summary.detections   = current_interaction.x1.size();
  
  return(interaction_summary);
};

// [[Rcpp::export]]
DataFrame  merge_interactions(DataFrame collisions, NumericVector pair_list, double max_distance_moved, double max_time_gap){
  /////////////////////////////////////////////////////////		
  //read variables from input collisions table
  collision_list          all_collisions;
  all_collisions.ant1   = collisions["ant1"];
  all_collisions.ant2   = collisions["ant2"];
  all_collisions.type   = collisions["types"];
  all_collisions.frame  = collisions["frame_number"];
  all_collisions.time   = collisions["time_second"];
  all_collisions.space  = collisions["space"];
  all_collisions.x1     = collisions["x1"];
  all_collisions.x2     = collisions["x2"];
  all_collisions.y1     = collisions["y1"];
  all_collisions.y2     = collisions["y2"];
  all_collisions.angle1 = collisions["angle1"];
  all_collisions.angle2 = collisions["angle2"];
  all_collisions.pair   = collisions["pair"];
  all_collisions.zone1  = collisions["zone1"];
  all_collisions.zone2  = collisions["zone2"];
  ///////////////////////////////////
  //get size of collisions table
  int collisions_size;
  collisions_size = all_collisions.ant1.size();
  
  ///////////////////////////////////
  //get number of pairs
  int pair_count;
  pair_count = pair_list.size();
  
  
  
  
  /////////////////////////////////////////////
  //initialise interactions per pair - developed
  vector <developed_interactions>  all_developed_interactions (pair_count);
  
  /////////////////////////////////////////////
  //initialise vector of all summarised interactions
  summarised_interactions   all_summarised_interactions ;
  
  
  ///////////////////////////////////
  ///Loop over collisions 
  interaction_developed last_interaction;
  double last_x1; double last_y1; double last_x2; double last_y2; double current_x1; double current_y1; double current_x2; double current_y2;
  double last_time; double current_time;
  int last_frame; int current_frame;
  int last_space; double current_space;
  double time_diff; double distance_diff1;double distance_diff2;
  for (int coll(0); coll < collisions_size;coll++){ //read all collisions
    ///first case: this is the first interaction between these two ants
    if (all_developed_interactions[all_collisions.pair[coll]].size() ==0 ){
      all_developed_interactions[all_collisions.pair[coll]].push_back(define_new_interaction(coll, all_collisions));
    }else{    ///second case: this is NOT the first interaction between these two ants
      ///in that case we want to look at the last interaction recorded between these two ants and evaluate if the current collisions is part of the same event
      last_interaction = all_developed_interactions[all_collisions.pair[coll]][all_developed_interactions[all_collisions.pair[coll]].size()-1];
      
      ///we have two criteria to define if two interactions are part of the same event: are the interactions temporally close and are the ants physically close to where they first were
      ///first, extract the last positions and last time of the last interaction
      last_x1   = last_interaction.x1[last_interaction.x1.size()-1];
      last_x2   = last_interaction.x2[last_interaction.x2.size()-1];
      last_y1   = last_interaction.y1[last_interaction.y1.size()-1];
      last_y2   = last_interaction.y2[last_interaction.y2.size()-1];
      last_time = last_interaction.end;
      last_frame = last_interaction.end_frame;
      last_space = last_interaction.space;
      
      ///second, extract the current positions and current time of collision
      current_x1   = all_collisions.x1[coll];
      current_x2   = all_collisions.x2[coll];
      current_y1   = all_collisions.y1[coll];
      current_y2   = all_collisions.y2[coll];
      current_time = all_collisions.time[coll];
      current_frame = all_collisions.frame[coll];
      current_space = all_collisions.space[coll];
      
      //use this to calculate distance diff
      distance_diff1 = sqrt(pow((current_x1-last_x1),2)+pow((current_y1-last_y1),2));
      distance_diff2 = sqrt(pow((current_x2-last_x2),2)+pow((current_y2-last_y2),2));
      time_diff      = current_time - last_time;
      
      //collision is part of the same interaction IF time_diff is < to time threshold AND both distance_diff are < to distance threshold AND collisions are in the same space
      if ( distance_diff1<=max_distance_moved && distance_diff2<=max_distance_moved && time_diff<=max_time_gap && current_space == last_space){
        //same interaction: update value of current interaction in all_developed_interactions; leave it open for now
        update_interaction(all_developed_interactions[all_collisions.pair[coll]][all_developed_interactions[all_collisions.pair[coll]].size()-1],coll,all_collisions);
        
      }else{
        //different interaction: close current interaction in all_developed_interactions; calculate summary values and push them to output variables
        all_summarised_interactions.push_back( close_interaction(all_developed_interactions[all_collisions.pair[coll]][all_developed_interactions[all_collisions.pair[coll]].size()-1]));
        // and create new interaction using current collision
        all_developed_interactions[all_collisions.pair[coll]].push_back(define_new_interaction(coll, all_collisions));
      }
    }
  }
  ///////////////////////////////////////////////
  ///finally, after reading all collisions, loop over all current interactions in all_developed_interactions, check if they have been officially closed, and if not, close current them; calculate summary values and push them to output variables
  for (int pair(0);pair<all_developed_interactions.size(); pair++){
    if (!all_developed_interactions[pair].empty()){
      if (all_developed_interactions[pair][all_developed_interactions[pair].size()-1].interaction_finished==0){
        all_summarised_interactions.push_back( close_interaction(all_developed_interactions[pair][all_developed_interactions[pair].size()-1]));
      }
    }
  }
  
  /////////////////////////////////////////////
  /// when this is done, prepare output
  
  /////////////////////////////////////////////
  //initialise interaction variables for Output
  vector <int>    ant1_interact;
  vector <int>    ant2_interact;
  vector <double> start_interact;
  vector <double> end_interact;
  vector <int>    start_interact_frame;
  vector <int>    end_interact_frame;
  vector <int>    space_interact;
  vector <string> type_interact;
  vector <double> meanx1_interact;
  vector <double> meany1_interact;
  vector <double> meanangle1_interact;
  vector <double> meanx2_interact;
  vector <double> meany2_interact;
  vector <double> meanangle2_interact;
  vector <string> zone1_interact;
  vector <string> zone2_interact;
  vector <int>    detections;
  /////////////////////////////////////////////
  ///loop over all_summarised_interactions
  for (int inter(0);inter < all_summarised_interactions.size();inter++){
    ant1_interact      .push_back(all_summarised_interactions[inter].ant1);
    ant2_interact      .push_back(all_summarised_interactions[inter].ant2);
    start_interact     .push_back(all_summarised_interactions[inter].start);
    end_interact       .push_back(all_summarised_interactions[inter].end);
    start_interact_frame     .push_back(all_summarised_interactions[inter].start_frame);
    end_interact_frame       .push_back(all_summarised_interactions[inter].end_frame);
    space_interact     .push_back(all_summarised_interactions[inter].space);
    type_interact      .push_back(all_summarised_interactions[inter].type);
    meanx1_interact    .push_back(all_summarised_interactions[inter].mean_x1);
    meany1_interact    .push_back(all_summarised_interactions[inter].mean_y1);
    meanangle1_interact.push_back(all_summarised_interactions[inter].mean_angle1);
    meanx2_interact    .push_back(all_summarised_interactions[inter].mean_x2);
    meany2_interact    .push_back(all_summarised_interactions[inter].mean_y2);
    meanangle2_interact.push_back(all_summarised_interactions[inter].mean_angle2);
    zone1_interact     .push_back(all_summarised_interactions[inter].zone1);
    zone2_interact     .push_back(all_summarised_interactions[inter].zone2);
    detections         .push_back(all_summarised_interactions[inter].detections);
  }
  
  
  
  //and export output variables
  return DataFrame::create(_["ant1"]             = ant1_interact,
                           _["ant2"]             = ant2_interact,
                           _["start"]            = start_interact,
                           _["end"]              = end_interact,
                           _["startframe"]       = start_interact_frame,
                           _["endframe"]         = end_interact_frame,
                           _["space"]            = space_interact,
                           _["types"]            = type_interact,
                           _["ant1.mean.x"]      = meanx1_interact,
                           _["ant1.mean.y"]      = meany1_interact,
                           _["ant1.mean.angle"]  = meanangle1_interact,
                           _["ant2.mean.x"]      = meanx2_interact,
                           _["ant2.mean.y"]      = meany2_interact,
                           _["ant2.mean.angle"]  = meanangle2_interact,
                           _["ant1.zones"]       = zone1_interact,
                           _["ant2.zones"]       = zone2_interact,
                           _["detections"]       = detections
  );
  
}


