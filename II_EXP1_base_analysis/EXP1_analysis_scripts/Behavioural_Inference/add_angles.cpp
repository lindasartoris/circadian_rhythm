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
#include <random>
#include <chrono>

double get_vector_angle( double ori_x, double ori_y, double target_x, double target_y){
  double attractor_direction_angle;
  
  //0. if no movement, angle should be set to NA
  if ((ori_x==target_x)&(ori_y==target_y)){
    attractor_direction_angle = NA_REAL;
  }else{
    //1. use sin to determine if angle is angle is positive or negative
    double dx (target_x-ori_x); double dy (target_y-ori_y);
    /////////////////////////////////////

    double asinus;
    asinus = asin( dy/sqrt(pow(dx,2)+pow(dy,2)) );
    int sign_angle;
    if (asinus<0) {sign_angle = -1;}
    if (asinus==0){sign_angle = 0;}
    if (asinus>0) {sign_angle = 1;}
    
    //2. in case asin is equal to zero (i.e. point and target have the same y coordinate), use the x coordinates to determine the direction
    if (sign_angle==0){
      if (dx<0){
        attractor_direction_angle = M_PI;
      }else{
        attractor_direction_angle = 0;
      }
      
    }else{
      //3. in case asin is not equal to zero, use the cosinus to define the angle
      attractor_direction_angle = sign_angle*acos( dx/sqrt(pow(dx,2)+pow(dy,2)) );
    }
    
  }
  return(attractor_direction_angle );
};

double pi_frame(double angle){
  while(angle < (-M_PI)){angle = angle + 2*M_PI;}
  while(angle >   M_PI) {angle = angle - 2*M_PI;} 
  return(angle);
}
  
// [[Rcpp::export]]
DataFrame  add_angles(DataFrame Trajectory){
  // Fixed input parameters 
  //int max_gap_frame(4);
  /////////////////////////////////////////////////////////		
  //define variables from input Trajectory Table
  vector <double> x       = Trajectory["x"];
  vector <double> y       = Trajectory["y"];
  vector <double> Angle   = Trajectory["angle"];
  vector <int> Frame      = Trajectory["frame"];
  ///////////////////////////////////
  //get size of  table
  int table_size;
  table_size = x.size();
  
  /////////////////////////////////////
  //initialise angle variables
  vector <double> Movement_direction (table_size);
  vector <double> Turn_Angle (table_size);
  vector <double> Movement_minus_body_angle (table_size);
  vector <double> Movement_Body_Inclination_angle (table_size);
  vector <double> Body_angle_change (table_size);

  double movement_angle ;
  double turn_angle ;
  double movement_minus_body_angle ;
  double movement_body_inclination_angle;
  double body_angle_change ;
  
  
  ///////////////////////////////////
  ///Loop over trajectory fixes
  
  for (int traj_i(0); traj_i < table_size;traj_i++){ //read all interaction lines
    //Part A: calculate the local movement angle (absolute) and the difference between movement angle and body angle
    //not defined for the last line because we are missing the previous fix
    if (traj_i<=(table_size-2)){// Reminder: indices in C++ go from 0 to (length -1). Here the condition should thus be that (traj_i) is inferior or equal to (length-2)
          movement_angle                  = get_vector_angle(x[traj_i],y[traj_i],x[traj_i+1],y[traj_i+1]);
          movement_minus_body_angle       = movement_angle - Angle[traj_i];
          movement_body_inclination_angle = min(abs(pi_frame(movement_minus_body_angle)),M_PI-abs(pi_frame(movement_minus_body_angle)));
          body_angle_change               = Angle[traj_i+1] - Angle[traj_i];
    }else{
          movement_angle                  = NA_REAL;
          movement_minus_body_angle       = NA_REAL;
          movement_body_inclination_angle = NA_REAL;
          body_angle_change               = NA_REAL;
    }
    //Part B: calculate turn angle
    //not defined for the first line because we are missing the previous fix, or for the last line because we are missing the next fix
    if (traj_i>0 && traj_i<=(table_size-2)){// Reminder: indices in C++ go from 0 to (length -1). Here the condition should thus be that (traj_i) is inferior or equal to (length-2)
      turn_angle = movement_angle - get_vector_angle(x[traj_i-1],y[traj_i-1],x[traj_i],y[traj_i]);
    }else{
      turn_angle = NA_REAL;
    }
    
    Movement_direction              [traj_i] = pi_frame(movement_angle);
    Turn_Angle                      [traj_i] = pi_frame(turn_angle);
    Movement_minus_body_angle       [traj_i] = pi_frame(movement_minus_body_angle);
    Movement_Body_Inclination_angle [traj_i] = movement_body_inclination_angle;
    Body_angle_change               [traj_i] = pi_frame(body_angle_change);
    
    
  }//traj_i
  return DataFrame::create(_["Movement_direction"]               = Movement_direction,
                           _["TurnAngle"]                        = Turn_Angle,
                           _["Movement_Body_angle_diff"]         = Movement_minus_body_angle,
                           _["Movement_Body_Inclination_angle"]  = Movement_Body_Inclination_angle,
                           _["Body_Rotation"]                    = Body_angle_change
  );
}


