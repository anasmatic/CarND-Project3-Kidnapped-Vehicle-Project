# Particle Filter Project 
Self-Driving Car Engineer Nanodegree Program

In this project , robot is kidnapped and has no idea where is it.
we have to detect its accurate possition , by first reciving gps coordinates , then reinforce it by dedecting map landmarks.

#### default project installation :
1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./particle_filter

#### my own Build Instructions
I was using Windows 10 and VisualStudio17

-to build this project using **Bash for window** :

    navigate to projet
    write cmd : mkdir build
    then navigate to build
    write cmd : cmake .. -G "Unix Makefiles" && make


## particle filter Steps

1. init function, recive first gps coordinates, create initial particles and add randome gaussian noise to each.
2. prediect robot position, using equations for x , y , theta
3. after predection we need to update particles weights according to map landmarks readings
   * we need to convert landmarks from car coordinates to map coordinates
   * then detect the nearst land marks
   * finaly update partical weight regarding how close it is to the map-landmark
4. final step is to resample the particles :
   * in this project we were adviced to use discrete_distribution
   * I was trying to use Resampling Wheel, but untill now I can't get good error value using it.

## trial and error

using normalization after prediction step give high error results.

number of particales :
if less than 50 , it gives highter error and takes more time to compute 
: `@ 25 Error x:0.137 y:0.129 Time:75s`

if more tha 100 gives better error, but takes long time
: `@ 200 Error x:0.110 y:0.103 Time:70s`

so , I use number of particales = 75 which gives these results :
: `@ 75 Error x:0.113 y:0.108 Time:53s`
