//! @file tsc_geometry.geo Definition of geometry by inclusion of geofuns.
//!

// OUTER BOUNDARY
Include "Rectangle.geofun"; // required for outer boundary

// OBSTACLE GEOMETRY
//Include "ObstaclesSPJF2004.geofun";  
//Include "ObstaclesThreeEllipses.geofun";  
//Include "ObstaclesOneEllipse.geofun"; // this has a bothersome control point at (.5,.5) 
Include "ObstaclesOneEllipse2.geofun";  
//Include "ObstaclesTwoRectangles.geofun";  

// The grid has to be scaled via Grid('*.geo', hmax)!





