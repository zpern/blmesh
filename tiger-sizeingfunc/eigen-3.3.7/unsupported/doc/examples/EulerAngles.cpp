#include <spdlog/spdlog.h> 
 #include <unsupported/Eigen/EulerAngles>
#include <iostream>

using namespace Eigen;

int main()
{
  // A common Euler system by many armies around the world,
  //  where the first one is the azimuth(the angle from the north -
  //   the same angle that is show in compass)
  //  and the second one is elevation(the angle from the horizon)
  //  and the third one is roll(the angle between the horizontal body
  //   direction and the plane ground surface)
  // Keep remembering we're using radian angles here!
  typedef EulerSystem<-EULER_Z, EULER_Y, EULER_X> MyArmySystem;
  typedef EulerAngles<double, MyArmySystem> MyArmyAngles;
  
  MyArmyAngles vehicleAngles(
    3.14/*PI*/ / 2, /* heading to east, notice that this angle is counter-clockwise */
    -0.3, /* going down from a mountain */
    0.1); /* slightly rolled to the right */
  
  // Some Euler angles representation that our plane use.
  EulerAnglesZYZd planeAngles(0.78474, 0.5271, -0.513794);
  
  MyArmyAngles planeAnglesInMyArmyAngles = MyArmyAngles::FromRotation<true, false, false>(planeAngles);
  
  spdlog::info("vehicle angles(MyArmy):     " << vehicleAngles);
  spdlog::info("plane angles(ZYZ):        " << planeAngles);
  spdlog::info("plane angles(MyArmy):     " << planeAnglesInMyArmyAngles);
  
  // Now lets rotate the plane a little bit
  spdlog::info("==========================================================\n");
  spdlog::info("rotating plane now!\n");
  spdlog::info("==========================================================\n");
  
  Quaterniond planeRotated = AngleAxisd(-0.342, Vector3d::UnitY()) * planeAngles;
  
  planeAngles = planeRotated;
  planeAnglesInMyArmyAngles = MyArmyAngles::FromRotation<true, false, false>(planeRotated);
  
  spdlog::info("new plane angles(ZYZ):     " << planeAngles);
  spdlog::info("new plane angles(MyArmy): " << planeAnglesInMyArmyAngles);
  
  return 0;
}
