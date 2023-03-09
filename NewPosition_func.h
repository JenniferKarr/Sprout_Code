#ifndef NewPosition_func_H_   /* Include guard */
#define NewPosition_func_H_

int NewPosition(double new_xyz[3], double xyz[3], double vel[3]);

int NewPosition4givenTau(double tau_goal, double new_xyz[3], double xyz[3], double vel[3]);
#endif // newPosition_func_H_
