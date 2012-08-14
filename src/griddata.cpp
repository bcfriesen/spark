#include <griddata.h>

GridDataClass::GridDataClass() { }

GridDataClass::~GridDataClass() { }

void GridDataClass::set_vel(double gridvel) { vel = gridvel; }

void GridDataClass::set_rad(double gridrad) { rad = gridrad; }

double GridDataClass::get_rad() { return rad; }

double GridDataClass::get_vel() { return vel; }
