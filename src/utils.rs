use parry2d_f64::math::Real;

pub const GRAVITATIONAL_ACCELERATION: Real = 9.81;
pub const COLLISION_EPS: Real = 1e-3;
pub const CONTACT_VELOCITY_EPS: Real = 2e-2;  // Valeur empirique, liée au dt