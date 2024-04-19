use nalgebra::{Isometry2, Point2, Rotation2, Vector2, Vector3};
use parry2d_f64::shape::{Shape, ShapeType};
use parry2d_f64::mass_properties::MassProperties;
use parry2d_f64::math::Real;
use parry2d_f64::shape::ConvexPolygon;
use piston_window::{Context, ellipse, Graphics, line, polygon};
use piston_window::math::Scalar;
use crate::utils::GRAVITATIONAL_ACCELERATION;

pub struct RigidBody {
    pub geometry: Box<dyn Shape>,
    pub transform: Isometry2<Real>,
    pub position: Vector3<Real>,
    pub velocity: Vector3<Real>,
    pub inv_mass: Vector3<Real>
}

impl RigidBody {
    pub fn new(shape: impl Shape, density: Real) -> Self {
        let mass_properties = shape.mass_properties(density);
        let (new_shape, transform) = center_shape(&shape, mass_properties);
        RigidBody {
            geometry: new_shape,
            position: Vector3::new(transform.translation.x, transform.translation.y, 0.0),
            velocity: Vector3::new(0.0, 0.0, 0.0),
            inv_mass: Vector3::new(mass_properties.inv_mass, mass_properties.inv_mass, mass_properties.inv_principal_inertia_sqrt),
            transform,
        }
    }

    pub fn set_transform(&mut self, transform: Isometry2<Real>) {
        self.position = Vector3::new(transform.translation.x, transform.translation.y, transform.rotation.angle());
        self.transform = transform;
    }

    pub fn get_state(&self) -> RigidBodyState {
        RigidBodyState {
            position: self.position,
            velocity: self.velocity
        }
    }

    pub fn apply_impulse(&mut self, impulse: Vector3<Real>) {
        self.velocity += impulse.component_mul(&self.inv_mass);
    }

    pub fn apply_physics(&mut self, velocity: Vector3<Real>, acceleration: Vector3<Real>, dt: Real) {
        self.position += dt * velocity;
        self.velocity += dt * acceleration;
        self.transform = Isometry2::new(self.position.xy(), self.position.z);
    }

    pub fn shape(&self) -> &dyn Shape {
        self.geometry.as_ref()
    }

    pub fn transform(&self) -> &Isometry2<Real> {
        &self.transform
    }

    pub fn get_energy(&self) -> Real {
        let kinetic_energy = Vector3::new(0.5, 0.5, 0.5)
            .component_div(&self.inv_mass)
            .component_mul(&self.velocity)
            .component_mul(&self.velocity)
            .sum();
        let potential_energy = GRAVITATIONAL_ACCELERATION / self.inv_mass.x * self.position.y;

        kinetic_energy + potential_energy
    }

    pub fn display(&self, context: Context, graphics: &mut impl Graphics, width: Scalar, height: Scalar) {
        match self.geometry.shape_type() {
            ShapeType::ConvexPolygon => {
                polygon(
                    [1.0, 0.0, 0.0, 0.5],
                    self.geometry
                        .as_convex_polygon()
                        .expect("Erreur forme !")
                        .points()
                        .iter()
                        .map(|v| self.transform * v)
                        .map(|v| [width / 2.0 + (v.x * 100.0) as Scalar, height - (v.y * 100.0) as Scalar])
                        .collect::<Vec<_>>()
                        .as_slice(),
                    context.transform,
                    graphics
                );
            }
            ShapeType::Ball => {
                let ball = self.geometry.as_ball().expect("Erreur forme !");
                ellipse(
                    [1.0, 0.0, 0.0, 0.5],
                    [
                        width / 2.0 + ((self.position.x - ball.radius) * 100.0) as Scalar,
                        height - ((self.position.y + ball.radius) * 100.0) as Scalar,
                        (200.0 * ball.radius) as Scalar,
                        (200.0 * ball.radius) as Scalar
                    ],
                    context.transform,
                    graphics
                );
                let rotation = Rotation2::new(self.position.z);
                let point_1: Vector2<Real> = self.position.xy();
                let point_2: Vector2<Real> = self.position.xy() + rotation * Vector2::<Real>::new(ball.radius, 0.0);
                line([1.0, 1.0, 1.0, 1.0],
                     1.0,
                     [
                         width / 2.0 + (point_1.x * 100.0) as Scalar,
                         height - (point_1.y * 100.0) as Scalar,
                         width / 2.0 + (point_2.x * 100.0) as Scalar,
                         height - (point_2.y * 100.0) as Scalar,
                     ],
                     context.transform,
                     graphics
                );
            }
            _ => panic!("Géométrie non gérée !")
        }
    }
}

fn center_shape(shape: &dyn Shape, mass_properties: MassProperties) -> (Box<dyn Shape>, Isometry2<Real>) {
    match shape.shape_type() {
        ShapeType::ConvexPolygon => {
            let convex_polygon = shape.as_convex_polygon().unwrap();
            let center = mass_properties.world_com(&Isometry2::identity()) - Point2::origin();
            let points = convex_polygon.points().iter().map(|p| p - center).collect();
            let polygon = ConvexPolygon::from_convex_polyline(points).expect("Polygone non convexe !");
            let transform = Isometry2::translation(center.x, center.y);
            (Box::new(polygon), transform)
        },
        ShapeType::Ball => {
            (shape.clone_box(), Isometry2::identity())
        },
        _ => panic!("Géométrie non gérée !")
    }
}

#[derive(Debug)]
pub struct RigidBodyState {
    pub position: Vector3<Real>,
    pub velocity: Vector3<Real>
}

impl RigidBodyState {

    pub fn apply_physics_step(&self, velocity: Vector3<Real>, acceleration: Vector3<Real>, dt: Real) -> RigidBodyState {
        RigidBodyState {
            position: self.position + dt * velocity,
            velocity: self.velocity + dt * acceleration
        }
    }
}