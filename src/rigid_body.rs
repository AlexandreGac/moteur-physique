use nalgebra::{Isometry2, Point2, Vector3};
use parry2d::shape::{Shape, ShapeType};
use parry2d::mass_properties::MassProperties;
use parry2d::math::Real;
use parry2d::shape::ConvexPolygon;

pub struct RigidBody {
    pub geometry: Box<dyn Shape>,
    pub transform: Isometry2<Real>,
    pub position: Vector3<Real>,
    pub velocity: Vector3<Real>,
    pub inv_mass: Vector3<Real>
}

impl RigidBody {
    pub fn new(shape: impl Shape, density: f32) -> Self {
        let mass_properties = shape.mass_properties(density);
        let (polygon, transform) = center_shape(&shape, mass_properties);
        RigidBody {
            geometry: Box::new(polygon),
            position: Vector3::new(0.0, 0.0, 0.0),
            velocity: Vector3::new(0.0, 0.0, 0.0),
            inv_mass: Vector3::new(mass_properties.inv_mass, mass_properties.inv_mass, mass_properties.inv_principal_inertia_sqrt),
            transform,
        }
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

    pub fn transform(&self) -> &Isometry2<f32> {
        &self.transform
    }
}

fn center_shape(shape: &dyn Shape, mass_properties: MassProperties) -> (impl Shape, Isometry2<f32>) {
    match shape.shape_type() {
        ShapeType::ConvexPolygon => {
            let convex_polygon = shape.as_convex_polygon().unwrap();
            let center = mass_properties.world_com(&Isometry2::identity()) - Point2::origin();
            let points = convex_polygon.points().iter().map(|p| p - center).collect();
            let polygon = ConvexPolygon::from_convex_polyline(points).expect("Polygone non convexe !");
            let transform = Isometry2::translation(center.x, center.y);
            (polygon, transform)
        },
        _ => panic!("Géométrie non gérée !")
    }
}

pub struct RigidBodyState {
    pub position: Vector3<f32>,
    pub velocity: Vector3<f32>
}

impl RigidBodyState {

    pub fn apply_physics_step(&self, velocity: Vector3<Real>, acceleration: Vector3<Real>, dt: Real) -> RigidBodyState {
        RigidBodyState {
            position: self.position + dt * velocity,
            velocity: self.velocity + dt * acceleration
        }
    }
}