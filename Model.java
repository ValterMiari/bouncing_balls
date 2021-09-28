package bouncing_balls;
import java.lang.Math;
import java.util.Vector;

/**
 * The physics model.
 * 
 * This class is where you should implement your bouncing balls model.
 * 
 * The code has intentionally been kept as simple as possible, but if you wish, you can improve the design.
 * 
 * @author Simon Robillard
 *
 */
class Model {

	double areaWidth, areaHeight;
	
	Ball [] balls;

	Model(double width, double height) {
		areaWidth = width;
		areaHeight = height;
		
		// Initialize the model with a few balls
		balls = new Ball[2];
		balls[0] = new Ball(width / 3, height * 0.9, 1.2, 1.6, 0.2);
		balls[1] = new Ball(2 * width / 3, height * 0.7, -0.6, 0.6, 0.3);
	}

	void step(double deltaT) {
		// this method implements one step of simulation with a step deltaT
		for (Ball b : balls) {
			// O(n²), not scalable, but works just fine for small amounts of balls 
			for (int i = 0; i < balls.length; i++) {
				if (b.collision(balls[i])) {	
					Vector2D[] newVelocity = handleCollision(b, balls[i], deltaT);
					b.vx = newVelocity[0].a;
					b.vy = newVelocity[0].b;
					balls[i].vx = newVelocity[1].a;
					balls[i].vy = newVelocity[1].b; 
				}
			}
			// Detect collision with the border			
			// Simulating 0.02 seconds ahead, such that a simulation with 60 fps avoids not registering a collision with a wall
			if (b.x + 0.02 * b.vx < b.radius || b.x + 0.02 * b.vx > areaWidth - b.radius) {
				b.vx *= -1; // change direction of ball
			}
			if (b.y + 0.02 * b.vy < b.radius || b.y + 0.02 * b.vy > areaHeight - b.radius) {
				b.vy *= -1;
			}
			if ((b.x + 0.02 * b.vx < b.radius || b.x + 0.02 * b.vx > areaWidth - b.radius) && (b.y + 0.02 * b.vy < b.radius || b.y + 0.02 * b.vy > areaHeight - b.radius)) {
				b.vy *= -1;
				b.vx *= -1;
			}
			
			// compute new position according to the speed of the ball
			b.x  += deltaT * b.vx;
			b.y  += deltaT * b.vy; 
			b.vy += deltaT * -9.8;
		}
	}

	Vector2D[] handleCollision(Ball b1, Ball b2, double deltaT) { 

		double dx, dy, theta, rB1vx, rB1vy, rB2vx, rB2vy, nB1vx, nB2vx, nB1vy, nB2vy;
		Vector2D impactOnCollision, b1Velocity, b2Velocity;

		// Distance between the both balls and their impact angle
		dx = (b1.x - b2.x);
		dy = (b1.y - b2.y);
		theta = Math.atan2(dy, dx);

		// Rotating the velocities with the impact angle theta
		rB1vx = b1.vx * Math.cos(theta) + b1.vy * Math.sin(theta);
		rB1vy = b1.vy * Math.cos(theta) - (b1.vx * Math.sin(theta));
		
		rB2vx = b2.vx * Math.cos(theta) + b2.vy * Math.sin(theta);
		rB2vy = b2.vy * Math.cos(theta) - (b2.vx * Math.sin(theta));	

		// Calculating the one dimensional collision impact, made possible by the rotation
		impactOnCollision = collisionImpact(rB1vx, b1.mass, rB2vx, b2.mass); 
		
		// The final velocity after rotation and impact has now been calculated
		// Now rotate the velocity back to original coordinate system
		nB1vx = impactOnCollision.a * Math.cos(theta) - rB1vy * Math.sin(theta);
		nB1vy = impactOnCollision.a * Math.sin(theta) + rB1vy * Math.cos(theta);

		nB2vx = impactOnCollision.b * Math.cos(theta) - rB2vy * Math.sin(theta); 
		nB2vy = impactOnCollision.b * Math.sin(theta) + rB2vy * Math.cos(theta);

 		b1Velocity = new Vector2D(nB1vx, nB1vy);
		b2Velocity = new Vector2D(nB2vx, nB2vy);

		return new Vector2D[] {b1Velocity, b2Velocity}; 
	}

	// Calculates the impact of the collision on the given balls in one dimension
	Vector2D collisionImpact(double u1, double m1, double u2, double m2) {
		double v1, v2;

		v1 = ((m1 - m2) / (m1 + m2)) * u1 
			 + ((2 * m2) / (m1 + m2)) * u2; 
		v2 = ((2 * m1) / (m1 + m2)) * u1 
			 + ((m2 - m1)/(m1 + m2)) * u2;

		return new Vector2D(v1, v2);
	}

	/**
	 * Simple inner class describing balls.
	 */
	class Ball {
		
		Ball(double x, double y, double vx, double vy, double r) {
			this.x = x;
			this.y = y;
			this.vx = vx;
			this.vy = vy;
			this.radius = r;
			// We assume that the desity of the ball is 1
			this.mass = Math.pow(r, 3);
		}

		// Checks if a ball collides with another ball
		boolean collision(Ball b) {
			if (this.x == b.x && this.y == b.y){
				return false;
			}
			if (this.euclideanDist(b) <= this.radius + b.radius){
				return true;
			}
			return false;
		}

		double length() {
			return Math.sqrt(x*x + y*y);
		}

		double euclideanDist(Ball b) {
			// Calculates the euclidean distance to the ball b
			// Simulating 0.02 seconds ahead, so that collisions get handled correctly in high speeds with 60 fps
			return Math.sqrt(Math.pow((x + 0.02 * vx) - b.x, 2) + Math.pow((y + 0.02 * vy) - b.y, 2));
		}

		/**
		 * Position, speed, and radius of the ball. You may wish to add other attributes.
		 */
		double x, y, vx, vy, radius, mass;
	}
	/**
	 * A two dimensional vector.
	 */
	class Vector2D {
		
		double a, b;
		
		Vector2D(double a, double b) {
			this.a = a;
			this.b = b;
		}
	}
}
