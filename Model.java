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
		// TODO this method implements one step of simulation with a step deltaT
		for (Ball b : balls) {
			// O(n²) not scalable, but works just fine for small amounts of balls 
			for (int i = 0; i < balls.length; i++) {
				if (b.collision(balls[i])) {
					System.out.println("Collision!");
					// handle collision
				}
				else{
					System.out.println("No collision");
				}
			}
			// detect collision with the border			
			if (b.x < b.radius || b.x > areaWidth - b.radius) {
				b.vx *= -1; // change direction of ball
			}
			if (b.y < b.radius || b.y > areaHeight - b.radius) {
				b.vy *= -1;
			}
			if ((b.x < b.radius || b.x > areaWidth - b.radius) && (b.y < b.radius || b.y > areaHeight - b.radius)) {
				b.vy *= -1;
				b.vx *= -1;
			}
			
			// compute new position according to the speed of the ball
			b.x  += deltaT * b.vx;
			b.y  += deltaT * b.vy; 
			b.vy += deltaT * -9.8;
		}
	}

	void handleCollision(Ball b1, Ball b2) { 
		collisionX = b1.x - b2.x;
		collisionY = b1.y - b2.y;
		Vector2D polarC = rectToPolar(new Vector2D(collisionX, collisionY));
		// using the formula below to translate the axis theta degrees
		// x' = x*cos(theta) + y*sin(theta)
		// y' = -x*sin(theta) + y*cos(theta)
		double theta = polarC.b;
		double newX = collisionX * Math.cos(theta) + collisionY * Math.sin(theta);
		double newY = -(collisionX * Math.sin(theta)) + collisionY * Math.cos(theta);

	

 
	}

	Vector2D rectToPolar(Vector2D v) {
		double length = Math.sqrt(v.a*v.a + v.b*v.b);
		double theta;
		switch (v) {
			case (v.a < 0 && v.b >= 0):
				theta = Math.arctan(v.b/v.a) + Math.PI;
				break;
			case (v.a < 0 && v.b < 0):
				theta = Math.arctan(v.b/v.a) - Math.PI;
				break;
			case (v.a == 0 && v.b > 0):
				theta = Math.PI/2;
				break;
			case (v.a == 0 && v.b < 0):
				theta = -(Math.PI/2);
				break;
			case (v.a == 0 && v.b == 0):
				throw Error("Undefined");
				break;
			default:
				theta = Math.arctan(v.b/v.a);
		}
		
		return Vector2D(length, theta); 
	}

	Vector2D polarToRect(Vector2D v) {
		double x = v.a * Math.cos(v.b);
		double y = v.a * Math.sin(v.b);
		return Vector2D(x, y); 
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

		boolean collision(Ball b) {
			if (this.euclideanDist(b) <= this.radius + b.radius){
				return true;
			}
			return false;
		}

		Vector2D collisionCoordinates(Ball b) {
			double x = Math.abs(this.x - b.x);
			double y = Math.abs(this.y - b.y);
			double angle = Math.atan2(y, x);
			double newX = Math.cos(angle);
			double newY = Math.sin(angle);
			return Vector2D(newX, newY);
		}

		double length() {
			return Math.sqrt(x*x + y*y);
		}

		double euclideanDist(Ball b) {
			// Calculates the euclidean distance to the ball b
			return Math.sqrt(Math.pow(x - b.x, 2) + Math.pow(y - b.y, 2));
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
			this.a = x;
			this.b = y;
		}
	}
}
