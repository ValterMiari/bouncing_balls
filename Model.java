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
					//Ball[] newBalls = handleCollision(b, balls[i]);
					Vector2D newVelocity = handleCollision(b, balls[i]);
					b.vx = newVelocity.a;
					balls[i].vx = newVelocity.b;
					break;
					// handle collision
					/*vx1 = b.x + balls[i].mass*balls[i].x;
					vy1 = b.y + balls[i].mass*balls[i].y;

					vx2 = b.x - balls[i].x;
					vy2 = b.y - balls[i].y;

					b.vx = vx1;
					b.vy = vy1;

					balls[i].vx = vx2;
					balls[i].vy = vy2;*/

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

	Vector2D handleCollision(Ball b1, Ball b2) { 
		double dx, dy, theta, pB1x, pB1y, pB2x, pB2y, rB1vx, rB1vy, rB2vx, rB2vy, nB1vx, nB2vx;
		Vector2D pB1, pB2;
		Ball transformedB1, transformedB2;

		// used to be abs
		dx = (b1.x - b2.x);
		dy = (b1.y - b2.y);
		theta = Math.atan2(dy, dx);

		pB1 = rectToPolar(new Vector2D(b1.x, b1.y));
		pB2 = rectToPolar(new Vector2D(b2.x, b2.y));
		// using the formula below to rotate the axis theta degrees
		// x' = r*cos(alpha)*cos(theta) + r*sin(alpha)*sin(theta)
		// y' = rsin(alpha)cos(theta) - rcos(alpha)sin(theta)
		pB1x = pB1.a * Math.cos(pB1.b) * Math.cos(theta) + pB1.a * Math.sin(pB1.b) * Math.sin(theta);
		pB1y = pB1.a * Math.sin(pB1.b) * Math.cos(theta) - pB1.a * Math.cos(pB1.b) * Math.sin(theta);
		pB2x = pB2.a * Math.cos(pB2.b) * Math.cos(theta) + pB2.a * Math.sin(pB2.b) * Math.sin(theta);
		pB2y = pB2.a * Math.sin(pB2.b) * Math.cos(theta) - pB2.a * Math.cos(pB2.b) * Math.sin(theta);

		//rotating velocities with the formula 
		// x' = x*cos(theta) + y*sin(theta)
		// y' = -x*sin(theta) + y*cos(theta)
		rB1vx = b1.vx * Math.cos(-theta) + b1.vy * Math.sin(-theta);
		rB1vy = -(b1.vx * Math.sin(-theta)) + b1.vy * Math.cos(-theta);
		
		rB2vx = b2.vx * Math.cos(-theta) + b2.vy * Math.sin(-theta);
		rB2vy = -(b2.vx * Math.sin(-theta)) + b2.vy * Math.cos(-theta);	
		

		transformedB1 = new Ball(pB1x, pB1y, rB1vx, rB1vy, b1.radius);//new Vector2D(pB1x, pB1y);
		transformedB2 = new Ball(pB2x, pB2y, rB2vx, rB2vy, b2.radius);

		Vector2D impactOnCollision = collisionImpact(transformedB1, transformedB2); 
		
		// final speed after rotation and impact, rotate the speed back to original system
		nB1vx = impactOnCollision.a * Math.cos(-theta) - rB1vy * Math.sin(-theta);
		nB2vx = impactOnCollision.b * Math.cos(-theta) - rB2vy * Math.sin(-theta); 

		Ball newB1 = new Ball(b1.x, b1.y, nB1vx, b1.vy, b1.radius); 
		Ball newB2 = new Ball(b2.x, b2.y, nB2vx, b2.vy, b2.radius);
		Ball[] newBalls = {newB1, newB2};

		//return newBalls;
		return new Vector2D(nB1vx, nB2vx);
		}

	Vector2D collisionImpact(Ball b1, Ball b2) {
		// note: I believe that only the velocity in the x directrion is the one that matters, since it is transformed to
		// a horizontal collision
		double v1 = ((b1.mass - b2.mass) / (b1.mass + b2.mass)) * b1.vx 
					+ ((2 * b2.mass) / (b1.mass + b2.mass)) * b2.vx; 
		double v2 = ((2 * b1.mass) / (b1.mass + b2.mass)) * b1.vx
					+ ((b2.mass - b1.mass)/(b1.mass + b2.mass)) * b2.vx;

		return new Vector2D(v1, v2);
	}

	Vector2D rectToPolar(Vector2D v) {
		double length = Math.sqrt(v.a*v.a + v.b*v.b);
		double theta = Math.atan2(v.b, v.a);

		return new Vector2D(length, theta); 
	}

	Vector2D polarToRect(Vector2D v) {
		double x = v.a * Math.cos(v.b);
		double y = v.a * Math.sin(v.b);
		return new Vector2D(x, y); 
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

		double length() {
			return Math.sqrt(x*x + y*y);
		}

		double euclideanDist(Ball b) {
			// Calculates the euclidean distance to the ball b
			return Math.sqrt(Math.pow(x - b.x + 1/100, 2) + Math.pow(y - b.y + 1/100, 2));
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
