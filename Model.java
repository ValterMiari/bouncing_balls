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
					//handleCol(b, balls[i], deltaT);
					System.out.println("Collision!");
				}
				else{
					System.out.println("No collision");
				}
			}
			// detect collision with the border			
			// foreshodowing with 0.01 seconds, such that a frame with 60 fps avoids not registering a collision
			if (b.x + 0.03 * b.vx < b.radius || b.x + 0.02 * b.vx > areaWidth - b.radius) {
				b.vx *= -1; // change direction of ball
			}
			if (b.y + 0.03 * b.vy < b.radius || b.y + 0.03 * b.vy > areaHeight - b.radius) {
				b.vy *= -1;
			}
			if ((b.x + 0.03 * b.vx < b.radius || b.x + 0.03 * b.vx > areaWidth - b.radius) && (b.y + 0.03 * b.vy < b.radius || b.y + 0.03 * b.vy > areaHeight - b.radius)) {
				b.vy *= -1;
				b.vx *= -1;
			}
			
			// compute new position according to the speed of the ball
			b.x  += deltaT * b.vx;
			b.y  += deltaT * b.vy; 
			b.vy += deltaT * -9.8;
		}
	}

	// make this function return void, and insted modify the values of the balls in this function
	Vector2D[] handleCollision(Ball b1, Ball b2, double deltaT) { 
		double dx, dy, theta, rB1x, rB1y, rB2x, rB2y, rB1vx, rB1vy, rB2vx, rB2vy, nB1vx, nB2vx, nB1vy, nB2vy;
		Vector2D pB1, pB2;
		Ball transformedB1, transformedB2;

		// used to be abs, and b1.x - b2.x
		dx = (b1.x - b2.x);
		dy = (b1.y - b2.y);
		theta = Math.atan2(dy, dx);

		/* pB1 = rectToPolar(new Vector2D(b1.x, b1.y));
		pB2 = rectToPolar(new Vector2D(b2.x, b2.y));
		
		// using the formula below to rotate the axis theta degrees
		// x' = r*cos(alpha)*cos(theta) + r*sin(alpha)*sin(theta)
		// y' = rsin(alpha)cos(theta) - rcos(alpha)sin(theta)
		rB1x = pB1.a * Math.cos(pB1.b) * Math.cos(theta) + pB1.a * Math.sin(pB1.b) * Math.sin(theta);
		rB1y = pB1.a * Math.sin(pB1.b) * Math.cos(theta) - pB1.a * Math.cos(pB1.b) * Math.sin(theta);
		rB2x = pB2.a * Math.cos(pB2.b) * Math.cos(theta) + pB2.a * Math.sin(pB2.b) * Math.sin(theta);
		rB2y = pB2.a * Math.sin(pB2.b) * Math.cos(theta) - pB2.a * Math.cos(pB2.b) * Math.sin(theta);
 */
		/* double x1 = 0, y1 = 0;
		double x2 = dx * Math.cos(theta) + dy * Math.sin(theta);
		double y2 = dy * Math.cos(theta) - dx * Math.sin(theta); */
		//double x2 = b2.x * Math.cos(theta) + b2.y * Math.sin(theta);
		//double y2 = b2.y * Math.cos(theta) - (b2.x * Math.sin(theta)); 
		//rotating velocities with the formula 
		// x' = x*cos(theta) + y*sin(theta)
		// y' = -x*sin(theta) + y*cos(theta)
		rB1vx = b1.vx * Math.cos(theta) + b1.vy * Math.sin(theta);
		rB1vy = b1.vy * Math.cos(theta) - (b1.vx * Math.sin(theta));
		
		rB2vx = b2.vx * Math.cos(theta) + b2.vy * Math.sin(theta);
		rB2vy = b2.vy * Math.cos(theta) - (b2.vx * Math.sin(theta));	
		
		// The balls rotated to the new axises
		transformedB1 = new Ball(b1.x, b1.y, rB1vx, rB1vy, b1.radius);
		transformedB2 = new Ball(b2.x, b2.y, rB2vx, rB2vy, b2.radius);

		Vector2D impactOnCollision = collisionImpact(transformedB1, transformedB2); 
		
		/* double totalSpeed = Math.abs(impactOnCollision.a + impactOnCollision.b);
		double overlap = (b1.radius + b2.radius) - Math.abs(b1.x- b2.x);//Math.abs(rB1x - rB2x);
		b1.x += (impactOnCollision.a / totalSpeed * overlap) * deltaT;
					//b.y += (overlap / 2) * b.vy * deltaT;
		b2.x += (impactOnCollision.b / totalSpeed * overlap) * deltaT;   */
					//balls[i].x += (overlap / 2) * balls[i].vy * deltaT; 

		/* double overlap = (b1.radius + b2.radius) - Math.abs(b1.x - b2.x);//Math.abs(rB1x - rB2x);	
		b1.x += overlap/2 * impactOnCollision.a * deltaT;
		b2.x += overlap/2 * impactOnCollision.b * deltaT; */
 
		/* // relative pos, rotate back
		double x1final = x1 * Math.cos(-theta) + y1 * Math.sin(-theta);
		double y1final = y1 * Math.cos(-theta) - x1 * Math.sin(-theta);
		double x2final = x2 * Math.cos(-theta) + y2 * Math.sin(-theta);
		double y2final = y2 * Math.cos(-theta) - x2 * Math.sin(-theta);

		// new abs position
		b1.x = b1.x + x1final;
		b1.y = b1.y + y1final;

		b2.x = b1.x + x2final;
		b2.y = b1.y + y2final; */

		// final speed after rotation and impact, rotate the speed back to original system
		nB1vx = impactOnCollision.a * Math.cos(-theta) - rB1vy * Math.sin(-theta);
		nB1vy = -impactOnCollision.a * Math.sin(-theta) + rB1vy * Math.cos(-theta);

		nB2vx = impactOnCollision.b * Math.cos(-theta) - rB2vy * Math.sin(-theta); 
		nB2vy = -impactOnCollision.b * Math.sin(-theta) + rB2vy * Math.cos(-theta);

		Vector2D b1Speed = new Vector2D(nB1vx, nB1vy);
		Vector2D b2Speed = new Vector2D(nB2vx, nB2vy);
		Vector2D[] result = {b1Speed, b2Speed};

		return result;
	}

	void handleCol(Ball b1, Ball b2, double deltaT) {
		double dx = b1.x - b2.x;
		double dy = b1.y - b2.y;
		double theta = Math.atan2(dy, dx);

		double dvx = b1.vx - b2.vx;
		double dvy = b1.vy - b2.vy;
		
		// avoid accidental overlapping
		//if (dvx * dx + dvy * dy <= 0) {
			b1.rotateVelocity(theta);
			b2.rotateVelocity(theta);

			Vector2D impactOnCollision = collisionImpact(b1, b2);
			b1.vx = impactOnCollision.a;
			b2.vx = impactOnCollision.b;
			/* b1.vy = b1.vy;
			b2.vy = b2.vy; */

			double x1 = 0, y1 = 0;
			double x2 = dx * Math.cos(theta) + dy * Math.sin(theta);
			double y2 = dy * Math.cos(theta) - dx * Math.sin(theta);

			double totalSpeed = Math.abs(impactOnCollision.a - impactOnCollision.b);
			double overlap = (b1.radius + b2.radius) - (b1.x - b2.x);//Math.abs(rB1x - rB2x);
			/* b1.x += impactOnCollision.a / totalSpeed * overlap * deltaT; 
			b2.x += impactOnCollision.b / totalSpeed * overlap * deltaT; */
			b1.x += (impactOnCollision.a * overlap * totalSpeed / 2) * deltaT;
			b2.x += (impactOnCollision.b * overlap * totalSpeed / 2) * deltaT;
			


			/* // relative pos, rotate back
		double x1final = x1 * Math.cos(-theta) + y1 * Math.sin(-theta);
		double y1final = y1 * Math.cos(-theta) - x1 * Math.sin(-theta);
		double x2final = x2 * Math.cos(-theta) + y2 * Math.sin(-theta);
		double y2final = y2 * Math.cos(-theta) - x2 * Math.sin(-theta);

		// new abs position
		b1.x = b1.x + x1final * deltaT;
		b1.y = b1.y + y1final * deltaT;

		b2.x = b1.x + x2final * deltaT;
		b2.y = b1.y + y2final * deltaT;
 */

			// rotate velocity back
			b1.rotateVelocity(-theta);
			b2.rotateVelocity(-theta);


			
		//}

		/*double x1 = 0, y1 = 0;
    	double x2 = dx * Math.cos(theta) + dy * Math.sin(theta);
    	double y2 = dy * Math.cos(theta) - dx * Math.sin(theta);*/
		
		/* b1.rotateVelocity(theta);
		b2.rotateVelocity(theta); 

		
		// handle overlap
		double totalSpeed = Math.abs(impactOnCollision.a) + Math.abs(impactOnCollision.b);
		double overlap = (b1.radius + b2.radius) - Math.abs(x1 - x2);//Math.abs(rB1x - rB2x);
		x1 += impactOnCollision.a / totalSpeed * overlap * deltaT; 
		x2 += impactOnCollision.b / totalSpeed * overlap * deltaT;  
		
		double x1final = x1 * Math.cos(-theta) + y1 * Math.sin(-theta); 
		double y1final = y1 * Math.cos(-theta) - x1 * Math.sin(theta); 
		double x2final = x2 * Math.cos(-theta) + y2 * Math.sin(-theta); 
		double y2final = y2 * Math.cos(-theta) - x2 * Math.sin(-theta); 

		// new abs position
		b2.x = b1.x + x2final; 
		b2.y = b1.y + y2final;
		b1.x = b1.x + x1final;
		b1.y = b1.y + y1final;


		// the collision is now handled, rotate the speed and position back to original axis
		b1.rotateVelocity(-theta);
		b2.rotateVelocity(-theta); */
	}

	Vector2D collisionImpact(Ball b1, Ball b2) {
		double m1, m2, u1, u2, v1, v2;

		m1 = b1.mass;
		m2 = b2.mass;
		u1 = b1.vx;
		u2 = b2.vx;

		v1 = ((m1 - m2) / (m1 + m2)) * u1 
			 + ((2 * m2) / (m1 + m2)) * u2; 
		v2 = ((2 * m1) / (m1 + m2)) * u1 
			 + ((m2 - m1)/(m1 + m2)) * u2;

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
			// Foreshadowing with 0.01 seconds, so that collisions get handled correctly in high speeds
			return Math.sqrt(Math.pow((x + 0.02 * vx) - b.x, 2) + Math.pow((y + 0.02 * vy) - b.y, 2));
		}

		void rotateVelocity(double angle) {
			/* x = x * Math.cos(angle) + y * Math.sin(angle);
			y = y * Math.cos(angle) - x * Math.sin(angle);
 */
			vx = vx * Math.cos(angle) + vy * Math.sin(angle);
			vy = vy * Math.cos(angle) - vx * Math.sin(angle);
		/* 	rB1vx = b1.vx * Math.cos(theta) + b1.vy * Math.sin(theta);
			rB1vy = b1.vy * Math.cos(theta) - (b1.vx * Math.sin(theta));
		
			rB2vx = b2.vx * Math.cos(theta) + b2.vy * Math.sin(theta);
			rB2vy = b2.vy * Math.cos(theta) - (b2.vx * Math.sin(theta));	
 */		
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
