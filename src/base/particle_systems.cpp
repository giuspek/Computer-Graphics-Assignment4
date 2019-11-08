#include "particle_systems.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

using namespace std;
using namespace FW;

boolean winder = false;

namespace {

	inline Vec3f fGravity(float mass) {
		return Vec3f(0, -9.8f * mass, 0);
	}

	// force acting on particle at pos1 due to spring attached to pos2 at the other end
	inline Vec3f fSpring(const Vec3f& pos1, const Vec3f& pos2, float k, float rest_length) {
		// YOUR CODE HERE (R2)
		Vec3f res = k * (rest_length - (pos1 - pos2).length()) * ((pos1 - pos2) / (pos1 - pos2).length());
		return res;
	}

	inline Vec3f fDrag(const Vec3f& v, float k) {
		// YOUR CODE HERE (R2)
		Vec3f res = -k * v;
		return res;
	}

	inline unsigned vel(unsigned index) {
		return (2 * index) + 1;
	}

	inline unsigned pos(unsigned index) {
		return 2 * index;
	}

} // namespace

void SimpleSystem::reset() {
	state_ = State(1, Vec3f(0, radius_, 0));
}

State SimpleSystem::evalF(const State& state) const {
	State f(1, Vec3f(-state[0].y, state[0].x, 0));
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H
// using the implicit Euler method, the simple system should converge towards origin -- as opposed to the explicit Euler, which diverges outwards from the origin.
void SimpleSystem::evalJ(const State&, SparseMatrix& result, bool initial) const {
	if (initial) {
		result.coeffRef(1, 0) = 1.0f;
		result.coeffRef(0, 1) = -1.0f;
	}
}
#endif

Points SimpleSystem::getPoints() {
	return Points(1, state_[0]);
}

Lines SimpleSystem::getLines() {
	static const auto n_lines = 50u;
	auto l = Lines(n_lines * 2);
	const auto angle_incr = 2 * FW_PI / n_lines;
	for (auto i = 0u; i < n_lines; ++i) {
		l[2 * i] = l[2 * i + 1] =
			Vec3f(radius_ * FW::sin(angle_incr * i), radius_ * FW::cos(angle_incr * i), 0);
	}
	rotate(l.begin(), l.begin() + 1, l.end());
	return l;
}

void SpringSystem::reset() {
	const auto start_pos = Vec3f(0.1f, -0.5f, 0.0f);
	
	const auto spring_k = 30.0f;
	state_ = State(4);
	// YOUR CODE HERE (R2)
	// Set the initial state for a particle system with one particle fixed
	// at origin and another particle hanging off the first one with a spring.
	// Place the second particle initially at start_pos.
	state_[0] = (0.0f, 0.0f, 0.0f);
	state_[1] = (0.0f, 0.0f, 0.0f);
	state_[2] = (start_pos);
	state_[3] = (0.0f, 0.0f, 0.0f);
	spring_.k = spring_k;
	spring_.rlen = start_pos.length();
	spring_.i1 = 0;
	spring_.i2 = 2;

}

State SpringSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	State f(4);
	// YOUR CODE HERE (R2)
	// Return a derivative for the system as if it was in state "state".
	// You can use the fGravity, fDrag and fSpring helper functions for the forces.
	f[0] = state[0];
	f[1] = state[1];
	f[2] = state[3];
	f[3] = fGravity(mass) + fSpring(state[spring_.i2], state[spring_.i1], spring_.k, spring_.rlen)/mass + fDrag(state[spring_.i2+1], drag_k);

	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

// This is a very useful read for the Jacobians of the spring forces. It deals with spring damping as well, we don't do that -- our drag is simply a linear damping of velocity (that results in some constants in the Jacobian).
// http://blog.mmacklin.com/2012/05/04/implicitsprings/

void SpringSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	// EXTRA: Evaluate the Jacobian into the 'result' matrix here. Only the free end of the spring should have any nonzero values related to it.
}
#endif

Points SpringSystem::getPoints() {
	auto p = Points(2);
	p[0] = state_[0]; p[1] = state_[2];
	return p;
}

Lines SpringSystem::getLines() {
	auto l = Lines(2);
	l[0] = state_[0]; l[1] = state_[2];
	return l;
}

void PendulumSystem::reset() {
	const auto spring_k = 1000.0f;
	const auto start_point = Vec3f(0);
	const auto end_point = Vec3f(0.05, -1.5, 0);
	state_ = State(2 * n_);
	auto mid_measure = Vec3f(0.05, -1.5, 0) / (n_ - 1);
	auto mid_length = mid_measure.length();
	//std::cout << "N: " << n_ << "\n";
	// YOUR CODE HERE (R4)
	// Set the initial state for a pendulum system with n_ particles
	state_[pos(0)] = start_point;
	springs_.clear();
	auto delta = (end_point - start_point) / (n_ - 1);
	for (int i = 1; i < n_; i++) {
		Spring s = Spring(i, i - 1, spring_k, delta.length());
		springs_.push_back(s);
		state_[pos(i)] = state_[pos(i - 1)] + delta;
	}
	//std::cout<< springs_.size() <<endl;
}

State PendulumSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 0.5f;
	auto f = State(2 * n_);
	// YOUR CODE HERE (R4)
	// As in R2, return a derivative of the system state "state".
	
	std::vector<Vec3f> sForce;
	for (int i = 0; i < n_; i++)
		sForce.push_back(0);
	for (int i = 0; i < n_ - 1; i++) {
		sForce[springs_[i].i1] += fSpring(state[pos(springs_[i].i1)], state[pos(springs_[i].i2)], springs_[i].k, springs_[i].rlen);
		sForce[springs_[i].i2] += fSpring(state[pos(springs_[i].i2)], state[pos(springs_[i].i1)], springs_[i].k, springs_[i].rlen);
	}
	for (int i = 1; i < n_; i++) {
		Vec3f force = fGravity(mass) + fDrag(state[vel(i)], drag_k) + sForce[i];
		f[pos(i)] = state[vel(i)];
		f[vel(i)] = force / mass;
	}
	f[pos(0)] = state[vel(0)];
	f[vel(0)] = Vec3f(0.0f, 0.0f, 0.0f);

	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void PendulumSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {

	const auto drag_k = 0.5f;
	const auto mass = 0.5f;

	// EXTRA: Evaluate the Jacobian here. Each spring has an effect on four blocks of the matrix -- both of the positions of the endpoints will have an effect on both of the velocities of the endpoints.
}
#endif


Points PendulumSystem::getPoints() {
	auto p = Points(n_);
	for (auto i = 0u; i < n_; ++i) {
		p[i] = state_[i * 2];
	}
	return p;
}

Lines PendulumSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(state_[2 * s.i1]);
		l.push_back(state_[2 * s.i2]);
	}
	return l;
}

unsigned ClothSystem::clothMapping(unsigned len, unsigned x, unsigned y) const {
	return len * y + x;
}

void ClothSystem::reset() {
	const auto spring_k = 300.0f;
	const auto width = 1.5f, height = 1.5f; // width and height of the whole grid
	state_ = State(2 * x_*y_);
	// YOUR CODE HERE (R5)
	// Construct a particle system with a x_ * y_ grid of particles,
	// connected with a variety of springs as described in the handout:
	// structural springs, shear springs and flex springs.
	auto delta_x = width / (x_ - 1);
	auto delta_y = height / (y_ - 1);
	auto delta_xy = FW::sqrt(delta_x * delta_x + delta_y * delta_y);
	int n = x_ * y_;
	springs_.clear();
	Spring s;
	for (int x = 0; x < x_; x++) {
		for (int y = 0; y < y_; y++) {
			state_[pos(clothMapping(x_, x, y))] = Vec3f((((float)x) / (x_ - 1)) * width - width / 2, 0, -(((float)y) / (y_ - 1)) * height);
			if (y < y_ - 1) {
				s = { clothMapping(x_, x, y), clothMapping(x_, x, y + 1), spring_k, delta_y };
				springs_.push_back(s);
			}
			if (x < x_ - 1) {
				s = { clothMapping(x_, x, y), clothMapping(x_, x + 1, y), spring_k, delta_x };
				springs_.push_back(s);
			}

			if (x < x_ - 1 && y > 0) {
				s = { clothMapping(x_, x, y), clothMapping(x_, x + 1, y - 1), spring_k, delta_xy };
				springs_.push_back(s);
			}

			if (x < x_ - 1 && y < y_ - 1) {
				s = { clothMapping(x_, x, y), clothMapping(x_, x + 1, y + 1), spring_k, delta_xy };
				springs_.push_back(s);
			}

			if (y < y_ - 2) {
				s = { clothMapping(x_, x, y), clothMapping(x_, x, y + 2), spring_k, delta_y * 2 };
				springs_.push_back(s);
			}

			if (x < x_ - 2) {
				s = { clothMapping(x_, x, y), clothMapping(x_, x + 2, y), spring_k, delta_x * 2 };
				springs_.push_back(s);
			}
		}
	}
	// std::cout << "N. of springs: " << springs_.size() << "\n";
}

State ClothSystem::evalF(const State& state) const {
	const auto drag_k = 0.08f;
	const auto n = x_ * y_;
	static const auto mass = 0.025f;
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	std::vector<Vec3f> sForce;
	Vec3f force;
	for (int i = 0; i < n; i++)
		sForce.push_back(0);

	for (int i = 0; i < springs_.size(); i++) {
		//calculationg forces taking into account both springs
		sForce[springs_[i].i1] += fSpring(state[pos(springs_[i].i1)], state[pos(springs_[i].i2)], springs_[i].k, springs_[i].rlen);
		sForce[springs_[i].i2] += fSpring(state[pos(springs_[i].i2)], state[pos(springs_[i].i1)], springs_[i].k, springs_[i].rlen);
	}
	auto gravity = fGravity(mass);
	float r2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 16));
	for (int i = 0; i < n; i++) {
		force = gravity + fDrag(state[vel(i)], drag_k) + sForce[i];
		if (winder) {
			force += Vec3f(0, 0, -8 + r2);
		}
			
		f[pos(i)] = state[vel(i)]; //derivate of position is velocity
		f[vel(i)] = force / mass; //derivate of velocity is acceleration
	}

	f[vel(clothMapping(x_, 0, 0))] = Vec3f(0);
	f[vel(clothMapping(x_, x_ - 1, 0))] = Vec3f(0);
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void ClothSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.08f;
	static const auto mass = 0.025f;

	// EXTRA: Evaluate the Jacobian here. The code is more or less the same as for the pendulum.
}

#endif

Points ClothSystem::getPoints() {
	auto n = x_ * y_;
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = state_[2 * i];
	}
	return p;
}

Lines ClothSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(state_[2 * s.i1]);
		l.push_back(state_[2 * s.i2]);
	}
	return l;
}



State FluidSystem::evalF(const State&) const {
	return State();
}

