#include <iostream>
#include <vector>

#include <SDL.h>
#include <GL/glew.h>
#include <Windows.h>
#include <gl\GL.h>
#include <gl\GLU.h>

#include "vec.h"
#include "constant.h"
#include "random.h"
#include "shader.h"
#include "rt.h"
#include "primitive.h"
#include "matrix.h"

static Shader depthShader;
static Shader blurXShader;
static Shader blurYShader;
static Shader shadeShader;

static RenderTarget *depthRT;
static RenderTarget *blurXRT;
static RenderTarget *blurYRT;

Random rnd(0);

const double SimulationSpaceWidth = 0.5; // m
const double SimulationSpaceHeight = 2.0; // m
const double SimulationSpaceDepth = 0.25; // m

struct Particle {
	Vec position_;
	Vec velocity_;
	Vec force_;
	double pressure_;
	double density_;

	std::vector<Particle*> neighbor_;

	Vec PRESSURE_FORCE;
};

void init();
int simulate();

static int fire = 0, del = 0;
const int width = 960;
const int height = 960;
double mouseX, mouseY;

int ProcessSDLEvents();

int main(int argc, char* argv[])
{    
	const SDL_VideoInfo* info = NULL;

	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		return -1;
	}

	info = SDL_GetVideoInfo( );

	if (!info) {
		return -1;
	}
	int bpp = info->vfmt->BitsPerPixel;

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

	if (SDL_SetVideoMode(width, height, bpp, SDL_OPENGL) == 0) {
		return -1;
	}
	
	if (glewInit() != GLEW_OK) {
		std::cerr << "Error: glewInit()" << std::endl;
		return -1;
	}
	std::cerr << "glewInit succeeded" << std::endl;


	init();

	while (true) {
		if (ProcessSDLEvents() < 0)
			break;


		simulate();
	}



	return 0;
}
void drawSphere(double r, int lats, int longs) {
	int i, j;
	for(i = 0; i <= lats; i++) {
		double lat0 = kPI * (-0.5 + (double) (i - 1) / lats);
		double z0  = sin(lat0);
		double zr0 =  cos(lat0);

		double lat1 = kPI * (-0.5 + (double) i / lats);
		double z1 = sin(lat1);
		double zr1 = cos(lat1);

		glBegin(GL_QUAD_STRIP);
		for(j = 0; j <= longs; j++) {
			double lng = 2 * kPI * (double) (j - 1) / longs;
			double x = cos(lng);
			double y = sin(lng);

			glNormal3f(x * zr0, y * zr0, z0);
			glVertex3f(r * x * zr0, r * y * zr0, r * z0);
			glNormal3f(x * zr1, y * zr1, z1);
			glVertex3f(r * x * zr1, r * y * zr1, r * z1);
		}
		glEnd();
	}
}


void draw_circle(const Vec &p, const double radius) {

	glPushMatrix();
	glTranslated(p.x_, p.y_, p.z_);
	drawSphere(radius, 10, 10);

	glPopMatrix();

	/*
	glPointSize(5); 
	glBegin(GL_POINTS);
	for (int i = 0; i < 1; i ++) {
		const double a = (i/16.0) * 2.0 * kPI;
		const double x = p.x_ + radius * cos(a);
		const double y = p.y_ + radius * sin(a);
		const double z = p.z_;
		glVertex3f(x, y, z);
	}
	glEnd();
	*/
}

// パーティクル初期化
/*
const int GridNumX = 12;
const int GridNumY = 72;
const int GridNumZ = 12;
const double GridWidth = 0.25; // m
const double GridHeight = 1.5; // m
const double GridDepth = 0.25; // m
*/
const int GridNumX = 16;
const int GridNumY = 64;
const int GridNumZ = 16;
const double GridWidth = 0.1; // m
const double GridHeight = 0.4; // m
const double GridDepth = 0.1; // m

const double Volume = GridWidth * GridHeight * GridDepth; // m^3
const double WaterDensity = 1000; // kg/(m^3)
const double ParticleMass = WaterDensity * Volume / (GridNumX * GridNumY * GridNumZ);
const double EffectiveRadius = 0.008;
const double ParticleRadius = 0.5 * EffectiveRadius;
const double InitDensity = WaterDensity;
const double GasStiffness = 5.0;
const double Viscosity = 1.0;
const Vec Gravity = Vec(0.0, -9.8, 0.0);
std::vector<Particle> particle_buffer[3];
int buffer_index = 0;

// カーネル定数
const double Wpoly6 = 315.0 / 64.0 / (kPI * pow(EffectiveRadius, 9));
const double GWpoly6 = -945.0 / 32.0 / (kPI * pow(EffectiveRadius, 9)); 
	
const double WSpiky = 15.0 / (kPI * pow(EffectiveRadius, 6));
const double GWSpiky = -45.0 / (kPI * pow(EffectiveRadius, 6));

const double LWVisc = 45.0 / (kPI * pow(EffectiveRadius, 6));


void init() {
	/*
	for (int i = 0; i < 3; i ++)
		particle_buffer[i].resize(GridNumX * GridNumY * GridNumZ);

	for (int ix = 0; ix < GridNumX; ++ix) {
		for (int iy = 0; iy < GridNumY; ++iy) {
			for (int iz = 0; iz < GridNumZ; ++iz) {
				particle_buffer[0][iz * GridNumX * GridNumY + iy * GridNumX + ix].position_ = 
					Vec(
					0.1 + GridWidth  * (double)(ix + 0.5) / GridNumX, 
					0.0 + GridHeight * (double)(iy + 0.5) / GridNumY, 
					0.1 + GridDepth  * (double)(iz + 0.5) / GridNumZ);
				particle_buffer[0][iz * GridNumX * GridNumY + iy * GridNumX + ix].velocity_ = Vec(0.1, -2.0, 0.0);
			}
		}
	}*/
	for (int ix = 0; ix < GridNumX; ++ix) {
		for (int iy = 0; iy < GridNumY; ++iy) {
			for (int iz = 0; iz < GridNumZ; ++iz) {

				const int cx = abs(ix - GridNumX / 2);
				const int cy = abs(iy - GridNumY / 2);
				const int cz = abs(iz - GridNumZ / 2);

				if (sqrt((double)(cx * cx + cz * cz)) >= GridNumX / 2)
					continue;

				Particle p;

				p.position_ = 
					Vec(
					0.05 + GridWidth  * (double)(ix + 0.5) / GridNumX, 
					0.0 + GridHeight * (double)(iy + 0.5) / GridNumY, 
					0.03 + GridDepth  * (double)(iz + 0.5) / GridNumZ);
				p.velocity_ = Vec(0.0, 0.0, 0.0);

				particle_buffer[0].push_back(p);
			}
		}
	}

	particle_buffer[1] = particle_buffer[0];
	particle_buffer[2] = particle_buffer[0];

	// OpenGL
	depthShader.CompileFromFile("depth.fs", "depth.vs");
	depthRT = new RenderTarget(1024, GL_CLAMP_TO_EDGE, GL_LINEAR);
	blurXShader.CompileFromFile("blurx.fs", "blur.vs");
	blurXRT = new RenderTarget(512, GL_CLAMP_TO_EDGE, GL_LINEAR);
	blurYShader.CompileFromFile("blury.fs", "blur.vs");
	blurYRT = new RenderTarget(512, GL_CLAMP_TO_EDGE, GL_LINEAR);

	shadeShader.CompileFromFile("shade.fs", "shade.vs");


	std::cout << "ParticleMass: " << ParticleMass << std::endl;

	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glEnable(GL_POINT_SPRITE);
}

const int NGridX = SimulationSpaceWidth / (EffectiveRadius * 2);
const int NGridY = SimulationSpaceHeight / (EffectiveRadius * 2);
const int NGridZ = SimulationSpaceDepth / (EffectiveRadius * 2);

std::vector<Particle*> *grid = new std::vector<Particle*>[NGridX * NGridY * NGridZ];

static double tm = 0.0;
void calc_pos(const Vec &v, int *x, int *y, int *z) {
	int ix = v.x_ / SimulationSpaceWidth * NGridX;
	int iy = v.y_ / SimulationSpaceHeight * NGridY;
	int iz = v.z_ / SimulationSpaceDepth * NGridZ;
	if (ix < 0) ix = 0;
	if (ix >= NGridX) ix = NGridX - 1;
	if (iy < 0) iy = 0;
	if (iy >= NGridY) iy = NGridY - 1;
	if (iz < 0) iz = 0;
	if (iz >= NGridZ) iz = NGridZ - 1;

	*x = ix;
	*y = iy;
	*z = iz;
}



int ProcessSDLEvents() {
	SDL_Event eve;

	SDL_PumpEvents();
	int x, y;
	SDL_GetMouseState(&x, &y);
	mouseX = (float)x / width;
	mouseY = 1.0f - (float)y / height;

	while (SDL_PollEvent(&eve)) {
		switch(eve.type) {
		case SDL_KEYDOWN: {
			int key = eve.key.keysym.sym;
			if (key == SDLK_SPACE) {
				FILE *fp = fopen("Ball.txt", "wt");
				const int particle_num = particle_buffer[0].size();

				for (int i = 0; i < particle_num; ++i) {
					Vec pos = particle_buffer[0][i].position_;
					fprintf(fp, "%f, %f, %f\n", pos.x_, pos.y_, pos.z_);
				}

				fclose(fp);
				
				/*
				const int particle_num = particle_buffer[0].size();

				for (int i = 0; i < NGridX * NGridY * NGridZ; ++i)
					grid[i].clear();
				for (int i = 0; i < particle_num; ++i) {
					int ix, iy, iz;
					calc_pos(particle_buffer[0][i].position_, &ix, &iy, &iz);
		
					grid[iz * NGridX * NGridY + iy * NGridX + ix].push_back(&particle_buffer[0][i]);
				}
				
				FILE *fp = fopen("particle.obj", "wt");
				const int Sz = 8;

				for (int ix = 0; ix < NGridX; ++ix) {
					for (int iy = 0; iy < NGridY; ++iy) {
							printf("*");
						for (int iz = 0; iz < NGridZ; ++iz) {

							for (int sx = 0; sx < Sz; ++ sx) {
								for (int sy = 0; sy < Sz; ++sy) {
									for (int sz = 0; sz < Sz; ++sz) {

										const Vec pos(
											((double)ix + ((double)sx / Sz)) * SimulationSpaceWidth / NGridX,
											((double)iy + ((double)sy / Sz)) * SimulationSpaceHeight / NGridY,
											((double)iz + ((double)sz / Sz)) * SimulationSpaceDepth / NGridZ);
										

										double minlen = kINF;

										for (int oy = -1; oy <= 1; oy ++) {
											for (int ox = -1; ox <= 1; ox ++) {
												for (int oz = -1; oz <= 1; oz ++) {
													const int nx = ox + ix;
													const int ny = oy + iy;
													const int nz = oz + iz;
													if (nx < 0 || NGridX <= nx)
														continue;
													if (ny < 0 || NGridY <= ny)
														continue;
													if (nz < 0 || NGridZ <= nz)
														continue;
				
													const int idx = nz * NGridX * NGridY + ny * NGridX + nx;

													for (int i = 0; i < grid[idx].size(); ++i) {
														const double len = (grid[idx][i]->position_ - pos).length();
														if (minlen > len)
															minlen = len;
													}
												}
											}
										}

										if (fabs(minlen - EffectiveRadius) <= 0.001)
											fprintf(fp, "v %f %f %f\n", pos.x_, pos.y_, pos.z_);


									}
								}
							}

						}
					}
				}
				fclose(fp);
				*/

				//int cell[NGridX * NGridY * NGridZ];

				/*
				FILE *fp = fopen("particle.obj", "wt");
				if (fp != NULL) {
					for (int i = 0; i < particle_buffer[0].size(); ++i) {
						fprintf(fp, "v %f %f %f\n", particle_buffer[0][i].position_.x_, particle_buffer[0][i].position_.y_, particle_buffer[0][i].position_.z_);
					}

					fclose(fp);
				}
				*/


			}
			}
			break;
			
		case SDL_MOUSEBUTTONDOWN:
			if(eve.button.button == SDL_BUTTON_LEFT){
				fire = 1;
			}
			if(eve.button.button == SDL_BUTTON_RIGHT){
				del = 1;
			}
			break;

		case SDL_QUIT:
			return -1;
		}

	}

	return 0;
}



int simulate() {

	char title[256];
	sprintf(title, "%f", tm);
	SDL_WM_SetCaption(title, NULL);

	if (fire) {
		fire = 0;
		std::vector<Particle> particle(GridNumX * GridNumY * GridNumZ);
		for (int ix = 0; ix < GridNumX; ++ix) {
			for (int iy = 0; iy < GridNumY; ++iy) {
				for (int iz = 0; iz < GridNumZ; ++iz) {
					particle[iz * GridNumX * GridNumY + iy * GridNumX + ix].position_ = 
						Vec(
						0.05 + GridWidth * (double)(ix + 0.5) / GridNumX,
						0.0 + GridHeight * (double)(iy + 0.5) / GridNumY,
						0.02 + GridDepth * (double)(iz + 0.5) / GridNumZ);
					particle[iz * GridNumX * GridNumY + iy * GridNumX + ix].velocity_ = Vec(0.0, 0.0, 0.0);
				}
			}
		}
		particle_buffer[0].insert(particle_buffer[0].end(), particle.begin(), particle.end()); 
		particle_buffer[1].insert(particle_buffer[1].end(), particle.begin(), particle.end()); 
		particle_buffer[2].insert(particle_buffer[2].end(), particle.begin(), particle.end()); 
	}


	// シミュレーション
	const double dt = 0.0005;
//	const double dt = 0.00025; // Tait

	tm += dt;

	std::vector<Particle> &new_particle = particle_buffer[buffer_index % 3];
	std::vector<Particle> &two_before_particle = particle_buffer[(buffer_index + 1) % 3];
	std::vector<Particle> &particle = particle_buffer[(buffer_index + 2) % 3];
	
	const int particle_num = new_particle.size();
	if (del) {
		del = 0;
		for (int i = 0; i < particle_num / 2; ++i) {
			particle_buffer[0].pop_back();
			particle_buffer[1].pop_back();
			particle_buffer[2].pop_back();
		}
	}
	// 近傍粒子計算
	for (int i = 0; i < NGridX * NGridY * NGridZ; ++i)
		grid[i].clear();
	for (int i = 0; i < particle_num; ++i) {
		int ix, iy, iz;
		calc_pos(particle[i].position_, &ix, &iy, &iz);
		
		grid[iz * NGridX * NGridY + iy * NGridX + ix].push_back(&particle[i]);
	}

	for (int i = 0; i < particle_num; ++i) {
		int ix0, iy0, iz0;
		calc_pos(particle[i].position_, &ix0, &iy0, &iz0);

		particle[i].neighbor_.clear();
		for (int oy = -1; oy <= 1; oy ++) {
			for (int ox = -1; ox <= 1; ox ++) {
				for (int oz = -1; oz <= 1; oz ++) {
					const int nx = ox + ix0;
					const int ny = oy + iy0;
					const int nz = oz + iz0;
					if (nx < 0 || NGridX <= nx)
						continue;
					if (ny < 0 || NGridY <= ny)
						continue;
					if (nz < 0 || NGridZ <= nz)
						continue;
				
					const int idx = nz * NGridX * NGridY + ny * NGridX + nx;
					particle[i].neighbor_.insert(particle[i].neighbor_.end(), grid[idx].begin(), grid[idx].end());
				}
			}
		}
	}

	// 各粒子位置での密度・圧力値計算
	for (int i = 0; i < particle_num; ++i) {
		particle[i].density_ = 0.0;
		// 密度計算
		const int size = particle[i].neighbor_.size();
		for (int j = 0; j < size; ++j) {
			const double r2 = (particle[i].position_ - particle[i].neighbor_[j]->position_).length_squared();
			const double q = EffectiveRadius * EffectiveRadius - r2;
			if (r2 < EffectiveRadius * EffectiveRadius) {
				particle[i].density_ += ParticleMass * Wpoly6 * q * q * q;
			}
		}

		// std::cout << particle[i].density_ << " ";
		// 圧力値計算
		particle[i].pressure_ = GasStiffness * (particle[i].density_ - InitDensity);
		
		// Tait
		/*
		const double Gamma = 7.0;
		const double Cs = 1;
		const double B = InitDensity * Cs * Cs / Gamma;
		particle[i].pressure_ = B * (pow(particle[i].density_ / InitDensity, Gamma)  - 1.0);
		*/
		
	}

	// 圧力項計算
#pragma omp parallel for schedule(dynamic, 1) num_threads(10)
	for (int i = 0; i < particle_num; ++i) {
		const int size = particle[i].neighbor_.size();
		const double pressure_i = particle[i].pressure_ / (particle[i].density_ * particle[i].density_);
		Vec pressure_force;
		Vec viscosity_force;
		for (int j = 0; j < size; ++j) {
			// 同じ粒子だったら考慮しない
			if (&particle[i] == particle[i].neighbor_[j])
				continue;

			const Vec rij = particle[i].position_ - particle[i].neighbor_[j]->position_;
			const Vec vji = particle[i].neighbor_[j]->velocity_ - particle[i].velocity_;
			const double r = rij.length();
			
			if (r > EffectiveRadius)
				continue;

			const double q = EffectiveRadius - r;
			const double pressure_j = particle[i].neighbor_[j]->pressure_ / (particle[i].neighbor_[j]->density_ * particle[i].neighbor_[j]->density_);

			Vec kernel;
			kernel = GWSpiky * q * q * rij / r;

			// 圧力項
			pressure_force = pressure_force + ParticleMass * (pressure_i + pressure_j) * kernel;
			particle[i].PRESSURE_FORCE = pressure_force;

			// 粘性
			viscosity_force = viscosity_force + ParticleMass * (vji / particle[i].neighbor_[j]->density_) * LWVisc * q;
		}

		particle[i].force_ = Vec();
		particle[i].force_  = particle[i].force_  + (-1.0 * particle[i].density_ * pressure_force);
		particle[i].force_  = particle[i].force_  + Viscosity * viscosity_force;

		// 外力
		particle[i].force_  = particle[i].force_  + Gravity * particle[i].density_;
	}


	// 位置更新
	for (int i = 0; i < particle_num; ++i) {
		new_particle[i].velocity_ = particle[i].velocity_ + dt * particle[i].force_ / particle[i].density_;
		new_particle[i].position_ = particle[i].position_ + dt * new_particle[i].velocity_;

		// シミュレーション空間との衝突判定
		/*
		double r = (new_particle[i].position_ - Vec(0.5, 0.5, 0.5)).length();
		
		Vec hit;
		if (r >= 0.7) {
			hit = normalize(new_particle[i].position_ - Vec(0.5, 0.5, 0.5));
			new_particle[i].position_ = new_particle[i].position_ -  (r - 0.7) * hit;
		}*/

		Vec hit;
		Vec repulson;
		const double R = 1000.0;
		const double DAMP = 0.0;
		if (new_particle[i].position_.x_ - ParticleRadius < 0.0) {
			repulson.x_ = R * fabs(new_particle[i].position_.x_ - ParticleRadius) + DAMP * new_particle[i].velocity_.x_;
			new_particle[i].position_.x_ = ParticleRadius;
			hit.x_ = -1.0;
		}
		if (new_particle[i].position_.x_ + ParticleRadius > SimulationSpaceWidth) {
			repulson.x_ = -R * fabs(new_particle[i].position_.x_ + ParticleRadius - SimulationSpaceWidth) - DAMP * new_particle[i].velocity_.x_;
			new_particle[i].position_.x_ = SimulationSpaceWidth - ParticleRadius;
			hit.x_ = 1.0;
		}
			
		if (new_particle[i].position_.y_ - ParticleRadius < 0.0) {
			repulson.y_ = R * fabs(new_particle[i].position_.y_ - ParticleRadius) + DAMP * new_particle[i].velocity_.y_;
			new_particle[i].position_.y_ = ParticleRadius;
			hit.y_ = -1.0;
		}
		if (new_particle[i].position_.y_ + ParticleRadius > SimulationSpaceHeight) {
			repulson.y_ = -R * fabs(new_particle[i].position_.y_ + ParticleRadius - SimulationSpaceHeight) - DAMP * new_particle[i].velocity_.y_;
			new_particle[i].position_.y_ = SimulationSpaceHeight - ParticleRadius;
			hit.y_ = 1.0;
		}

		
		if (new_particle[i].position_.z_ - ParticleRadius < 0.0) {
			repulson.z_ = R * fabs(new_particle[i].position_.z_ - ParticleRadius) + DAMP * new_particle[i].velocity_.z_;
			new_particle[i].position_.z_ = ParticleRadius;
			hit.z_ = -1.0;
		}
		if (new_particle[i].position_.z_ + ParticleRadius > SimulationSpaceDepth) {
			repulson.z_ = -R * fabs(new_particle[i].position_.z_ + ParticleRadius - SimulationSpaceDepth) - DAMP * new_particle[i].velocity_.z_;
			new_particle[i].position_.z_ = SimulationSpaceDepth - ParticleRadius;
			hit.z_ = 1.0;
		}


		if (hit.x_ == 0.0 && hit.y_ == 0.0 && hit.z_ == 0.0)
			continue;
		
		//new_particle[i].velocity_ = (hit * dt * -10 + new_particle[i].velocity_);

		//new_particle[i].velocity_ = new_particle[i].velocity_ + repulson;
		//new_particle[i].velocity_ = 0.1 * multiply(hit, new_particle[i].velocity_);
		/*
		normalize(hit);
		new_particle[i].velocity_ = new_particle[i].velocity_ - 0.2 * dot(hit, new_particle[i].velocity_) * hit;
		*/
		new_particle[i].velocity_ = new_particle[i].velocity_ - dot(hit, new_particle[i].velocity_) * hit;
	}

	buffer_index ++;
	 
	 
	// 描画
	depthRT->Bind();
	depthShader.Bind();

	const double scale = 1;
	/*
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, 1.0, 0.1, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	const double angle = mouseX * 2.0 * kPI;
	gluLookAt(2.5 * cos(angle), 1.0, 2.5 * sin(angle),
			  0.0, 0.0, 0.0,
			  0.0, 1.0, 0.0);
			  */

	
	Matrix uCamProj = Matrix::PerspectiveMatrix(60.0f, (float)width / height, 0.001f, 100.0f);
	const double angle = mouseX * 2.0 * kPI;
	Matrix uCamView = Matrix::LookAt(0.4 * sin(angle) + 0.1, 0.2, 0.4 * cos(angle) + 0.1,
					   0.1f, 0.0f, 0.1f,
					   0.0f, 1.0f, 0.0f);
	depthShader.SetUniformMatrix4x4("uCamProj", uCamProj.m);
	depthShader.SetUniformMatrix4x4("uCamView", uCamView.m);
	float uSize = 1500 * 2 * scale * ParticleRadius / 2;
	depthShader.SetUniform("uSize", uSize);

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	glBegin(GL_POINTS);

	for (int i = 0; i < particle_num; i ++) {
		// const Vec p = Vec(-scale, -scale, -scale) / 2.0 + scale * particle[i].position_;
		const Vec p = particle[i].position_;

		double pr = particle[i].pressure_ /1500.0;
		double pg = particle[i].density_ /200.0;

		glColor3f(pr + 0.4, 0.4, 0.4);
//		depthShader.SetUniformMatrix4x4("uModel", Matrix::TranslateMatrix(p.x_, p.y_, p.z_).m);
		//draw_circle(p, 2 * scale * ParticleRadius / 2);
		glVertex3f(p.x_, p.y_, p.z_);

		/*
		glColor3f(0.5, 0.0, 0.0);
		glBegin(GL_LINE_LOOP);
		glVertex3f(p.x_, p.y_, p.z_);
		glVertex3f(p.x_ + particle[i].force_.x_ * 0.0001, p.y_+ particle[i].force_.y_ * 0.0001, p.z_+ particle[i].force_.z_ * 0.0001);
		glEnd();
		glColor3f(0.0, 0.5, 0.0);
		glBegin(GL_LINE_LOOP);
		glVertex3f(p.x_, p.y_, p.z_);
		glVertex3f(p.x_ + particle[i].PRESSURE_FORCE.x_ * 0.0001, p.y_+ particle[i].PRESSURE_FORCE.y_ * 0.0001, p.z_+ particle[i].PRESSURE_FORCE.z_ * 0.0001);
		glEnd();
		*/
	}
	glEnd();


	depthShader.Unbind();
	depthRT->Unbind();
	
	// ブラー処理
	blurXRT->Bind();
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindTexture(GL_TEXTURE_2D, depthRT->Texture());
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	DrawTexture(0, 0, 1, 1);
	blurXRT->Unbind();

	for (int i = 0; i < 3; ++i) {
		blurYRT->Bind();
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		blurYShader.Bind();
		blurYShader.SetUniform("resolution", (float)blurYRT->Size(), (float)blurYRT->Size());
		blurYShader.SetUniform("texture", (int)0);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, blurXRT->Texture());
		DrawTexture(0, 0, 1, 1);
		blurYShader.Unbind();
		blurYRT->Unbind();
		
		blurXRT->Bind();
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		blurXShader.Bind();
		blurXShader.SetUniform("resolution", (float)blurXRT->Size(), (float)blurXRT->Size());
		blurXShader.SetUniform("texture", (int)0);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, blurYRT->Texture());
		DrawTexture(0, 0, 1, 1);
		blurXShader.Unbind();
		blurXRT->Unbind();
	}
	

	shadeShader.Bind();
	shadeShader.SetUniform("resolution", (float)width, (float)height);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, depthRT->Texture());
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, blurXRT->Texture());
	shadeShader.SetUniform("tex0", (int)0);
	shadeShader.SetUniform("tex1", (int)1);

	glViewport(0, 0, width, height);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	

	// テクスチャ描画
	const float aspect = (float)width / height;
	//DrawTexture(0.5f, -0.5f, 0.5f, 0.5f * aspect);
//	glBindTexture(GL_TEXTURE_2D, blurXRT->Texture());
	DrawTexture(0, 0, 1, 1);

	shadeShader.Unbind();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, 1.0, 0.1, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(2.5 * cos(angle), 1.0, 2.5 * sin(angle),
			  0.0, 0.0, 0.0,
			  0.0, 1.0, 0.0);
	
	glColor3f(1.0, 1.0, 1.0);
	const double w = scale / 2.0;
	glBegin(GL_LINE_LOOP);
		glVertex2f(-w, w);
		glVertex2f( w, w);
		glVertex2f( w, -w);
		glVertex2f( -w, -w);
	glEnd();
	
	SDL_GL_SwapBuffers();
	return 0;
}